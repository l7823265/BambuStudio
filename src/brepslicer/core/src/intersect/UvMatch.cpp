#include <intersect/UvMatch.h>
#include <intersect/UvMarch.h>
#include <intersect/BSplineFit.h>
#include <geom/GeomUtil.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <utility>
#include <vector>

namespace brepslicer {
namespace {

struct UvPt {
    double u = 0;
    double v = 0;
    Vec3 p;
};

struct FitResult {
    Segment geom;
    bool closed_loop = false;
    std::vector<Vec3> sample_pts;
    Vec3 constraint_start;
    Vec3 constraint_end;
    bool has_open_constraints = false;
};

bool inOrOn(PointClass c) { return c == PointClass::In || c == PointClass::On; }

double fval(const IFace& face, const Plane& pln, double u, double v) {
    return signedPlaneDist(pln, face.evalUV(u, v));
}

void uvAt(bool fix_v, double fixed, double free, double& u, double& v) {
    u = fix_v ? free : fixed;
    v = fix_v ? fixed : free;
}

// Trim face: UV classifier + 3D point classifier (extension patch can be In at UV but Out in 3D).
bool uvOnTrimFace(const IFace& face, double u, double v, double tol) {
    if (!inOrOn(face.classifyUV(u, v, tol))) return false;
    return inOrOn(face.classify(face.evalUV(u, v), tol));
}

void densifySamples(std::vector<double>& s, double lo, double hi, int min_count) {
    if (s.size() < 2) {
        s = {lo, hi};
    }
    std::sort(s.begin(), s.end());
    s.erase(std::unique(s.begin(), s.end(),
                        [](double a, double b) { return std::abs(a - b) < 1e-12; }),
            s.end());
    if (s.front() > lo + 1e-14) s.insert(s.begin(), lo);
    if (s.back() < hi - 1e-14) s.push_back(hi);

    while (static_cast<int>(s.size()) < min_count) {
        const size_t before = s.size();
        std::vector<double> mid;
        mid.reserve(s.size() * 2);
        mid.push_back(s.front());
        for (size_t i = 0; i + 1 < s.size(); ++i) {
            mid.push_back(0.5 * (s[i] + s[i + 1]));
            mid.push_back(s[i + 1]);
        }
        s.swap(mid);
        s.erase(std::unique(s.begin(), s.end(),
                            [](double a, double b) { return std::abs(a - b) < 1e-12; }),
                s.end());
        if (s.size() <= before) break;
    }
}

std::vector<double> freeProbe(double lo, double hi, int n) {
    std::vector<double> p;
    p.reserve(n + 1);
    for (int i = 0; i <= n; ++i) p.push_back(lo + (hi - lo) * (static_cast<double>(i) / n));
    return p;
}

void pushHit(std::vector<UvPt>& out, UvPt hit, double dedup) {
    for (const UvPt& p : out) {
        if (dist(p.p, hit.p) < dedup) return;
    }
    out.push_back(hit);
}

bool boundaryConstraintPt(const IFace& face, const Vec3& p, double tol, UvPt& out);
void mergeBoundaryIntoChain(std::vector<UvPt>& chain, const IFace& face,
                            const std::vector<Vec3>& boundary3d, double tol, double link3d);
std::vector<std::vector<UvPt>> chainConstraintSeeds(const IFace& face, const Plane& pln,
                                                      const UVBox& dom,
                                                      const std::vector<UvPt>& constraints,
                                                      const SliceOptions& opt, double link3d);
bool nearBoundaryHit(const Vec3& p, const std::vector<Vec3>& boundary, double tol);
double boundarySnapDist(double link3d, double tol);
std::vector<std::vector<UvPt>> splitChainAtTrim(const std::vector<UvPt>& chain, const IFace& face,
                                                double tol, FaceSeedStats* stats, SeedLayer* seed_out,
                                                int face_id);

int countConstraintsOnChain(const std::vector<UvPt>& chain, const std::vector<Vec3>& boundary3d,
                            double tol, std::vector<int>* hit_indices = nullptr) {
    int n = 0;
    for (size_t i = 0; i < boundary3d.size(); ++i) {
        for (const UvPt& q : chain) {
            if (dist(q.p, boundary3d[i]) <= tol) {
                ++n;
                if (hit_indices) hit_indices->push_back(static_cast<int>(i));
                break;
            }
        }
    }
    return n;
}

double minDistToConstraints(const Vec3& p, const std::vector<Vec3>& boundary3d) {
    double best = 1e300;
    for (const Vec3& b : boundary3d) best = std::min(best, dist(p, b));
    return best;
}

void auditFaceConstraintCoverage(int face_id, const std::vector<Vec3>& boundary3d,
                                 const std::vector<std::vector<UvPt>>& chains,
                                 const std::vector<RawSegment>& segs, double tol,
                                 std::vector<std::string>* logs) {
    if (!logs) return;
    std::ostringstream os;
    os << std::fixed << std::setprecision(6);
    const bool chains_only = !chains.empty();
    const bool segs_only = chains.empty() && !segs.empty();
    if (!chains_only && !segs_only) return;

    if (chains_only) {
        os << "constraint_audit face=" << face_id << " constraints=" << boundary3d.size()
           << " kept_chains=" << chains.size();
        if (boundary3d.empty()) {
            os << " (no edge∩plane hits)";
            logs->push_back(os.str());
            return;
        }
        logs->push_back(os.str());

        std::vector<char> constraint_used(boundary3d.size(), 0);
        for (size_t ci = 0; ci < chains.size(); ++ci) {
            const std::vector<UvPt>& ch = chains[ci];
            std::vector<int> hits;
            const int on_chain = countConstraintsOnChain(ch, boundary3d, tol, &hits);
            for (int idx : hits) {
                if (idx >= 0 && idx < static_cast<int>(constraint_used.size()))
                    constraint_used[static_cast<size_t>(idx)] = 1;
            }
            const double start_gap =
                ch.empty() ? 0.0 : minDistToConstraints(ch.front().p, boundary3d);
            const double end_gap = ch.empty() ? 0.0 : minDistToConstraints(ch.back().p, boundary3d);
            os.str("");
            os << "  chain[" << ci << "] pts=" << ch.size() << " on_constraints=" << on_chain
               << "/" << boundary3d.size() << " start_gap=" << start_gap << " end_gap=" << end_gap;
            if (on_chain < 2) os << " WARN<2_constraints";
            if (start_gap > tol) os << " WARN_start_off";
            if (end_gap > tol) os << " WARN_end_off";
            logs->push_back(os.str());
        }
        int missing = 0;
        for (size_t i = 0; i < constraint_used.size(); ++i) {
            if (!constraint_used[i]) ++missing;
        }
        if (missing > 0) {
            os.str("");
            os << "  uncovered_constraints=" << missing << "/" << boundary3d.size();
            logs->push_back(os.str());
        }
        return;
    }

    for (size_t si = 0; si < segs.size(); ++si) {
        const Segment& g = segs[si].geom;
        Vec3 p0 = g.start;
        Vec3 p1 = g.end;
        if (g.type == SegmentType::BSpline && !g.ctrl_pts.empty()) {
            p0 = g.ctrl_pts.front();
            p1 = g.ctrl_pts.back();
        }
        const double g0 = minDistToConstraints(p0, boundary3d);
        const double g1 = minDistToConstraints(p1, boundary3d);
        os.str("");
        os << "  seg[" << si << "] fit_start_gap=" << g0 << " fit_end_gap=" << g1;
        if (g0 > tol || g1 > tol) os << " WARN_fit_off_constraint";
        if (g0 > 0.002 || g1 > 0.002) os << " GAP>0.002";
        logs->push_back(os.str());
    }
}

// ---------------------------------------------------------------------------
// 1. Seeding — Wei et al. iso bisection on the full surface patch (no trim).
// ---------------------------------------------------------------------------

bool acceptIsoSeed(const IFace& face, const Plane& pln, const UVBox& dom, bool fix_v, double fixed,
                   double free, double geomTol, double classTol, UvPt& hit) {
    double u = fix_v ? free : fixed;
    double v = fix_v ? fixed : free;
    if (dom.periodic_u || dom.periodic_v) wrapUV(dom, u, v);
    hit = {u, v, face.evalUV(u, v)};
    return std::abs(signedPlaneDist(pln, hit.p)) <= 10.0 * geomTol + classTol;
}

void bisectIsoSeed(const IFace& face, const Plane& pln, const UVBox& dom, bool fix_v, double fixed,
                   const std::vector<double>& free_probe, double geomTol, double classTol,
                   std::vector<UvPt>& out) {
    if (free_probe.size() < 2) return;
    auto evalAt = [&](double free) {
        return fval(face, pln, fix_v ? free : fixed, fix_v ? fixed : free);
    };

    for (size_t i = 0; i + 1 < free_probe.size(); ++i) {
        const double a = free_probe[i], b = free_probe[i + 1];
        double fa = evalAt(a), fb = evalAt(b);
        if (fa * fb > 0.0 && std::abs(fa) > geomTol && std::abs(fb) > geomTol) continue;
        if (std::abs(fa) <= geomTol) {
            UvPt hit;
            if (acceptIsoSeed(face, pln, dom, fix_v, fixed, a, geomTol, classTol, hit))
                pushHit(out, hit, classTol);
            continue;
        }
        if (fa * fb > 0.0) continue;

        double lo = a, hi = b, flo = fa, fhi = fb;
        for (int it = 0; it < 40; ++it) {
            const double c = 0.5 * (lo + hi);
            const double fc = evalAt(c);
            if (std::abs(fc) <= geomTol || 0.5 * (hi - lo) <= 1e-14) {
                lo = hi = c;
                break;
            }
            if (flo * fc <= 0.0) {
                hi = c;
                fhi = fc;
            } else {
                lo = c;
                flo = fc;
            }
        }
        UvPt hit;
        if (acceptIsoSeed(face, pln, dom, fix_v, fixed, 0.5 * (lo + hi), geomTol, classTol, hit))
            pushHit(out, hit, classTol);
    }
}

void addBoundaryIsoSeeds(const IFace& face, const Plane& pln, const UVBox& dom, bool fix_v,
                         double fixed, const std::vector<Vec3>& boundary3d, double classTol,
                         std::vector<UvPt>& out) {
    const double du = std::max(1e-16, dom.umax - dom.umin);
    const double dv = std::max(1e-16, dom.vmax - dom.vmin);
    const double band = std::max(classTol, 0.01 * (du + dv));
    for (const Vec3& p : boundary3d) {
        if (std::abs(signedPlaneDist(pln, p)) > 10.0 * classTol) continue;
        double u = 0, v = 0;
        if (!face.invertUV(p, u, v, classTol)) continue;
        if (dom.periodic_u || dom.periodic_v) wrapUV(dom, u, v);
        const double free = fix_v ? u : v;
        const double fixd = fix_v ? v : u;
        if (std::abs(fixd - fixed) > band) continue;
        pushHit(out, {u, v, p}, classTol);
    }
}

int scoreIsoDirection(const IFace& face, const Plane& pln, bool fix_v,
                      const std::vector<double>& fixed, const std::vector<double>& probe,
                      double geomTol) {
    int singles = 0, total = 0;
    for (double fxd : fixed) {
        int br = 0;
        for (size_t i = 0; i + 1 < probe.size(); ++i) {
            const double f0 = fval(face, pln, fix_v ? probe[i] : fxd, fix_v ? fxd : probe[i]);
            const double f1 =
                fval(face, pln, fix_v ? probe[i + 1] : fxd, fix_v ? fxd : probe[i + 1]);
            if (f0 * f1 <= 0.0 || std::abs(f0) <= geomTol || std::abs(f1) <= geomTol) ++br;
        }
        total += br;
        if (br == 1) ++singles;
    }
    return singles * 1000 + total;
}

void chainHits(const std::vector<std::vector<UvPt>>& per_iso, double link3d,
               std::vector<std::vector<UvPt>>& chains) {
    chains.clear();
    struct Node {
        size_t iso = 0;
        size_t idx = 0;
        bool used = false;
    };
    std::vector<Node> nodes;
    for (size_t i = 0; i < per_iso.size(); ++i)
        for (size_t j = 0; j < per_iso[i].size(); ++j) nodes.push_back({i, j, false});
    auto pt = [&](const Node& n) -> const UvPt& { return per_iso[n.iso][n.idx]; };

    for (size_t seed = 0; seed < nodes.size(); ++seed) {
        if (nodes[seed].used) continue;
        std::vector<UvPt> forward;
        Node cur = nodes[seed];
        nodes[seed].used = true;
        forward.push_back(pt(cur));
        while (true) {
            int best = -1;
            double bestD = link3d;
            for (size_t k = 0; k < nodes.size(); ++k) {
                if (nodes[k].used || nodes[k].iso <= cur.iso || nodes[k].iso > cur.iso + 2)
                    continue;
                const double d = dist(pt(cur).p, pt(nodes[k]).p);
                if (d < bestD) {
                    bestD = d;
                    best = static_cast<int>(k);
                }
            }
            if (best < 0) break;
            nodes[best].used = true;
            cur = nodes[best];
            forward.push_back(pt(cur));
        }
        std::vector<UvPt> backward;
        cur = nodes[seed];
        while (true) {
            int best = -1;
            double bestD = link3d;
            for (size_t k = 0; k < nodes.size(); ++k) {
                if (nodes[k].used || nodes[k].iso >= cur.iso || cur.iso > nodes[k].iso + 2)
                    continue;
                const double d = dist(pt(cur).p, pt(nodes[k]).p);
                if (d < bestD) {
                    bestD = d;
                    best = static_cast<int>(k);
                }
            }
            if (best < 0) break;
            nodes[best].used = true;
            cur = nodes[best];
            backward.push_back(pt(cur));
        }
        std::vector<UvPt> chain;
        for (int i = static_cast<int>(backward.size()) - 1; i >= 0; --i) chain.push_back(backward[i]);
        chain.insert(chain.end(), forward.begin(), forward.end());
        if (chain.size() >= 2) chains.push_back(std::move(chain));
    }
}

struct SeedData {
    std::vector<std::vector<UvPt>> per_iso;
    std::vector<std::vector<UvPt>> chains;
    int nhit = 0;
    double link3d = 1e-2;
};

void recordRawSeeds(const FaceRecord& iface, const SeedData& seeds, SeedLayer* out) {
    if (!out) return;
    for (const auto& iso : seeds.per_iso) {
        for (const UvPt& q : iso) {
            out->points.push_back({q.p, SeedPointKind::Raw, iface.face_id});
        }
    }
    for (const auto& ch : seeds.chains) {
        std::vector<Vec3> pts;
        pts.reserve(ch.size());
        for (const UvPt& q : ch) pts.push_back(q.p);
        out->seed_chains.push_back(std::move(pts));
    }
}

SeedData seedIsoIntersections(const IFace& face, const Plane& pln, const UVBox& dom,
                              const std::vector<Vec3>& boundary3d, const SliceOptions& opt,
                              int probe_mul = 1, int iso_mul = 1) {
    SeedData data;
    const BBox bb = face.bbox();

    std::vector<double> u_samp, v_samp;
    face.uvIsoSamples(u_samp, v_samp);
    densifySamples(u_samp, dom.umin, dom.umax, 24 * iso_mul);
    densifySamples(v_samp, dom.vmin, dom.vmax, 24 * iso_mul);

    const int probe_n = 40 * probe_mul;
    const auto probe_u = freeProbe(dom.umin, dom.umax, probe_n);
    const auto probe_v = freeProbe(dom.vmin, dom.vmax, probe_n);
    const bool fix_v =
        scoreIsoDirection(face, pln, true, v_samp, probe_u, opt.geom_tolerance) >=
        scoreIsoDirection(face, pln, false, u_samp, probe_v, opt.geom_tolerance);
    const std::vector<double>& fixed_vals = fix_v ? v_samp : u_samp;
    const std::vector<double>& probe = fix_v ? probe_u : probe_v;

    std::vector<std::vector<UvPt>> per_iso;
    for (double fixed : fixed_vals) {
        std::vector<UvPt> hits;
        bisectIsoSeed(face, pln, dom, fix_v, fixed, probe, opt.geom_tolerance, opt.tolerance,
                      hits);
        addBoundaryIsoSeeds(face, pln, dom, fix_v, fixed, boundary3d, opt.tolerance, hits);
        std::sort(hits.begin(), hits.end(),
                  [&](const UvPt& a, const UvPt& b) { return fix_v ? (a.u < b.u) : (a.v < b.v); });
        data.nhit += static_cast<int>(hits.size());
        per_iso.push_back(std::move(hits));
    }

    data.link3d = std::max(1e-2, 50.0 * opt.tolerance);
    if (bb.valid) {
        const double diag =
            std::hypot(bb.xmax - bb.xmin, std::hypot(bb.ymax - bb.ymin, bb.zmax - bb.zmin));
        data.link3d = std::max(data.link3d, 0.04 * diag);
    }

    chainHits(per_iso, data.link3d, data.chains);
    data.per_iso = std::move(per_iso);

    // Edge∩plane constraints always participate in chaining (not only when on iso grid).
    std::vector<UvPt> constraints;
    for (const Vec3& p : boundary3d) {
        UvPt hit;
        if (!boundaryConstraintPt(face, p, opt.tolerance, hit)) continue;
        pushHit(constraints, hit, opt.tolerance);
    }
    if (!constraints.empty()) {
        if (data.chains.empty() && constraints.size() >= 2) {
            data.chains = chainConstraintSeeds(face, pln, dom, constraints, opt, data.link3d);
            data.nhit += static_cast<int>(constraints.size());
        } else if (!data.chains.empty()) {
            for (std::vector<UvPt>& ch : data.chains)
                mergeBoundaryIntoChain(ch, face, boundary3d, opt.tolerance, data.link3d);
        }
    }
    return data;
}

std::vector<UvPt> collectConstraintUvPts(const IFace& face, const std::vector<Vec3>& boundary3d,
                                           double tol) {
    std::vector<UvPt> pts;
    for (const Vec3& p : boundary3d) {
        UvPt hit;
        if (!boundaryConstraintPt(face, p, tol, hit)) continue;
        pushHit(pts, hit, tol);
    }
    return pts;
}

void addPlaneHitsAlongUv(const IFace& face, const Plane& pln, const UVBox& dom, const UvPt& a,
                         const UvPt& b, const SliceOptions& opt, std::vector<UvPt>& chain) {
    const int seg_steps = 16;
    auto evalT = [&](double t) {
        double u = a.u + t * (b.u - a.u);
        double v = a.v + t * (b.v - a.v);
        if (dom.periodic_u || dom.periodic_v) wrapUV(dom, u, v);
        return fval(face, pln, u, v);
    };
    for (int s = 0; s < seg_steps; ++s) {
        const double t0 = static_cast<double>(s) / seg_steps;
        const double t1 = static_cast<double>(s + 1) / seg_steps;
        double f0 = evalT(t0);
        double f1 = evalT(t1);
        if (f0 * f1 > 0.0 && std::abs(f0) > opt.geom_tolerance && std::abs(f1) > opt.geom_tolerance)
            continue;
        if (std::abs(f0) <= opt.geom_tolerance) {
            UvPt hit{a.u + t0 * (b.u - a.u), a.v + t0 * (b.v - a.v),
                     face.evalUV(a.u + t0 * (b.u - a.u), a.v + t0 * (b.v - a.v))};
            if (dom.periodic_u || dom.periodic_v) wrapUV(dom, hit.u, hit.v);
            hit.p = face.evalUV(hit.u, hit.v);
            if (uvOnTrimFace(face, hit.u, hit.v, opt.tolerance)) pushHit(chain, hit, opt.tolerance);
            continue;
        }
        if (f0 * f1 > 0.0) continue;
        double lo = t0, hi = t1, flo = f0, fhi = f1;
        for (int it = 0; it < 32; ++it) {
            const double mid = 0.5 * (lo + hi);
            const double fm = evalT(mid);
            if (std::abs(fm) <= opt.geom_tolerance || 0.5 * (hi - lo) < 1e-14) {
                lo = hi = mid;
                break;
            }
            if (flo * fm <= 0.0) {
                hi = mid;
                fhi = fm;
            } else {
                lo = mid;
                flo = fm;
            }
        }
        const double t = 0.5 * (lo + hi);
        double u = a.u + t * (b.u - a.u);
        double v = a.v + t * (b.v - a.v);
        if (dom.periodic_u || dom.periodic_v) wrapUV(dom, u, v);
        UvPt hit{u, v, face.evalUV(u, v)};
        if (uvOnTrimFace(face, hit.u, hit.v, opt.tolerance)) pushHit(chain, hit, opt.tolerance);
    }
}

std::vector<std::vector<UvPt>> chainConstraintSeeds(const IFace& face, const Plane& pln,
                                                      const UVBox& dom,
                                                      const std::vector<UvPt>& constraints,
                                                      const SliceOptions& opt, double link3d) {
    std::vector<std::vector<UvPt>> chains;
    if (constraints.size() < 2) return chains;

    std::vector<char> used(constraints.size(), 0);
    for (size_t seed = 0; seed < constraints.size(); ++seed) {
        if (used[seed]) continue;
        std::vector<UvPt> chain{constraints[seed]};
        used[seed] = 1;
        bool grew = true;
        while (grew) {
            grew = false;
            int best = -1;
            double bestD = link3d * 2.0;
            for (size_t j = 0; j < constraints.size(); ++j) {
                if (used[j]) continue;
                const double d = dist(chain.back().p, constraints[j].p);
                if (d < bestD) {
                    bestD = d;
                    best = static_cast<int>(j);
                }
            }
            if (best < 0) break;
            addPlaneHitsAlongUv(face, pln, dom, chain.back(), constraints[static_cast<size_t>(best)],
                                opt, chain);
            chain.push_back(constraints[static_cast<size_t>(best)]);
            used[static_cast<size_t>(best)] = 1;
            grew = true;
        }
        if (chain.size() >= 2) chains.push_back(std::move(chain));
    }
    return chains;
}

std::vector<std::vector<UvPt>> recoverKeptFromConstraints(
    const IFace& face, const FaceRecord& iface, const Plane& pln, const UVBox& dom,
    const std::vector<Vec3>& boundary3d, const SliceOptions& opt, double link3d,
    FaceSeedStats* stats) {
    std::vector<std::vector<UvPt>> kept;
    const std::vector<UvPt> constraints = collectConstraintUvPts(face, boundary3d, opt.tolerance);
    if (constraints.size() < 2) return kept;

    std::vector<std::vector<UvPt>> chains =
        chainConstraintSeeds(face, pln, dom, constraints, opt, link3d);
    for (std::vector<UvPt>& ch : chains) {
        std::vector<std::vector<UvPt>> trim_runs =
            splitChainAtTrim(ch, face, opt.tolerance, stats, opt.seed_out, iface.face_id);
        for (std::vector<UvPt>& run : trim_runs) {
            if (run.size() < 2) continue;
            mergeBoundaryIntoChain(run, face, boundary3d, opt.tolerance, link3d);
            const bool closed = dist(run.front().p, run.back().p) <= opt.tolerance && run.size() >= 4;
            if (!closed) {
                const double snap = boundarySnapDist(link3d, opt.tolerance);
                if (!nearBoundaryHit(run.front().p, boundary3d, snap) &&
                    !nearBoundaryHit(run.back().p, boundary3d, snap))
                    continue;
            }
            if (stats) {
                stats->kept += static_cast<int>(run.size());
                ++stats->kept_chains;
            }
            if (opt.seed_out) {
                std::vector<Vec3> pts;
                pts.reserve(run.size());
                for (const UvPt& q : run) pts.push_back(q.p);
                opt.seed_out->kept_chains.push_back(std::move(pts));
                opt.seed_out->kept_chain_face_ids.push_back(iface.face_id);
            }
            kept.push_back(std::move(run));
        }
    }
    return kept;
}

void fillSeedStats(const SeedData& seeds, FaceSeedStats& stats) {
    stats.raw = seeds.nhit;
    for (const auto& iso : seeds.per_iso) {
        if (!iso.empty()) ++stats.iso_hits;
    }
    stats.chains = static_cast<int>(seeds.chains.size());
    stats.chain_pts = 0;
    for (const auto& ch : seeds.chains) stats.chain_pts += static_cast<int>(ch.size());
}

// ---------------------------------------------------------------------------
// 2. Filter — drop seeds outside the trimmed face; split broken chains.
// ---------------------------------------------------------------------------

bool nearBoundaryHit(const Vec3& p, const std::vector<Vec3>& boundary, double tol) {
    for (const Vec3& b : boundary) {
        if (dist(p, b) <= tol) return true;
    }
    return false;
}

bool isBoundaryConstraint(const UvPt& q, const std::vector<Vec3>& boundary, double tol) {
    return nearBoundaryHit(q.p, boundary, tol);
}

bool boundaryUvPt(const IFace& face, const Vec3& p, double tol, UvPt& out) {
    double u = 0, v = 0;
    if (!face.invertUV(p, u, v, tol)) return false;
    if (!uvOnTrimFace(face, u, v, tol)) return false;
    out = {u, v, p};
    return true;
}

int nearestBoundaryIndex(const Vec3& p, const std::vector<Vec3>& boundary, double maxDist) {
    int ibest = -1;
    double best = maxDist;
    for (size_t i = 0; i < boundary.size(); ++i) {
        const double d = dist(p, boundary[i]);
        if (d <= best) {
            best = d;
            ibest = static_cast<int>(i);
        }
    }
    return ibest;
}

// Edge∩plane hit for boundary snap (allow UV on trim boundary).
bool boundaryConstraintPt(const IFace& face, const Vec3& p, double tol, UvPt& out) {
    double u = 0, v = 0;
    if (!face.invertUV(p, u, v, tol)) return false;
    const UVBox dom = face.uvDomain();
    if (dom.periodic_u || dom.periodic_v) wrapUV(dom, u, v);
    out = {u, v, p};
    return true;
}

double faceOwnershipScore(const IFace& face, const Vec3& p, double tol) {
    double u = 0, v = 0;
    if (!face.invertUV(p, u, v, tol)) return 1e100;
    if (!inOrOn(face.classifyUV(u, v, tol))) return 1e100;
    const Vec3 s = face.evalUV(u, v);
    double penalty = 0;
    if (!inOrOn(face.classify(s, tol))) penalty = 10.0;
    Vec3 Su, Sv;
    if (!face.derivUV(u, v, Su, Sv)) return penalty + dist(p, s);
    Vec3 n = cross(Su, Sv);
    if (!unitize(n, n, 1e-18)) return penalty + dist(p, s);
    return penalty + std::abs(dot({p.x - s.x, p.y - s.y, p.z - s.z}, n));
}

std::vector<const FaceRecord*> neighborFaces(int self_id, const SliceOptions& opt) {
    std::vector<const FaceRecord*> out;
    if (!opt.plane_faces) return out;
    out.reserve(opt.plane_faces->size());
    for (const FaceRecord* fr : *opt.plane_faces) {
        if (fr && fr->face_id != self_id) out.push_back(fr);
    }
    return out;
}

int ownerFaceId(const IFace& self, int self_id, const std::vector<const FaceRecord*>& neighbors,
                const Vec3& p, double tol) {
    const double self_score = faceOwnershipScore(self, p, tol);
    int best_id = self_id;
    double best = self_score;
    for (const FaceRecord* fr : neighbors) {
        if (!fr || !fr->face) continue;
        const double sc = faceOwnershipScore(*fr->face, p, tol);
        if (sc + tol < best) {
            best = sc;
            best_id = fr->face_id;
        }
    }
    // Prefer self when ownership is ambiguous (common at mirror seams, e.g. faces 38/39).
    if (best_id != self_id && self_score <= best + tol) best_id = self_id;
    return best_id;
}

UvPt interpolateUvPt(const UvPt& a, const UvPt& b, double t) {
    return {a.u + (b.u - a.u) * t, a.v + (b.v - a.v) * t,
            {a.p.x + (b.p.x - a.p.x) * t, a.p.y + (b.p.y - a.p.y) * t, a.p.z + (b.p.z - a.p.z) * t}};
}

bool seedOwnedBySelf(const IFace& self, int self_id, const std::vector<const FaceRecord*>& neighbors,
                     const UvPt& q, double tol) {
    return ownerFaceId(self, self_id, neighbors, q.p, tol) == self_id;
}

void appendUniquePt(std::vector<UvPt>& chain, const UvPt& hit, double tol) {
    if (!chain.empty() && dist(chain.back().p, hit.p) <= tol) {
        chain.back() = hit;
        return;
    }
    chain.push_back(hit);
}

double distPointSegment(const Vec3& p, const Vec3& a, const Vec3& b) {
    const Vec3 ab{b.x - a.x, b.y - a.y, b.z - a.z};
    const double L2 = length2(ab);
    double t = L2 > 0 ? dot({p.x - a.x, p.y - a.y, p.z - a.z}, ab) / L2 : 0.0;
    t = std::max(0.0, std::min(1.0, t));
    const Vec3 proj{a.x + ab.x * t, a.y + ab.y * t, a.z + ab.z * t};
    return dist(p, proj);
}

struct GapConstraints {
    bool has_left = false;
    bool has_right = false;
    UvPt left_end;
    UvPt right_start;
};

GapConstraints findGapSplitConstraints(const IFace& face, const UvPt& prev, const UvPt& next,
                                     const std::vector<Vec3>& boundary3d, double tol,
                                     double link3d) {
    GapConstraints gc;
    const double gap_span = dist(prev.p, next.p);
    const double near_dist = std::max(link3d, gap_span);
    const double seg_tol = std::max(tol, 0.15 * gap_span);
    double best_left = near_dist;
    double best_right = near_dist;

    for (const Vec3& bp : boundary3d) {
        UvPt hit;
        if (!boundaryConstraintPt(face, bp, tol, hit)) continue;
        const double ds = distPointSegment(hit.p, prev.p, next.p);
        if (ds > seg_tol && dist(hit.p, prev.p) > near_dist && dist(hit.p, next.p) > near_dist)
            continue;

        const double dl = dist(hit.p, prev.p);
        const double dr = dist(hit.p, next.p);
        if (dl <= dr) {
            if (dl < best_left) {
                best_left = dl;
                gc.left_end = hit;
                gc.has_left = true;
            }
        } else if (dr < best_right) {
            best_right = dr;
            gc.right_start = hit;
            gc.has_right = true;
        }
    }
    return gc;
}

UvPt splitConstraintAtOwnershipChange(const UvPt& a, const UvPt& b, bool a_owned,
                                      const IFace& self, int self_id,
                                      const std::vector<const FaceRecord*>& neighbors,
                                      const std::vector<Vec3>& boundary3d, double tol,
                                      double link3d) {
    const Vec3 mid{(a.p.x + b.p.x) * 0.5, (a.p.y + b.p.y) * 0.5, (a.p.z + b.p.z) * 0.5};
    const int ib = nearestBoundaryIndex(mid, boundary3d, link3d);
    if (ib >= 0) {
        UvPt hit;
        if (boundaryConstraintPt(self, boundary3d[static_cast<size_t>(ib)], tol, hit)) return hit;
    }

    UvPt lo = a;
    UvPt hi = b;
    for (int k = 0; k < 32; ++k) {
        const UvPt m = interpolateUvPt(lo, hi, 0.5);
        if (seedOwnedBySelf(self, self_id, neighbors, m, tol) == a_owned)
            lo = m;
        else
            hi = m;
    }
    UvPt out;
    if (boundaryConstraintPt(self, hi.p, tol, out)) return out;
    return hi;
}

std::vector<std::vector<UvPt>> splitChainAtNeighborFaces(const std::vector<UvPt>& chain,
                                                         int self_id, const IFace& self,
                                                         const std::vector<const FaceRecord*>& neighbors,
                                                         const std::vector<Vec3>& boundary3d,
                                                         double tol, double link3d,
                                                         int& split_count) {
    std::vector<std::vector<UvPt>> parts;
    split_count = 0;
    if (chain.size() < 2) return parts;

    std::vector<UvPt> cur;
    if (seedOwnedBySelf(self, self_id, neighbors, chain.front(), tol)) cur.push_back(chain.front());

    for (size_t i = 1; i < chain.size(); ++i) {
        const UvPt& prev = chain[i - 1];
        const UvPt& next = chain[i];
        const bool prev_self = seedOwnedBySelf(self, self_id, neighbors, prev, tol);
        const bool next_self = seedOwnedBySelf(self, self_id, neighbors, next, tol);

        if (prev_self != next_self) {
            GapConstraints gc =
                findGapSplitConstraints(self, prev, next, boundary3d, tol, link3d);
            if (!gc.has_left && !gc.has_right) {
                const UvPt split = splitConstraintAtOwnershipChange(
                    prev, next, prev_self, self, self_id, neighbors, boundary3d, tol, link3d);
                if (dist(split.p, prev.p) <= dist(split.p, next.p)) {
                    gc.left_end = split;
                    gc.has_left = true;
                } else {
                    gc.right_start = split;
                    gc.has_right = true;
                }
            }
            if (prev_self) {
                if (gc.has_left) appendUniquePt(cur, gc.left_end, tol);
                if (cur.size() >= 2) parts.push_back(std::move(cur));
                cur.clear();
            } else if (next_self && gc.has_right) {
                cur.clear();
                appendUniquePt(cur, gc.right_start, tol);
            }
        }

        if (next_self) appendUniquePt(cur, next, tol);
    }

    if (cur.size() >= 2) parts.push_back(std::move(cur));
    if (parts.size() > 1) split_count = static_cast<int>(parts.size()) - 1;
    return parts;
}

double gapSplitDist(double link3d, double tol) {
    // Hole/cap gaps are ~1.2 mm on pos2; keep below typical link3d from face bbox.
    return std::min(link3d, std::max(0.001, 1.25));
}

double boundarySnapDist(double link3d, double tol) {
    (void)tol;
    return std::max(link3d, 1.25);
}

// Insert edge∩plane boundary constraint points at open chain endpoints.
void mergeBoundaryIntoChain(std::vector<UvPt>& chain, const IFace& face,
                            const std::vector<Vec3>& boundary3d, double tol, double link3d) {
    if (chain.size() < 2 || boundary3d.empty()) return;
    const bool closed = dist(chain.front().p, chain.back().p) <= tol && chain.size() >= 4;
    if (closed) return;

    const double snap = boundarySnapDist(link3d, tol);

    auto ensureEnd = [&](bool at_front) {
        const UvPt& tip = at_front ? chain.front() : chain.back();
        const int ib = nearestBoundaryIndex(tip.p, boundary3d, snap);
        if (ib < 0) return;
        UvPt hit;
        if (!boundaryConstraintPt(face, boundary3d[static_cast<size_t>(ib)], tol, hit)) return;
        if (dist(tip.p, hit.p) <= tol) {
            if (at_front)
                chain.front() = hit;
            else
                chain.back() = hit;
            return;
        }
        if (at_front)
            chain.insert(chain.begin(), hit);
        else
            chain.push_back(hit);
    };

    ensureEnd(true);
    ensureEnd(false);
}

void snapToBoundary(UvPt& q, const std::vector<Vec3>& boundary3d, double tol) {
    const int ib = nearestBoundaryIndex(q.p, boundary3d, tol);
    if (ib >= 0) q.p = boundary3d[static_cast<size_t>(ib)];
}

std::vector<std::vector<UvPt>> splitChainAtGaps(const IFace& face, const std::vector<UvPt>& chain,
                                                double link3d, const std::vector<Vec3>& boundary3d,
                                                double tol) {
    const double gap_dist = gapSplitDist(link3d, tol);
    std::vector<std::vector<UvPt>> parts;
    if (chain.empty()) return parts;
    std::vector<UvPt> cur;
    cur.push_back(chain.front());
    for (size_t i = 1; i < chain.size(); ++i) {
        const UvPt& prev = chain[i - 1];
        const UvPt& next = chain[i];
        if (dist(prev.p, next.p) > gap_dist) {
            const GapConstraints gc =
                findGapSplitConstraints(face, prev, next, boundary3d, tol, link3d);
            if (gc.has_left) appendUniquePt(cur, gc.left_end, tol);
            if (cur.size() >= 2) parts.push_back(std::move(cur));
            cur.clear();
            if (gc.has_right) appendUniquePt(cur, gc.right_start, tol);
            appendUniquePt(cur, next, tol);
            continue;
        }
        cur.push_back(next);
    }
    if (cur.size() >= 2) parts.push_back(std::move(cur));
    return parts;
}

std::vector<std::vector<UvPt>> splitChainAtTrim(const std::vector<UvPt>& chain, const IFace& face,
                                                double tol, FaceSeedStats* stats, SeedLayer* seed_out,
                                                int face_id) {
    std::vector<std::vector<UvPt>> runs;
    std::vector<UvPt> cur;
    auto flush = [&]() {
        if (cur.size() >= 2) runs.push_back(std::move(cur));
        cur.clear();
    };
    for (const UvPt& q : chain) {
        if (uvOnTrimFace(face, q.u, q.v, tol)) {
            cur.push_back(q);
        } else {
            if (stats) ++stats->dropped;
            if (seed_out) seed_out->points.push_back({q.p, SeedPointKind::Dropped, face_id});
            flush();
        }
    }
    flush();
    return runs;
}

void keepTrimChainPart(std::vector<UvPt>&& part, const FaceRecord& iface, const SliceOptions& opt,
                       const std::vector<Vec3>& boundary3d, double link3d, double tol,
                       FaceSeedStats* stats, std::vector<std::vector<UvPt>>& kept) {
    if (part.size() < 2) return;
    const double snap = boundarySnapDist(link3d, tol);
    const bool closed = dist(part.front().p, part.back().p) <= tol && part.size() >= 4;
    if (!closed) {
        if (!nearBoundaryHit(part.front().p, boundary3d, snap) ||
            !nearBoundaryHit(part.back().p, boundary3d, snap))
            return;
    }
    if (stats) {
        stats->kept += static_cast<int>(part.size());
        ++stats->kept_chains;
    }
    if (opt.seed_out) {
        std::vector<Vec3> pts;
        pts.reserve(part.size());
        for (const UvPt& q : part) pts.push_back(q.p);
        opt.seed_out->kept_chains.push_back(std::move(pts));
        opt.seed_out->kept_chain_face_ids.push_back(iface.face_id);
    }
    kept.push_back(std::move(part));
}

std::vector<std::vector<UvPt>> filterTrimChains(const IFace& face, const FaceRecord& iface,
                                                const SliceOptions& opt,
                                                const std::vector<Vec3>& boundary3d,
                                                double link3d, SeedData& seeds,
                                                FaceSeedStats* stats = nullptr) {
    std::vector<std::vector<UvPt>> kept;
    const std::vector<const FaceRecord*> neighbors = neighborFaces(iface.face_id, opt);
    auto& chains = seeds.chains;
    for (auto& ch : chains) {
        std::vector<std::vector<UvPt>> trim_runs =
            splitChainAtTrim(ch, face, opt.tolerance, stats, opt.seed_out, iface.face_id);
        for (std::vector<UvPt>& run : trim_runs) {
            int split_count = 0;
            std::vector<std::vector<UvPt>> clipped = splitChainAtNeighborFaces(
                run, iface.face_id, face, neighbors, boundary3d, opt.tolerance, link3d,
                split_count);
            if (stats) stats->neighbor_splits += split_count;
            if (clipped.empty()) clipped.push_back(std::move(run));
            for (std::vector<UvPt>& sub : clipped) {
                for (std::vector<UvPt>& part :
                     splitChainAtGaps(face, sub, link3d, boundary3d, opt.tolerance)) {
                    mergeBoundaryIntoChain(part, face, boundary3d, opt.tolerance, link3d);
                    keepTrimChainPart(std::move(part), iface, opt, boundary3d, link3d,
                                      opt.tolerance, stats, kept);
                }
            }
        }
    }
    return kept;
}

// ---------------------------------------------------------------------------
// 3. Fitting — boundary constraint points + B-spline through interior seeds.
// ---------------------------------------------------------------------------

int nearestBoundary(const Vec3& p, const std::vector<Vec3>& boundary, double maxDist) {
    return nearestBoundaryIndex(p, boundary, maxDist);
}

bool resolveBoundaryConstraint(const Vec3& p, const std::vector<Vec3>& boundary, double tol,
                               double widen, Vec3& out) {
    int ib = nearestBoundaryIndex(p, boundary, tol);
    if (ib < 0) ib = nearestBoundaryIndex(p, boundary, widen);
    if (ib < 0) return false;
    out = boundary[static_cast<size_t>(ib)];
    return true;
}

void enforceOpenEndpointConstraints(Segment& geom, const Vec3& start, const Vec3& end) {
    geom.start = start;
    geom.end = end;
    if (geom.type == SegmentType::BSpline && !geom.ctrl_pts.empty()) {
        geom.ctrl_pts.front() = start;
        geom.ctrl_pts.back() = end;
    }
}

FitResult fitChain(const FaceRecord& iface, std::vector<UvPt> chain, const SliceOptions& opt,
                   const std::vector<Vec3>& boundary, double snapMax) {
    FitResult fit;
    if (chain.size() < 2) return fit;

    const bool closed0 = dist(chain.front().p, chain.back().p) <= opt.tolerance && chain.size() >= 4;
    if (!closed0 && !boundary.empty()) {
        Vec3 cstart, cend;
        const bool has_start =
            resolveBoundaryConstraint(chain.front().p, boundary, opt.tolerance, snapMax, cstart);
        const bool has_end =
            resolveBoundaryConstraint(chain.back().p, boundary, opt.tolerance, snapMax, cend);
        if (has_start) {
            chain.front().p = cstart;
            if (iface.face) {
                iface.face->invertUV(chain.front().p, chain.front().u, chain.front().v, opt.tolerance);
            }
            fit.constraint_start = cstart;
        }
        if (has_end) {
            chain.back().p = cend;
            if (iface.face) {
                iface.face->invertUV(chain.back().p, chain.back().u, chain.back().v, opt.tolerance);
            }
            fit.constraint_end = cend;
        }
        fit.has_open_constraints = has_start && has_end;
    }
    if (closed0 && dist(chain.front().p, chain.back().p) > opt.tolerance) {
        chain.push_back(chain.front());
    }

    fit.sample_pts.reserve(chain.size());
    for (const UvPt& q : chain) fit.sample_pts.push_back(q.p);
    const bool closed =
        dist(fit.sample_pts.front(), fit.sample_pts.back()) <= opt.tolerance &&
        fit.sample_pts.size() >= 4;
    fit.closed_loop = closed;

    fit.geom = fitCubicBSpline(fit.sample_pts, closed, opt.tolerance);
    if (fit.geom.ctrl_pts.size() < 4) {
        fit.geom = cubicFromPolyline(fit.sample_pts, closed);
    }
    if (!closed && fit.has_open_constraints) {
        enforceOpenEndpointConstraints(fit.geom, fit.constraint_start, fit.constraint_end);
    } else {
        fit.geom.start = fit.sample_pts.front();
        fit.geom.end = fit.sample_pts.back();
    }
    return fit;
}

// ---------------------------------------------------------------------------
// 4. makeSeg — package geometry into RawSegment.
// ---------------------------------------------------------------------------

RawSegment makeSeg(const FaceRecord& iface, FitResult&& fit) {
    RawSegment rs;
    rs.solid_id = iface.solid_id;
    rs.shell_id = iface.shell_id;
    rs.face_id = iface.face_id;

    if (fit.sample_pts.size() < 2) {
        rs.degenerate = true;
        return rs;
    }

    if (!fit.closed_loop && fit.has_open_constraints) {
        fit.sample_pts.front() = fit.constraint_start;
        fit.sample_pts.back() = fit.constraint_end;
    }

    if (fit.geom.ctrl_pts.size() < 4 && fit.geom.type != SegmentType::Line) {
        fit.geom = cubicFromPolyline(fit.sample_pts, fit.closed_loop);
    }
    if (fit.geom.ctrl_pts.size() < 4 && fit.geom.type != SegmentType::Line) {
        fit.geom.type = SegmentType::Line;
        fit.geom.ctrl_pts.clear();
    }

    if (!fit.closed_loop && fit.has_open_constraints) {
        enforceOpenEndpointConstraints(fit.geom, fit.constraint_start, fit.constraint_end);
    } else {
        fit.geom.start = fit.sample_pts.front();
        fit.geom.end = fit.sample_pts.back();
        if (fit.geom.type == SegmentType::BSpline && !fit.geom.ctrl_pts.empty()) {
            fit.geom.ctrl_pts.front() = fit.geom.start;
            fit.geom.ctrl_pts.back() = fit.geom.end;
        }
    }

    rs.closed_loop = fit.closed_loop;
    rs.geom = std::move(fit.geom);
    return rs;
}

}  // namespace

std::vector<RawSegment> intersectNurbsFaceWithPlaneUvMatch(const FaceRecord& iface, const Plane& pln,
                                                           const SliceFrame& frame,
                                                           const SliceOptions& opt,
                                                           const std::vector<Vec3>& boundary3d) {
    if (!iface.face) return {};
    const IFace& face = *iface.face;
    const UVBox dom = face.uvDomain();

    auto logSetup = [&](const char* stage, int nhit, size_t nchains, bool recovered) {
        if (!opt.constraint_audit) return;
        std::ostringstream os;
        os << "constraint_setup face=" << iface.face_id << " stage=" << stage
           << " boundary_hits=" << boundary3d.size() << " iso_raw=" << nhit
           << " kept_chains=" << nchains;
        if (boundary3d.size() < 2)
            os << " recovery=skipped(need>=2_boundary_hits)";
        else if (recovered)
            os << " recovery=constraint_reseed";
        else if (nchains == 0 && nhit == 0)
            os << " recovery=boost_iso_pending";
        opt.constraint_audit->push_back(os.str());
    };

    // 1. Seeding
    SeedData seeds = seedIsoIntersections(face, pln, dom, boundary3d, opt);
    recordRawSeeds(iface, seeds, opt.seed_out);
    logSetup("after_seed", seeds.nhit, 0, false);

    FaceSeedStats stats;
    stats.face_id = iface.face_id;
    fillSeedStats(seeds, stats);

    // 2. Filter
    std::vector<std::vector<UvPt>> chains =
        filterTrimChains(face, iface, opt, boundary3d, seeds.link3d, seeds, &stats);

    // Faces that hit the plane but lost all seeds after trim: constraint re-seed or denser iso.
    bool constraint_recovered = false;
    if (chains.empty() && boundary3d.size() >= 2) {
        if (seeds.nhit == 0) {
            seeds = seedIsoIntersections(face, pln, dom, boundary3d, opt, 2, 2);
            fillSeedStats(seeds, stats);
            chains = filterTrimChains(face, iface, opt, boundary3d, seeds.link3d, seeds, &stats);
            logSetup("after_boost_iso", seeds.nhit, chains.size(), false);
        }
        if (chains.empty()) {
            auto recovered = recoverKeptFromConstraints(face, iface, pln, dom, boundary3d, opt,
                                                        seeds.link3d, &stats);
            if (!recovered.empty()) constraint_recovered = true;
            chains.insert(chains.end(), std::make_move_iterator(recovered.begin()),
                          std::make_move_iterator(recovered.end()));
            logSetup("after_constraint_reseed", seeds.nhit, chains.size(), constraint_recovered);
        }
    } else if (chains.empty() && opt.constraint_audit) {
        logSetup("after_filter", seeds.nhit, 0, false);
    }

    // 3. Fitting + 4. makeSeg — kept chains must each yield a segment.
    if (opt.constraint_audit) {
        auditFaceConstraintCoverage(iface.face_id, boundary3d, chains, {}, opt.tolerance,
                                    opt.constraint_audit);
    }
    std::vector<RawSegment> segs;
    for (auto& ch : chains) {
        FitResult fit = fitChain(iface, std::move(ch), opt, boundary3d, seeds.link3d);
        RawSegment rs = makeSeg(iface, std::move(fit));
        if (rs.degenerate) continue;
        segs.push_back(std::move(rs));
    }
    stats.segments = static_cast<int>(segs.size());

    if (opt.constraint_audit) {
        auditFaceConstraintCoverage(iface.face_id, boundary3d, {}, segs, opt.tolerance,
                                    opt.constraint_audit);
    }

    if (chains.empty() && seeds.nhit == 0 && boundary3d.size() < 2) {
        stats.uvmarch_fallback = true;
        if (opt.constraint_audit) {
            std::ostringstream os;
            os << "constraint_setup face=" << iface.face_id
               << " stage=uvmarch_fallback reason=no_boundary_hits_for_recovery";
            opt.constraint_audit->push_back(os.str());
        }
        if (opt.face_seed_stats) opt.face_seed_stats->push_back(stats);
        return intersectNurbsFaceWithPlane(iface, pln, frame, opt, boundary3d);
    }
    if (segs.empty()) {
        stats.uvmarch_fallback = true;
        if (opt.face_seed_stats) opt.face_seed_stats->push_back(stats);
        return intersectNurbsFaceWithPlane(iface, pln, frame, opt, boundary3d);
    }
    if (opt.face_seed_stats) opt.face_seed_stats->push_back(stats);
    return segs;
}

}  // namespace brepslicer
