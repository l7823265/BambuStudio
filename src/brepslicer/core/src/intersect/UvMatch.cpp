#include <intersect/UvMatch.h>
#include <intersect/BSplineFit.h>
#include <geom/GeomUtil.h>

#include <algorithm>
#include <cmath>
#include <functional>
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

// Tip–tip chord interior band for rim-support / chord-skip / rebuild filters.
// Near-tip iso on short fillet corners (ROBOT F380) often lands at t≈0.07/0.97 —
// keep the band wide so those samples still count as interior support.
constexpr double kTipChordInteriorLo = 0.005;
constexpr double kTipChordInteriorHi = 0.995;

double fval(const IFace& face, const Plane& pln, double u, double v) {
    return signedPlaneDist(pln, face.evalUV(u, v));
}

// Tip clusters inflate raw size; count well-spaced samples along the chord.
int uniqueSpanCount(const std::vector<UvPt>& span, double tol) {
    if (span.size() < 2) return static_cast<int>(span.size());
    const double chord = dist(span.front().p, span.back().p);
    const double merge = std::max({5.0 * tol, 1e-3, 0.05 * chord});
    int n = 1;
    Vec3 last = span.front().p;
    for (size_t i = 1; i < span.size(); ++i) {
        if (dist(span[i].p, last) > merge) {
            ++n;
            last = span[i].p;
        }
    }
    return n;
}

int countCoveredConstraints(const std::vector<std::vector<UvPt>>& chains,
                            const std::vector<Vec3>& boundary3d, double tol) {
    int covered = 0;
    for (const Vec3& bp : boundary3d) {
        bool hit = false;
        for (const auto& ch : chains) {
            for (const UvPt& q : ch) {
                if (dist(q.p, bp) <= tol) {
                    hit = true;
                    break;
                }
            }
            if (hit) break;
        }
        if (hit) ++covered;
    }
    return covered;
}

// Post-filter "need": enough samples to span constraints (onespan ≥6 well-spaced;
// N≥3/4: every edge∩plane tip covered by a kept chain).
bool keptChainsMeetNeed(const std::vector<std::vector<UvPt>>& chains,
                        const std::vector<Vec3>& boundary3d, double tol, int nhit) {
    if (chains.empty()) return nhit == 0 && boundary3d.size() < 2;
    if (boundary3d.size() >= 2) {
        if (countCoveredConstraints(chains, boundary3d, tol) <
            static_cast<int>(boundary3d.size()))
            return false;
        if (boundary3d.size() == 2) {
            int best = 0;
            for (const auto& ch : chains)
                best = std::max(best, uniqueSpanCount(ch, tol));
            return best >= 6;
        }
    }
    return true;
}

void uvAt(bool fix_v, double fixed, double free, double& u, double& v) {
    u = fix_v ? free : fixed;
    v = fix_v ? fixed : free;
}

// Project a 3D seed onto the face to obtain its parametric identity (u,v).
bool resolveUvIdentity(const IFace& face, const Vec3& p, double tol, double& u, double& v) {
    double iu = 0, iv = 0;
    if (!face.invertUV(p, iu, iv, tol)) return false;
    const UVBox dom = face.uvDomain();
    if (dom.periodic_u || dom.periodic_v) wrapUV(dom, iu, iv);
    u = iu;
    v = iv;
    return true;
}

double uvParamDist(double u0, double v0, double u1, double v1, const UVBox& dom) {
    double du = u0 - u1;
    double dv = v0 - v1;
    if (dom.periodic_u && dom.period_u > 1e-14) {
        du = std::fmod(du, dom.period_u);
        if (du > 0.5 * dom.period_u) du -= dom.period_u;
        if (du < -0.5 * dom.period_u) du += dom.period_u;
    }
    if (dom.periodic_v && dom.period_v > 1e-14) {
        dv = std::fmod(dv, dom.period_v);
        if (dv > 0.5 * dom.period_v) dv -= dom.period_v;
        if (dv < -0.5 * dom.period_v) dv += dom.period_v;
    }
    return std::hypot(du, dv);
}

// Seed in/out vs the face's trimmed UV domain (projected pcurves / trim wire).
// Uses the seed's identity UV from 3D→surface projection — not a 3D FaceClassifier
// probe (that misclassifies vase/fillet NURBS and lets MST bridge wrong contours).
// Edge∩plane constraint UVs support On / near-boundary membership.
bool seedInTrimUvDomain(const IFace& face, UvPt& q, double tol,
                        const std::vector<UvPt>* constraints = nullptr) {
    // Prefer the seed's known UV (iso bisection already evaluated it). Only re-project
    // when that UV is Out — invertUV dominates UVMatch cost on dense isos.
    if (inOrOn(face.classifyUV(q.u, q.v, tol))) return true;

    double u = q.u, v = q.v;
    if (resolveUvIdentity(face, q.p, tol, u, v)) {
        q.u = u;
        q.v = v;
        if (inOrOn(face.classifyUV(u, v, tol))) return true;
    } else {
        u = q.u;
        v = q.v;
    }

    // Constraint support: on/near an edge∩plane hit ⇒ treat as On the trim wire.
    if (constraints && !constraints->empty()) {
        const UVBox dom = face.uvDomain();
        const double du = std::max(1e-16, dom.umax - dom.umin);
        const double dv = std::max(1e-16, dom.vmax - dom.vmin);
        const double snap3 = std::max(tol, 1e-3);
        const double snap_uv = std::max(tol, 1e-4 * (du + dv));
        for (const UvPt& c : *constraints) {
            if (dist(q.p, c.p) <= snap3) return true;
            if (uvParamDist(u, v, c.u, c.v, dom) <= snap_uv) return true;
        }
    }
    return false;
}

// Convenience when only (u,v) is available (no 3D seed to re-identify).
bool uvOnTrimFace(const IFace& face, double u, double v, double tol) {
    return inOrOn(face.classifyUV(u, v, tol));
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
bool synthesizeConstraintUv(const Vec3& p, const std::vector<std::vector<UvPt>>& per_iso,
                            const std::vector<UvPt>& known, UvPt& out);
void mergeBoundaryIntoChain(std::vector<UvPt>& chain, const IFace& face,
                            const std::vector<Vec3>& boundary3d, double tol, double link3d);
struct SeedData;
std::vector<std::vector<UvPt>> chainConstraintSeeds(const IFace& face, const Plane& pln,
                                                      const UVBox& dom,
                                                      const std::vector<UvPt>& constraints,
                                                      const SliceOptions& opt, double link3d,
                                                      const SeedData* seeds = nullptr);
bool nearBoundaryHit(const Vec3& p, const std::vector<Vec3>& boundary, double tol);
double boundarySnapDist(double link3d, double tol);
int nearestBoundaryIndex(const Vec3& p, const std::vector<Vec3>& boundary, double maxDist);
std::vector<std::vector<UvPt>> splitChainAtTrim(const std::vector<UvPt>& chain, const IFace& face,
                                                double tol, FaceSeedStats* stats, SeedLayer* seed_out,
                                                int face_id,
                                                const std::vector<UvPt>* constraints = nullptr);
std::vector<std::vector<UvPt>> splitChainAtNeighborFaces(const std::vector<UvPt>& chain,
                                                         int self_id, const IFace& self,
                                                         const std::vector<const FaceRecord*>& neighbors,
                                                         const std::vector<Vec3>& boundary3d,
                                                         double tol, double link3d,
                                                         int& split_count);
std::vector<const FaceRecord*> neighborFaces(int self_id, int self_solid, const SliceOptions& opt,
                                             const BBox& self_box);

// Drop open-chain samples that project past the tip–tip chord. UV classify can still
// say In on surface-extension iso chatter beyond a shared tip (huapingdun F19/F38
// at z≈-11.62 / z≈11.49: spur to y≈75.61 past the shared tip, then back).
void clipChainBeyondTipChord(std::vector<UvPt>& chain, double tol, const IFace* face = nullptr,
                             std::vector<std::string>* audit = nullptr, int face_id = -1) {
    if (chain.size() < 3) return;
    const bool closed = dist(chain.front().p, chain.back().p) <= tol && chain.size() >= 4;
    if (closed) return;
    const Vec3 a = chain.front().p;
    const Vec3 b = chain.back().p;
    const Vec3 ab{b.x - a.x, b.y - a.y, b.z - a.z};
    const double ab2 = length2(ab);
    if (ab2 < tol * tol) return;
    std::vector<UvPt> out;
    out.reserve(chain.size());
    out.push_back(chain.front());
    for (size_t i = 1; i + 1 < chain.size(); ++i) {
        const Vec3& p = chain[i].p;
        const double t = dot({p.x - a.x, p.y - a.y, p.z - a.z}, ab) / ab2;
        // Small margin for tip snap / numeric chatter on the endpoint itself.
        if (t < -0.02 || t > 1.02) {
            if (audit && face) {
                const PointClass c = face->classifyUV(chain[i].u, chain[i].v, tol);
                const char* cs = c == PointClass::In ? "In" : (c == PointClass::On ? "On" : "Out");
                std::ostringstream os;
                os << "tip_chord_clip face=" << face_id << " t=" << t << " uv_class=" << cs
                   << " uv=(" << chain[i].u << "," << chain[i].v << ") p=(" << p.x << "," << p.y
                   << "," << p.z << ")";
                audit->push_back(os.str());
            }
            continue;
        }
        if (dist(out.back().p, p) <= tol) {
            out.back() = chain[i];
            continue;
        }
        out.push_back(chain[i]);
    }
    if (dist(out.back().p, chain.back().p) <= tol)
        out.back() = chain.back();
    else
        out.push_back(chain.back());
    if (out.size() >= 2) chain = std::move(out);
}

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
            if (!ch.empty()) {
                os << " start=(" << ch.front().p.x << "," << ch.front().p.y << "," << ch.front().p.z
                   << ") end=(" << ch.back().p.x << "," << ch.back().p.y << "," << ch.back().p.z
                   << ")";
            }
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
            for (size_t i = 0; i < constraint_used.size(); ++i) {
                if (constraint_used[i]) continue;
                const Vec3& p = boundary3d[i];
                os << " miss=(" << p.x << "," << p.y << "," << p.z << ")";
            }
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
        for (int it = 0; it < 24; ++it) {
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
                         double fixed, const std::vector<UvPt>& constraints, double classTol,
                         std::vector<UvPt>& out) {
    // Attach edge∩plane points to isos. Prefer isos that already have a nearby plane hit
    // (avoids copying one boundary onto every iso on large domains — huapingdun).
    // When the iso is empty, still inject on-band edge hits so narrow trim ribbons that
    // miss the uniform iso grid still get seeds near the bound edges.
    if (constraints.empty()) return;

    const double du = std::max(1e-16, dom.umax - dom.umin);
    const double dv = std::max(1e-16, dom.vmax - dom.vmin);
    const double band = std::max(classTol, 0.002 * (du + dv));
    const double near3d = std::max(1.0, 50.0 * classTol);
    for (const UvPt& c : constraints) {
        if (std::abs(signedPlaneDist(pln, c.p)) > 10.0 * classTol) continue;
        double u = c.u, v = c.v;
        if (dom.periodic_u || dom.periodic_v) wrapUV(dom, u, v);
        const double free = fix_v ? u : v;
        const double fixd = fix_v ? v : u;
        (void)free;
        if (std::abs(fixd - fixed) > band) continue;
        if (!out.empty()) {
            double best = 1e100;
            for (const UvPt& h : out) best = std::min(best, dist(h.p, c.p));
            if (best > near3d) continue;
        }
        UvPt q = c;
        q.u = u;
        q.v = v;
        if (!seedInTrimUvDomain(face, q, classTol, &constraints)) continue;
        pushHit(out, q, classTol);
    }
}

// Extra fixed-iso values at edge∩plane UVs (and small offsets) so thin trim faces are probed.
void densifyIsosNearBoundary(const UVBox& dom, bool fix_v, const std::vector<UvPt>& constraints,
                             double tol, std::vector<double>& fixed_vals) {
    if (constraints.empty()) return;
    const double du = std::max(1e-16, dom.umax - dom.umin);
    const double dv = std::max(1e-16, dom.vmax - dom.vmin);
    const double step = std::max(tol, 0.001 * (fix_v ? dv : du));
    auto pushFixed = [&](double f) {
        const double lo = fix_v ? dom.vmin : dom.umin;
        const double hi = fix_v ? dom.vmax : dom.umax;
        if (f < lo - 1e-14 || f > hi + 1e-14) return;
        f = std::max(lo, std::min(hi, f));
        for (double x : fixed_vals) {
            if (std::abs(x - f) <= 1e-12) return;
        }
        fixed_vals.push_back(f);
    };
    for (const UvPt& c : constraints) {
        double u = c.u, v = c.v;
        if (dom.periodic_u || dom.periodic_v) wrapUV(dom, u, v);
        const double f0 = fix_v ? v : u;
        pushFixed(f0);
        pushFixed(f0 - step);
        pushFixed(f0 + step);
    }
    std::sort(fixed_vals.begin(), fixed_vals.end());
}

int scoreIsoDirection(const IFace& face, const Plane& pln, bool fix_v,
                      const std::vector<double>& fixed, const std::vector<double>& probe,
                      double geomTol) {
    int singles = 0, total = 0;
    // Subsample fixed isos for direction pick — full grid is paid again in bisectIsoSeed.
    const size_t stride = fixed.size() > 16 ? 2 : 1;
    for (size_t fi = 0; fi < fixed.size(); fi += stride) {
        const double fxd = fixed[fi];
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

// Two edge∩plane hits ⇒ one open section arc. Iso tip densification makes short
// tip-local stubs; mergeBoundary's full tip–tip snap then turns each stub into an
// overlay A↔B digon and leaves mid-span iso seeds unchained (ROBOT_4 z≈197.56 F27/F258).
bool buildTwoConstraintSpan(const IFace& face, const std::vector<UvPt>& constraints,
                            const std::vector<std::vector<UvPt>>& per_iso, double tol,
                            std::vector<UvPt>& out) {
    if (constraints.size() != 2) return false;
    const UvPt& A = constraints[0];
    const UvPt& B = constraints[1];
    const Vec3 ab{B.p.x - A.p.x, B.p.y - A.p.y, B.p.z - A.p.z};
    const double ab2 = length2(ab);
    if (ab2 < tol * tol) return false;

    struct Node {
        double t;
        UvPt q;
    };
    std::vector<Node> nodes;
    nodes.reserve(per_iso.size() + 2);
    nodes.push_back({0.0, A});
    for (const auto& iso : per_iso) {
        for (UvPt h : iso) {
            // Exact tip duplicates only — constraints own the ends.
            if (dist(h.p, A.p) <= tol || dist(h.p, B.p) <= tol) continue;
            if (!seedInTrimUvDomain(face, h, tol, &constraints)) continue;
            double t = dot({h.p.x - A.p.x, h.p.y - A.p.y, h.p.z - A.p.z}, ab) / ab2;
            // Keep every In/On sample (incl. near-tip / slight past-tip iso chatter).
            // Clamp the sort key into (0,1) so past-tip In seeds still sit next to
            // their tip instead of being dropped by a t-band (ROBOT_4 F258).
            if (t < 0.0) t = 1e-9;
            else if (t > 1.0) t = 1.0 - 1e-9;
            nodes.push_back({t, std::move(h)});
        }
    }
    nodes.push_back({1.0, B});
    std::sort(nodes.begin(), nodes.end(),
              [](const Node& a, const Node& b) { return a.t < b.t; });

    out.clear();
    out.reserve(nodes.size());
    for (const Node& n : nodes) {
        if (!out.empty() && dist(out.back().p, n.q.p) <= tol) continue;
        out.push_back(n.q);
    }
    if (out.size() < 2) return false;
    out.front() = A;
    out.back() = B;
    return true;
}

void chainHits(const std::vector<std::vector<UvPt>>& per_iso, double link3d,
               std::vector<std::vector<UvPt>>& chains) {
    chains.clear();
    // Sweep along iso index. Each open chain has a tip; each hit on the next iso
    // attaches to at most one nearest tip (and each tip takes at most one hit).
    // This keeps left/right (multi-hit) branches from cross-linking.
    struct OpenChain {
        std::vector<UvPt> pts;
        size_t last_iso = 0;
    };
    std::vector<OpenChain> open;
    auto flushStale = [&](size_t iso) {
        for (auto& ch : open) {
            if (ch.pts.size() < 2) continue;
            if (iso > ch.last_iso + 2) {
                chains.push_back(std::move(ch.pts));
                ch.pts.clear();
            }
        }
        open.erase(std::remove_if(open.begin(), open.end(),
                                  [](const OpenChain& c) { return c.pts.empty(); }),
                   open.end());
    };

    for (size_t iso = 0; iso < per_iso.size(); ++iso) {
        flushStale(iso);
        const auto& hits = per_iso[iso];
        if (hits.empty()) continue;

        std::vector<char> hit_used(hits.size(), 0);
        std::vector<char> tip_used(open.size(), 0);
        struct Cand {
            double d;
            size_t tip;
            size_t hit;
        };
        std::vector<Cand> cands;
        for (size_t ti = 0; ti < open.size(); ++ti) {
            if (open[ti].pts.empty()) continue;
            if (iso > open[ti].last_iso + 2) continue;
            const Vec3& tip = open[ti].pts.back().p;
            for (size_t hi = 0; hi < hits.size(); ++hi) {
                const double d = dist(tip, hits[hi].p);
                if (d < 1e-9 || d > link3d) continue;
                cands.push_back({d, ti, hi});
            }
        }
        std::sort(cands.begin(), cands.end(),
                  [](const Cand& a, const Cand& b) { return a.d < b.d; });
        for (const Cand& c : cands) {
            if (tip_used[c.tip] || hit_used[c.hit]) continue;
            tip_used[c.tip] = 1;
            hit_used[c.hit] = 1;
            open[c.tip].pts.push_back(hits[c.hit]);
            open[c.tip].last_iso = iso;
        }
        for (size_t hi = 0; hi < hits.size(); ++hi) {
            if (hit_used[hi]) continue;
            OpenChain ch;
            ch.pts.push_back(hits[hi]);
            ch.last_iso = iso;
            open.push_back(std::move(ch));
        }
    }
    for (auto& ch : open) {
        if (ch.pts.size() >= 2) chains.push_back(std::move(ch.pts));
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

    // Edge∩plane constraints early: support seed In/On decisions via UV identity.
    std::vector<UvPt> constraints;
    for (const Vec3& p : boundary3d) {
        UvPt hit;
        if (!boundaryConstraintPt(face, p, opt.tolerance, hit)) continue;
        pushHit(constraints, hit, opt.tolerance);
    }

    std::vector<double> u_samp, v_samp;
    face.uvIsoSamples(u_samp, v_samp);
    // Adaptive iso density: faces with ≥2 edge hits need less uniform probing.
    const int iso_base = (boundary3d.size() >= 2) ? 16 : 24;
    const int probe_base = (boundary3d.size() >= 2) ? 24 : 40;
    densifySamples(u_samp, dom.umin, dom.umax, iso_base * iso_mul);
    densifySamples(v_samp, dom.vmin, dom.vmax, iso_base * iso_mul);

    const int probe_n = probe_base * probe_mul;
    const auto probe_u = freeProbe(dom.umin, dom.umax, probe_n);
    const auto probe_v = freeProbe(dom.vmin, dom.vmax, probe_n);
    const bool fix_v =
        scoreIsoDirection(face, pln, true, v_samp, probe_u, opt.geom_tolerance) >=
        scoreIsoDirection(face, pln, false, u_samp, probe_v, opt.geom_tolerance);
    std::vector<double> fixed_vals = fix_v ? v_samp : u_samp;
    densifyIsosNearBoundary(dom, fix_v, constraints, opt.tolerance, fixed_vals);
    const std::vector<double>& probe = fix_v ? probe_u : probe_v;

    std::vector<std::vector<UvPt>> per_iso;
    for (double fixed : fixed_vals) {
        std::vector<UvPt> hits;
        bisectIsoSeed(face, pln, dom, fix_v, fixed, probe, opt.geom_tolerance, opt.tolerance,
                      hits);
        addBoundaryIsoSeeds(face, pln, dom, fix_v, fixed, constraints, opt.tolerance, hits);
        // Drop extension-patch / foreign-lobe hits before chaining. Identity UV vs
        // projected trim (+ constraint support) — not a 3D FaceClassifier probe.
        hits.erase(std::remove_if(hits.begin(), hits.end(),
                                  [&](UvPt& h) {
                                      return !seedInTrimUvDomain(face, h, opt.tolerance,
                                                                &constraints);
                                  }),
                   hits.end());
        std::sort(hits.begin(), hits.end(),
                  [&](const UvPt& a, const UvPt& b) { return fix_v ? (a.u < b.u) : (a.v < b.v); });
        data.nhit += static_cast<int>(hits.size());
        per_iso.push_back(std::move(hits));
    }

    data.link3d = std::max(1e-2, 50.0 * opt.tolerance);
    if (bb.valid) {
        const double dx = bb.xmax - bb.xmin;
        const double dy = bb.ymax - bb.ymin;
        const double dz = bb.zmax - bb.zmin;
        const double diag = std::hypot(dx, std::hypot(dy, dz));
        data.link3d = std::max(data.link3d, 0.04 * diag);
        // Slab / short-cylinder faces: section chords can span the thin extent (~2–3 mm)
        // while 0.04*diag stays smaller — then constraint pairing (link3d*2) fails (R20 face 45).
        const double thin = std::min(dx, std::min(dy, dz));
        const double thick = std::max(dx, std::max(dy, dz));
        if (thin > 1e-6 && thin < 0.35 * thick)
            data.link3d = std::max(data.link3d, 1.25 * thin);
    }

    chainHits(per_iso, data.link3d, data.chains);
    data.per_iso = std::move(per_iso);

    // Recover edge∩plane tips that failed invertUV on this face (shared fillet corners
    // can Extrema-fail on one side — ROBOT_4 z≈204.58 F27 tip at (-150,-43.1)). Borrow
    // UV from the nearest in-trim iso seed so the 2-hit onespan path can run.
    if (boundary3d.size() >= 2 && constraints.size() < boundary3d.size()) {
        for (const Vec3& p : boundary3d) {
            bool have = false;
            for (const UvPt& c : constraints) {
                if (dist(c.p, p) <= std::max(opt.tolerance, 1e-3)) {
                    have = true;
                    break;
                }
            }
            if (have) continue;
            UvPt synth;
            if (!synthesizeConstraintUv(p, data.per_iso, constraints, synth)) continue;
            pushHit(constraints, synth, opt.tolerance);
            // Also place into the nearest iso so RAW seed export includes the tip
            // (recordRawSeeds only dumps per_iso points as kind=raw).
            size_t best_iso = 0;
            double best_d = 1e300;
            for (size_t ii = 0; ii < data.per_iso.size(); ++ii) {
                for (const UvPt& h : data.per_iso[ii]) {
                    const double d = dist(h.p, synth.p);
                    if (d < best_d) {
                        best_d = d;
                        best_iso = ii;
                    }
                }
            }
            if (!data.per_iso.empty()) {
                const size_t before = data.per_iso[best_iso].size();
                pushHit(data.per_iso[best_iso], synth, opt.tolerance);
                if (data.per_iso[best_iso].size() > before) ++data.nhit;
            }
        }
    }

    // Edge∩plane constraints always participate in chaining (not only when on iso grid).
    if (!constraints.empty()) {
        std::vector<UvPt> two_hit_span;
        if (constraints.size() == 2 &&
            buildTwoConstraintSpan(face, constraints, data.per_iso, opt.tolerance,
                                  two_hit_span)) {
            // Fitting a tip–tip spline needs enough samples; sparse spans leave an
            // open gap vs the straight constraint wall (ROBOT_4 F258). Densify isos
            // once inside the primary seed path — no recovery overlay.
            const int span_n = uniqueSpanCount(two_hit_span, opt.tolerance);
            if (span_n < 6 && iso_mul <= 1 && probe_mul <= 1) {
                SeedData denser =
                    seedIsoIntersections(face, pln, dom, boundary3d, opt, /*probe_mul=*/2,
                                         /*iso_mul=*/2);
                if (denser.chains.size() == 1 &&
                    uniqueSpanCount(denser.chains.front(), opt.tolerance) >= span_n &&
                    denser.chains.front().size() >= two_hit_span.size())
                    return denser;
            }
            data.chains.clear();
            data.chains.push_back(std::move(two_hit_span));
        } else if (constraints.size() >= 3) {
            // Blind-hole tops / multi-hit rims: chainHits+mergeBoundary invents short
            // tip–tip digons across the hole (ROBOT_4 z≈278 F256/F257). Prefer
            // seed-supported rim pairing (NN / perfect match) instead.
            std::vector<std::vector<UvPt>> rim =
                chainConstraintSeeds(face, pln, dom, constraints, opt, data.link3d, &data);
            if (!rim.empty()) {
                data.chains = std::move(rim);
            } else if (data.chains.empty()) {
                data.chains =
                    chainConstraintSeeds(face, pln, dom, constraints, opt, data.link3d);
                data.nhit += static_cast<int>(constraints.size());
            } else {
                for (std::vector<UvPt>& ch : data.chains)
                    mergeBoundaryIntoChain(ch, face, boundary3d, opt.tolerance, data.link3d);
            }
        } else if (data.chains.empty() && constraints.size() >= 2) {
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
    const double chord = dist(a.p, b.p);
    // UV chords on folded/periodic domains can cross other plane∩surface zeros far from
    // the tip segment (WithSphere z≈65.7 faces 9/10/21: mid sample at z≈47). Keep only
    // hits near the 3D tip chord for short/local pairs.
    const double near_chord =
        chord <= 8.0 ? std::max(1.0, chord + 0.5) : std::max(8.0, 0.35 * chord);
    const double plane_tol = std::max(10.0 * opt.geom_tolerance, opt.tolerance * 5.0);
    auto nearTipChord = [&](const Vec3& p) {
        if (std::abs(signedPlaneDist(pln, p)) > plane_tol) return false;
        const Vec3 ab = b.p - a.p;
        const double ab2 = dot(ab, ab);
        if (ab2 < 1e-18) return dist(p, a.p) <= near_chord;
        double t = dot(p - a.p, ab) / ab2;
        t = std::max(0.0, std::min(1.0, t));
        const Vec3 proj{a.p.x + t * ab.x, a.p.y + t * ab.y, a.p.z + t * ab.z};
        return dist(p, proj) <= near_chord;
    };
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
            if (seedInTrimUvDomain(face, hit, opt.tolerance, nullptr) && nearTipChord(hit.p))
                pushHit(chain, hit, opt.tolerance);
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
        if (seedInTrimUvDomain(face, hit, opt.tolerance, nullptr) && nearTipChord(hit.p))
            pushHit(chain, hit, opt.tolerance);
    }
}

std::vector<std::vector<UvPt>> chainConstraintSeeds(const IFace& face, const Plane& pln,
                                                      const UVBox& dom,
                                                      const std::vector<UvPt>& constraints,
                                                      const SliceOptions& opt, double link3d,
                                                      const SeedData* seeds) {
    std::vector<std::vector<UvPt>> chains;
    if (constraints.size() < 2) return chains;

    auto emitPair = [&](size_t i, size_t j) {
        std::vector<UvPt> chain{constraints[i]};
        const double pair_d = dist(constraints[i].p, constraints[j].p);
        const double cover = 1.25;
        const double max_jump_cap = std::min(16.0, std::max({12.0, 0.85 * pair_d}));
        // Walk / chord-collect iso seeds whenever available. Short torus corners
        // (ROBOT z≈290.6 F380 ~7 mm) used to take only UV lerp hits (pair_d≤8), which
        // exit fillet trim and leave tip–tip stubs (41 raw → 4 kept).
        if (seeds && pair_d > 3.0) {
            std::vector<UvPt> pool;
            for (const auto& iso : seeds->per_iso)
                for (const UvPt& h : iso) pool.push_back(h);
            for (const auto& ch : seeds->chains)
                for (const UvPt& h : ch) pool.push_back(h);
            const double step = std::min(14.0, std::max({link3d * 3.0, 5.0}));
            std::vector<char> used_p(pool.size(), 0);
            for (int guard = 0; guard < 96; ++guard) {
                const Vec3 tip = chain.back().p;
                if (dist(tip, constraints[j].p) <= cover) break;
                int best = -1;
                double bestD = step;
                double bestTo = 1e300;
                for (size_t pi = 0; pi < pool.size(); ++pi) {
                    if (used_p[pi]) continue;
                    const double d_tip = dist(tip, pool[pi].p);
                    if (d_tip > step || d_tip < 1e-9) continue;
                    const double d_tgt = dist(pool[pi].p, constraints[j].p);
                    if (d_tgt >= dist(tip, constraints[j].p) - 1e-9) continue;
                    if (d_tip < bestD - 1e-12 ||
                        (std::abs(d_tip - bestD) <= 1e-12 && d_tgt < bestTo)) {
                        bestD = d_tip;
                        bestTo = d_tgt;
                        best = static_cast<int>(pi);
                    }
                }
                if (best < 0) break;
                used_p[static_cast<size_t>(best)] = 1;
                chain.push_back(pool[static_cast<size_t>(best)]);
            }
            const double rem = dist(chain.back().p, constraints[j].p);
            // Refuse pairs that leap across a seed gap (ROBOT_4 z≈281 faces 256/257).
            // Allow short/medium fillet closeout (ROBOT z≈262/305 ~12 mm corners).
            // Sparse iso spacing can exceed walk step on real rim arcs (F423 ~20 mm) —
            // rebuild from chord-sorted in-band seeds instead of dropping the pair.
            // Also rebuild short corners when the walk stalls (F380 ~7 mm torus arcs).
            if (rem > max_jump_cap || (rem > cover && chain.size() < 3 && pair_d <= 16.0)) {
                const UvPt& A = constraints[i];
                const UvPt& B = constraints[j];
                const Vec3 ab{B.p.x - A.p.x, B.p.y - A.p.y, B.p.z - A.p.z};
                const double ab2 = length2(ab);
                const double chord_tol = std::max({opt.tolerance, 0.35, 0.2 * pair_d});
                struct Node {
                    double t;
                    UvPt q;
                };
                std::vector<Node> nodes;
                nodes.push_back({0.0, A});
                auto consider = [&](const UvPt& h) {
                    if (dist(h.p, A.p) <= opt.tolerance || dist(h.p, B.p) <= opt.tolerance)
                        return;
                    const double t =
                        ab2 > 1e-18
                            ? dot({h.p.x - A.p.x, h.p.y - A.p.y, h.p.z - A.p.z}, ab) / ab2
                            : 0.0;
                    if (t <= kTipChordInteriorLo || t >= kTipChordInteriorHi) return;
                    const Vec3 proj{A.p.x + t * ab.x, A.p.y + t * ab.y, A.p.z + t * ab.z};
                    if (dist(h.p, proj) > chord_tol) return;
                    nodes.push_back({t, h});
                };
                for (const UvPt& h : pool) consider(h);
                if (nodes.size() >= 2 || pair_d > 16.0) {
                    if (nodes.size() < 2) return;  // long span, A only — no in-band support
                    nodes.push_back({1.0, B});
                    std::sort(nodes.begin(), nodes.end(),
                              [](const Node& a, const Node& b) { return a.t < b.t; });
                    chain.clear();
                    for (const Node& n : nodes) {
                        if (!chain.empty() && dist(chain.back().p, n.q.p) <= opt.tolerance) continue;
                        chain.push_back(n.q);
                    }
                    if (chain.size() >= 2) chains.push_back(std::move(chain));
                    return;
                }
                chain.assign(1, constraints[i]);
                addPlaneHitsAlongUv(face, pln, dom, constraints[i], constraints[j], opt, chain);
            }
        } else {
            addPlaneHitsAlongUv(face, pln, dom, constraints[i], constraints[j], opt, chain);
        }
        chain.push_back(constraints[j]);
        if (chain.size() >= 2) chains.push_back(std::move(chain));
    };

    // Multi-hit faces: only emit seed-supported (or short in-trim UV) pairs. Greedy
    // A→B→C with pair_link from max_nn joins disconnected clusters (z≈281 3-hit
    // fillets). Blind-hole tops (z≈278 F256/F257): reject hole-crossing tip–tip
    // chords that leave trim / have no interior rim seeds.
    if (constraints.size() >= 3 && seeds) {
        auto tipPairHasRimSupport = [&](size_t i, size_t j) -> bool {
            const UvPt& A = constraints[i];
            const UvPt& B = constraints[j];
            const double span = dist(A.p, B.p);
            if (span < opt.tolerance) return false;
            const Vec3 ab{B.p.x - A.p.x, B.p.y - A.p.y, B.p.z - A.p.z};
            const double ab2 = length2(ab);
            const double chord_tol = std::max({opt.tolerance, 0.35, 0.2 * span});

            // Prefer iso/rim seeds along the 3D tip chord. Fillet UV is often non-linear, so
            // a straight UV lerp can exit trim even when the 3D arc is seed-supported
            // (ROBOT_4 z≈307.6 F178 right corner: uv_mid_out despite iso hits on the arc).
            int interior = 0;
            auto consider = [&](const UvPt& h) {
                if (dist(h.p, A.p) <= opt.tolerance || dist(h.p, B.p) <= opt.tolerance) return;
                const double t =
                    ab2 > 1e-18
                        ? dot({h.p.x - A.p.x, h.p.y - A.p.y, h.p.z - A.p.z}, ab) / ab2
                        : 0.0;
                // Match chord-skip / rebuild bands: short torus corners (ROBOT F380
                // z≈303.67 ~10 mm) often only sample iso near the tips (t≈0.07/0.97);
                // a narrow interior band false-rejects those arcs → one corner dropped.
                if (t <= kTipChordInteriorLo || t >= kTipChordInteriorHi) return;
                const Vec3 proj{A.p.x + t * ab.x, A.p.y + t * ab.y, A.p.z + t * ab.z};
                if (dist(h.p, proj) > chord_tol) return;
                ++interior;
            };
            for (const auto& iso : seeds->per_iso)
                for (const UvPt& h : iso) consider(h);
            for (const auto& ch : seeds->chains)
                for (const UvPt& h : ch) consider(h);

            // Slot / dual-arc rims (ROBOT F423 z≈277.55): a tip–tip chord that skips
            // another edge∩plane hit bridges the gap (R_gap on L_gap↔R_outer) while UV
            // mid samples stay In and one-sided rim seeds still set interior>=1. Also
            // covers huapingdun F4 groove false edges (middle tip on the long chord).
            // Prefer this over UV-exterior alone: F423's real rim arcs have interior
            // seeds yet their straight 3D chords sample UV Out.
            auto chordSkipsConstraint = [&]() -> bool {
                const double tip_tol = std::max({opt.tolerance, 1.25, 0.15 * span});
                for (size_t k = 0; k < constraints.size(); ++k) {
                    if (k == i || k == j) continue;
                    const UvPt& C = constraints[k];
                    const double t =
                        ab2 > 1e-18
                            ? dot({C.p.x - A.p.x, C.p.y - A.p.y, C.p.z - A.p.z}, ab) / ab2
                            : 0.0;
                    if (t <= kTipChordInteriorLo || t >= kTipChordInteriorHi) continue;
                    const Vec3 proj{A.p.x + t * ab.x, A.p.y + t * ab.y, A.p.z + t * ab.z};
                    if (dist(C.p, proj) <= tip_tol) return true;
                }
                return false;
            };
            if (chordSkipsConstraint()) return false;
            if (interior >= 1) return true;

            // No interior seeds: reject chords that leave trim (blind-hole diameters).
            auto chordCrossesExterior = [&]() -> bool {
                for (double t : {0.35, 0.5, 0.65}) {
                    const Vec3 p3{A.p.x + t * ab.x, A.p.y + t * ab.y, A.p.z + t * ab.z};
                    double u = 0, v = 0;
                    if (face.invertUV(p3, u, v, opt.tolerance)) {
                        UvPt m{u, v, face.evalUV(u, v)};
                        if (!seedInTrimUvDomain(face, m, opt.tolerance, &constraints))
                            return true;
                        continue;
                    }
                    // invertUV failed: fall back to UV lerp (legacy F256 path).
                    UvPt m{A.u + t * (B.u - A.u), A.v + t * (B.v - A.v), p3};
                    m.p = face.evalUV(m.u, m.v);
                    if (!seedInTrimUvDomain(face, m, opt.tolerance, &constraints)) return true;
                }
                return false;
            };
            if (chordCrossesExterior()) return false;

            // Short/medium fillet corners may lack iso samples (ROBOT F380 z≈290≈7 mm,
            // z≈303≈10 mm; emitPair already allows ~12 mm closeouts). Require in-trim UV.
            if (span > 12.0) return false;
            std::vector<UvPt> probe{A};
            addPlaneHitsAlongUv(face, pln, dom, A, B, opt, probe);
            return probe.size() >= 2;
        };

        struct Edge {
            size_t a, b;
            double d;
        };
        std::vector<Edge> edges;
        for (size_t i = 0; i < constraints.size(); ++i) {
            for (size_t j = i + 1; j < constraints.size(); ++j) {
                if (!tipPairHasRimSupport(i, j)) continue;
                // Admission is rim/skip gated above. Do not require a greedy NN seed walk
                // here: sparse iso spacing (> step) false-rejects real rim arcs (F423
                // ~20 mm gaps, rem≈28) and falls back to gap-bridging iso chains.
                edges.push_back({i, j, dist(constraints[i].p, constraints[j].p)});
            }
        }
        std::sort(edges.begin(), edges.end(),
                  [](const Edge& a, const Edge& b) { return a.d < b.d; });

        // 4-hit blind-hole / slot tops: two disjoint rim arcs (perfect matching),
        // not an MST that would still prefer a short hole-gap chord if mis-gated.
        if (constraints.size() == 4) {
            const int idx[4] = {0, 1, 2, 3};
            const int perms[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
            auto edgeOk = [&](int a, int b) {
                const size_t lo = static_cast<size_t>(std::min(a, b));
                const size_t hi = static_cast<size_t>(std::max(a, b));
                for (const Edge& e : edges)
                    if (e.a == lo && e.b == hi) return true;
                return false;
            };
            struct Pairing {
                int a0, a1, b0, b1;
                double cost;
            };
            std::vector<Pairing> opts;
            for (const auto& pr : perms) {
                const int a0 = idx[pr[0]], a1 = idx[pr[1]], b0 = idx[pr[2]], b1 = idx[pr[3]];
                if (!edgeOk(a0, a1) || !edgeOk(b0, b1)) continue;
                const double c =
                    dist(constraints[static_cast<size_t>(a0)].p,
                         constraints[static_cast<size_t>(a1)].p) +
                    dist(constraints[static_cast<size_t>(b0)].p,
                         constraints[static_cast<size_t>(b1)].p);
                opts.push_back({a0, a1, b0, b1, c});
            }
            if (!opts.empty()) {
                std::sort(opts.begin(), opts.end(),
                          [](const Pairing& a, const Pairing& b) { return a.cost < b.cost; });
                const Pairing& best = opts.front();
                emitPair(static_cast<size_t>(best.a0), static_cast<size_t>(best.a1));
                emitPair(static_cast<size_t>(best.b0), static_cast<size_t>(best.b1));
                return chains;
            }
        }

        // Kruskal forest on supported edges (3-hit fillet paths; shared tips OK).
        std::vector<int> parent(constraints.size());
        for (size_t i = 0; i < parent.size(); ++i) parent[i] = static_cast<int>(i);
        std::function<int(int)> find = [&](int x) {
            return parent[static_cast<size_t>(x)] == x
                       ? x
                       : (parent[static_cast<size_t>(x)] = find(parent[static_cast<size_t>(x)]));
        };
        for (const Edge& e : edges) {
            const int ra = find(static_cast<int>(e.a));
            const int rb = find(static_cast<int>(e.b));
            if (ra == rb) continue;
            parent[static_cast<size_t>(ra)] = rb;
            emitPair(e.a, e.b);
        }
        return chains;
    }

    // 4-hit without seeds: min-weight perfect matching (legacy path).
    if (constraints.size() == 4) {
        const int idx[4] = {0, 1, 2, 3};
        struct Pairing {
            int a0, a1, b0, b1;
            double cost;
        };
        std::vector<Pairing> opts;
        const int perms[3][4] = {{0, 1, 2, 3}, {0, 2, 1, 3}, {0, 3, 1, 2}};
        for (const auto& pr : perms) {
            const int a0 = idx[pr[0]], a1 = idx[pr[1]], b0 = idx[pr[2]], b1 = idx[pr[3]];
            const double c =
                dist(constraints[static_cast<size_t>(a0)].p, constraints[static_cast<size_t>(a1)].p) +
                dist(constraints[static_cast<size_t>(b0)].p, constraints[static_cast<size_t>(b1)].p);
            opts.push_back({a0, a1, b0, b1, c});
        }
        std::sort(opts.begin(), opts.end(),
                  [](const Pairing& a, const Pairing& b) { return a.cost < b.cost; });
        const Pairing& best = opts.front();
        emitPair(static_cast<size_t>(best.a0), static_cast<size_t>(best.a1));
        emitPair(static_cast<size_t>(best.b0), static_cast<size_t>(best.b1));
        return chains;
    }

    // Pair radius from nearest-neighbor scale so 2-hit faces (~8 mm mid-fillet) and
    // multi-hit connectors both chain, without linking opposite long sides.
    double max_nn = 0.0;
    for (size_t i = 0; i < constraints.size(); ++i) {
        double nn = 1e300;
        for (size_t j = 0; j < constraints.size(); ++j) {
            if (i == j) continue;
            nn = std::min(nn, dist(constraints[i].p, constraints[j].p));
        }
        if (nn < 1e299) max_nn = std::max(max_nn, nn);
    }
    double pair_link = std::max({link3d * 2.0, 5.0, max_nn * 1.05});

    std::vector<char> used(constraints.size(), 0);
    auto growGreedy = [&](double link) {
        for (size_t seed = 0; seed < constraints.size(); ++seed) {
            if (used[seed]) continue;
            std::vector<UvPt> chain{constraints[seed]};
            used[seed] = 1;
            bool grew = true;
            while (grew) {
                grew = false;
                int best = -1;
                double bestD = link;
                for (size_t j = 0; j < constraints.size(); ++j) {
                    if (used[j]) continue;
                    const double d = dist(chain.back().p, constraints[j].p);
                    if (d < bestD) {
                        bestD = d;
                        best = static_cast<int>(j);
                    }
                }
                if (best < 0) break;
                addPlaneHitsAlongUv(face, pln, dom, chain.back(),
                                    constraints[static_cast<size_t>(best)], opt, chain);
                chain.push_back(constraints[static_cast<size_t>(best)]);
                used[static_cast<size_t>(best)] = 1;
                grew = true;
            }
            if (chain.size() >= 2) chains.push_back(std::move(chain));
            else
                used[seed] = 0;  // singleton — free for a looser second pass
        }
    };
    growGreedy(pair_link);

    // Second pass: leftover tips (unpaired long spans) with a looser link.
    int leftover = 0;
    for (char u : used)
        if (!u) ++leftover;
    if (leftover >= 2) {
        double span = 0;
        for (size_t i = 0; i < constraints.size(); ++i) {
            if (used[i]) continue;
            for (size_t j = i + 1; j < constraints.size(); ++j) {
                if (used[j]) continue;
                span = std::max(span, dist(constraints[i].p, constraints[j].p));
            }
        }
        const double loose = std::max(pair_link, span * 1.05);
        growGreedy(loose);
    }
    return chains;
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
    out = {u, v, p};
    return seedInTrimUvDomain(face, out, tol, nullptr);
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

Vec3 segTip(const Segment& g, bool start) {
    if (g.type == SegmentType::BSpline && !g.ctrl_pts.empty())
        return start ? g.ctrl_pts.front() : g.ctrl_pts.back();
    return start ? g.start : g.end;
}

void replaceSeedOutKeptForFace(SeedLayer* seed_out, int face_id,
                               const std::vector<std::vector<UvPt>>& chains) {
    if (!seed_out) return;
    std::vector<std::vector<Vec3>> kept;
    std::vector<int> faces;
    kept.reserve(seed_out->kept_chains.size());
    faces.reserve(seed_out->kept_chain_face_ids.size());
    for (size_t i = 0; i < seed_out->kept_chains.size(); ++i) {
        const int fid =
            i < seed_out->kept_chain_face_ids.size() ? seed_out->kept_chain_face_ids[i] : -1;
        if (fid == face_id) continue;
        kept.push_back(std::move(seed_out->kept_chains[i]));
        faces.push_back(fid);
    }
    for (const auto& ch : chains) {
        if (ch.size() < 2) continue;
        std::vector<Vec3> pts;
        pts.reserve(ch.size());
        for (const UvPt& q : ch) pts.push_back(q.p);
        kept.push_back(std::move(pts));
        faces.push_back(face_id);
    }
    seed_out->kept_chains = std::move(kept);
    seed_out->kept_chain_face_ids = std::move(faces);
}

bool boundaryConstraintPt(const IFace& face, const Vec3& p, double tol, UvPt& out) {
    double u = 0, v = 0;
    // Shared fillet corners sometimes need a wider Extrema acceptance than the slice tol
    // (invertUV accepts √d ≤ 10*tol; ROBOT_4 z≈197.56 F27/F258 far tip ~8 mm).
    if (!face.invertUV(p, u, v, tol) &&
        !face.invertUV(p, u, v, std::max(tol * 20.0, 0.05)) &&
        !face.invertUV(p, u, v, std::max(tol * 200.0, 1.0)))
        return false;
    const UVBox dom = face.uvDomain();
    if (dom.periodic_u || dom.periodic_v) wrapUV(dom, u, v);
    out = {u, v, p};
    return true;
}

// When Extrema cannot identify a tip on this face (shared F27/F258 corner at z≈204.58),
// borrow UV from the nearest in-trim iso seed so onespan can still use the 3D tip.
bool synthesizeConstraintUv(const Vec3& p, const std::vector<std::vector<UvPt>>& per_iso,
                            const std::vector<UvPt>& known, UvPt& out) {
    UvPt best{};
    double best_d = 1e300;
    bool found = false;
    auto consider = [&](const UvPt& q) {
        const double d = dist(q.p, p);
        if (d >= best_d) return;
        best_d = d;
        best = q;
        found = true;
    };
    for (const auto& iso : per_iso)
        for (const UvPt& q : iso) consider(q);
    for (const UvPt& q : known) consider(q);
    if (!found) return false;
    out = {best.u, best.v, p};
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

std::vector<const FaceRecord*> neighborFaces(int self_id, int self_solid, const SliceOptions& opt,
                                             const BBox& self_box) {
    std::vector<const FaceRecord*> out;
    if (!opt.plane_faces) return out;
    out.reserve(opt.plane_faces->size());
    // Same solid only — cross-solid faces must not steal seeds (solid-first assemble).
    const double pad = 5.0;  // mm
    for (const FaceRecord* fr : *opt.plane_faces) {
        if (!fr || fr->face_id == self_id) continue;
        if (fr->solid_id != self_solid) continue;
        if (self_box.valid && fr->box.valid) {
            if (fr->box.xmax < self_box.xmin - pad || fr->box.xmin > self_box.xmax + pad ||
                fr->box.ymax < self_box.ymin - pad || fr->box.ymin > self_box.ymax + pad ||
                fr->box.zmax < self_box.zmin - pad || fr->box.zmin > self_box.zmax + pad)
                continue;
        }
        out.push_back(fr);
    }
    return out;
}

int ownerFaceId(const IFace& self, int self_id, const std::vector<const FaceRecord*>& neighbors,
                const Vec3& p, double tol) {
    const double self_score = faceOwnershipScore(self, p, tol);
    // Already clearly on self — skip neighbor Extrema.
    if (self_score <= tol) return self_id;
    int best_id = self_id;
    double best = self_score;
    for (const FaceRecord* fr : neighbors) {
        if (!fr || !fr->face) continue;
        // Cheap reject: skip faces whose bbox is far from this seed.
        if (fr->box.valid) {
            const double pad = std::max(2.0, 10.0 * tol);
            if (p.x < fr->box.xmin - pad || p.x > fr->box.xmax + pad ||
                p.y < fr->box.ymin - pad || p.y > fr->box.ymax + pad ||
                p.z < fr->box.zmin - pad || p.z > fr->box.zmax + pad)
                continue;
        }
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
                     const UvPt& q, double tol,
                     const std::vector<Vec3>* boundary3d = nullptr) {
    // Shared edge∩plane tips belong to this face even when Extrema fails here and a
    // neighbor scores better (ROBOT_4 z≈206.93 F27 tip at (-150,-45.9) stolen by F258).
    if (boundary3d && nearBoundaryHit(q.p, *boundary3d, std::max(tol, 1e-3))) return true;
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
    // Tip densify leaves samples just inside the ends; the tip constraint then sits
    // near prev/next and must not count as a hole-in-gap (F258 onespan → A↔A digon).
    const double end_eps = std::max({tol, 0.05, 0.02 * gap_span});
    double best_left = near_dist;
    double best_right = near_dist;
    const Vec3 pn{next.p.x - prev.p.x, next.p.y - prev.p.y, next.p.z - prev.p.z};
    const double pn2 = length2(pn);

    for (const Vec3& bp : boundary3d) {
        UvPt hit;
        if (!boundaryConstraintPt(face, bp, tol, hit)) continue;
        if (dist(hit.p, prev.p) <= end_eps || dist(hit.p, next.p) <= end_eps) continue;
        if (pn2 > 1e-18) {
            const double t =
                dot({hit.p.x - prev.p.x, hit.p.y - prev.p.y, hit.p.z - prev.p.z}, pn) / pn2;
            if (t <= kTipChordInteriorLo || t >= kTipChordInteriorHi) continue;
        }
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
    if (seedOwnedBySelf(self, self_id, neighbors, chain.front(), tol, &boundary3d))
        cur.push_back(chain.front());

    for (size_t i = 1; i < chain.size(); ++i) {
        const UvPt& prev = chain[i - 1];
        const UvPt& next = chain[i];
        const bool prev_self = seedOwnedBySelf(self, self_id, neighbors, prev, tol, &boundary3d);
        const bool next_self = seedOwnedBySelf(self, self_id, neighbors, next, tol, &boundary3d);

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

// Effective gap threshold for one chain: never treat the chain's own typical
// iso spacing as a hole (face 11/31 had ~1.39 mm steps vs a hard 1.25 mm cap).
double gapSplitDistForChain(const std::vector<UvPt>& chain, double link3d, double tol) {
    double base = gapSplitDist(link3d, tol);
    if (chain.size() < 3) return base;
    std::vector<double> steps;
    steps.reserve(chain.size() - 1);
    for (size_t i = 1; i < chain.size(); ++i) {
        const double d = dist(chain[i - 1].p, chain[i].p);
        if (d > tol) steps.push_back(d);
    }
    if (steps.empty()) return base;
    const size_t mid = steps.size() / 2;
    std::nth_element(steps.begin(), steps.begin() + static_cast<std::ptrdiff_t>(mid), steps.end());
    const double median = steps[mid];
    return std::max(base, 1.5 * median);
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
    if (closed) {
        // Real closed section loops cover ≥2 edge hits; tip-loop stubs cover one and
        // must still merge the uncovered far tip (ROBOT_4 z≈197.56 F27/F258).
        const double cover = std::max(tol, 1.25);
        if (countConstraintsOnChain(chain, boundary3d, cover) >= 2) return;
    }

    // Narrow trim (few edge hits): widen tip snap so interior iso tips can still attach
    // to rim constraints (ROBOT_4 cone/torus UVMatch fallbacks).
    double snap = boundarySnapDist(link3d, tol);
    double path = 0;
    for (size_t i = 1; i < chain.size(); ++i) path += dist(chain[i - 1].p, chain[i].p);
    // Tip-local stubs on multi-hit faces must not widen enough to reach the opposite
    // hole tip (~4 mm) and invent A↔B digons (ROBOT_4 z≈278 F256/F257).
    if (boundary3d.size() <= 4 && path >= std::max(1.0, 2.0 * link3d))
        snap = std::max(snap, std::min(5.0, 4.0 * link3d));
    // 2-hit fillets: widen to full tip–tip only when the chain already spans most of the
    // arc. Blind full-span snap turns tip-local digons into overlay A↔B chains (F27/F258).
    if (boundary3d.size() == 2) {
        const double span = dist(boundary3d[0], boundary3d[1]);
        if (path >= 0.35 * span) snap = std::max(snap, span * 1.05);
    }

    int used_ib = -1;
    auto ensureEnd = [&](bool at_front) {
        const UvPt& tip = at_front ? chain.front() : chain.back();
        // Prefer a constraint distinct from the other tip (avoid collapsing both ends
        // onto the same fillet corner — ROBOT_4 face 150 tip loops).
        int ib = -1;
        double best = snap;
        for (size_t i = 0; i < boundary3d.size(); ++i) {
            if (static_cast<int>(i) == used_ib) continue;
            const double d = dist(tip.p, boundary3d[i]);
            if (d <= best) {
                best = d;
                ib = static_cast<int>(i);
            }
        }
        if (ib < 0) {
            ib = nearestBoundaryIndex(tip.p, boundary3d, snap);
            if (ib < 0 || ib == used_ib) return;
        }
        UvPt hit;
        if (!boundaryConstraintPt(face, boundary3d[static_cast<size_t>(ib)], tol, hit)) {
            // Shared fillet corners can fail Extrema on one face (ROBOT_4 z≈204.58 F27
            // tip at -150). Keep the authoritative 3D constraint; borrow UV from the tip.
            hit = {tip.u, tip.v, boundary3d[static_cast<size_t>(ib)]};
        }
        used_ib = ib;
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
    std::vector<std::vector<UvPt>> parts;
    if (chain.empty()) return parts;
    // Two-hit tip–tip spans are already complete; tip densify micro-steps shrink the
    // median gap threshold and must not shatter the arc (ROBOT_4 F27/F258).
    if (boundary3d.size() == 2 && chain.size() >= 3) {
        const double cover = std::max({tol, 1.25, link3d});
        if (countConstraintsOnChain(chain, boundary3d, cover) >= 2 &&
            dist(chain.front().p, chain.back().p) > cover) {
            parts.push_back(chain);
            return parts;
        }
    }
    const double gap_dist = gapSplitDistForChain(chain, link3d, tol);
    std::vector<UvPt> cur;
    cur.push_back(chain.front());
    for (size_t i = 1; i < chain.size(); ++i) {
        const UvPt& prev = chain[i - 1];
        const UvPt& next = chain[i];
        if (dist(prev.p, next.p) > gap_dist) {
            const GapConstraints gc =
                findGapSplitConstraints(face, prev, next, boundary3d, tol, link3d);
            // Irregular fillet iso spacing can exceed gap_dist without a real hole.
            // Only split when an edge∩plane hit sits in/near the gap.
            if (!gc.has_left && !gc.has_right) {
                cur.push_back(next);
                continue;
            }
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
                                                int face_id,
                                                const std::vector<UvPt>* constraints) {
    std::vector<std::vector<UvPt>> runs;
    std::vector<UvPt> cur;
    auto flush = [&]() {
        if (cur.size() >= 2) runs.push_back(std::move(cur));
        cur.clear();
    };
    for (UvPt q : chain) {
        if (seedInTrimUvDomain(face, q, tol, constraints)) {
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
    double snap = boundarySnapDist(link3d, tol);
    if (boundary3d.size() <= 4) snap = std::max(snap, std::min(5.0, 4.0 * link3d));
    const bool closed = dist(part.front().p, part.back().p) <= tol && part.size() >= 4;
    if (!closed) {
        const bool a_on = nearBoundaryHit(part.front().p, boundary3d, snap);
        const bool b_on = nearBoundaryHit(part.back().p, boundary3d, snap);
        if (!(a_on && b_on)) {
            // Long interior run on the trimmed face: keep if at least one end is on a
            // boundary hit (mergeBoundary may not have snapped the other yet).
            double path = 0;
            for (size_t i = 1; i < part.size(); ++i) path += dist(part[i - 1].p, part[i].p);
            const bool long_run = part.size() >= 6 && path >= std::max(1.0, 5.0 * snap);
            // Tiny tip clusters (z≈281): short in-trim runs that touch one edge hit.
            const bool tiny_tip =
                (a_on || b_on) && path >= 0.05 && path < std::max(1.0, 5.0 * snap) &&
                part.size() >= 2 && boundary3d.size() >= 3;
            if (!(((a_on || b_on) && long_run) || tiny_tip)) return;
        }
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
    const std::vector<const FaceRecord*> neighbors =
        neighborFaces(iface.face_id, iface.solid_id, opt, iface.box);
    std::vector<UvPt> constraints = collectConstraintUvPts(face, boundary3d, opt.tolerance);
    // Match seedIsoIntersections: tips that fail invertUV must still count as On during
    // trim split, or onespan endpoints get dropped (ROBOT_4 z≈204.58 F27 → end_gap≈3 mm).
    if (boundary3d.size() >= 2 && constraints.size() < boundary3d.size()) {
        for (const Vec3& p : boundary3d) {
            bool have = false;
            for (const UvPt& c : constraints) {
                if (dist(c.p, p) <= std::max(opt.tolerance, 1e-3)) {
                    have = true;
                    break;
                }
            }
            if (have) continue;
            UvPt synth;
            if (!synthesizeConstraintUv(p, seeds.per_iso, constraints, synth)) continue;
            pushHit(constraints, synth, opt.tolerance);
        }
    }
    auto& chains = seeds.chains;
    for (auto& ch : chains) {
        std::vector<std::vector<UvPt>> trim_runs =
            splitChainAtTrim(ch, face, opt.tolerance, stats, opt.seed_out, iface.face_id,
                             &constraints);
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
                    // Primary path: drop iso chatter past tip–tip (UV may still say In).
                    clipChainBeyondTipChord(part, opt.tolerance, &face, opt.constraint_audit,
                                           iface.face_id);
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

RawSegment makeSeg(const FaceRecord& iface, FitResult&& fit, int chain_idx) {
    RawSegment rs;
    rs.solid_id = iface.solid_id;
    rs.shell_id = iface.shell_id;
    rs.face_id = iface.face_id;
    rs.chain_idx = chain_idx;

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
    rs.geom.face_id = iface.face_id;
    rs.geom.chain_idx = chain_idx;
    return rs;
}

}  // namespace

std::vector<RawSegment> intersectNurbsFaceWithPlaneUvMatch(const FaceRecord& iface, const Plane& pln,
                                                           const SliceFrame& frame,
                                                           const SliceOptions& opt,
                                                           const std::vector<Vec3>& boundary3d) {
    (void)frame;
    if (!iface.face) return {};
    const IFace& face = *iface.face;
    const UVBox dom = face.uvDomain();

    auto logSetup = [&](const char* stage, int nhit, size_t nchains) {
        if (!opt.constraint_audit) return;
        std::ostringstream os;
        os << "constraint_setup face=" << iface.face_id << " stage=" << stage
           << " boundary_hits=" << boundary3d.size() << " iso_raw=" << nhit
           << " kept_chains=" << nchains;
        opt.constraint_audit->push_back(os.str());
    };

    // 1. Seeding
    SeedData seeds = seedIsoIntersections(face, pln, dom, boundary3d, opt);
    recordRawSeeds(iface, seeds, opt.seed_out);
    logSetup("after_seed", seeds.nhit, 0);

    FaceSeedStats stats;
    stats.face_id = iface.face_id;
    fillSeedStats(seeds, stats);

    // 2. Filter
    std::vector<std::vector<UvPt>> chains =
        filterTrimChains(face, iface, opt, boundary3d, seeds.link3d, seeds, &stats);

    // If all seeds drop (or kept/span samples < need), densify iso and rebuild once–twice.
    // Covers 2-hit onespan (<6 well-spaced) and N≥3/4 fillet rims with uncovered tips.
    auto clearFaceSeedArtifacts = [&]() {
        if (!opt.seed_out) return;
        auto& pts = opt.seed_out->points;
        pts.erase(std::remove_if(pts.begin(), pts.end(),
                                 [&](const SeedPoint& sp) {
                                     return sp.face_id == iface.face_id;
                                 }),
                  pts.end());
        replaceSeedOutKeptForFace(opt.seed_out, iface.face_id, {});
    };
    for (int pass = 0;
         pass < 2 &&
         !keptChainsMeetNeed(chains, boundary3d, opt.tolerance, seeds.nhit);
         ++pass) {
        clearFaceSeedArtifacts();
        const int mul = 2 * (pass + 1);
        seeds = seedIsoIntersections(face, pln, dom, boundary3d, opt, mul, mul);
        recordRawSeeds(iface, seeds, opt.seed_out);
        logSetup(pass == 0 ? "after_seed_refine1" : "after_seed_refine2", seeds.nhit, 0);
        stats = FaceSeedStats{};
        stats.face_id = iface.face_id;
        fillSeedStats(seeds, stats);
        chains = filterTrimChains(face, iface, opt, boundary3d, seeds.link3d, seeds, &stats);
    }

    // Diagnose faces that lose every iso seed on the trim filter.
    if (opt.constraint_audit && chains.empty() && seeds.nhit > 0) {
        int in_trim = 0, on_trim = 0, out_trim = 0;
        int plane_in = 0;
        double best_abs_f = 1e300;
        UvPt best{};
        const int nu = 32, nv = 32;
        for (int iu = 0; iu <= nu; ++iu) {
            for (int iv = 0; iv <= nv; ++iv) {
                const double u =
                    dom.umin + (dom.umax - dom.umin) * static_cast<double>(iu) / nu;
                const double v =
                    dom.vmin + (dom.vmax - dom.vmin) * static_cast<double>(iv) / nv;
                const PointClass c = face.classifyUV(u, v, opt.tolerance);
                if (c == PointClass::In) ++in_trim;
                else if (c == PointClass::On) ++on_trim;
                else ++out_trim;
                if (!inOrOn(c)) continue;
                const double f = std::abs(fval(face, pln, u, v));
                if (f < best_abs_f) {
                    best_abs_f = f;
                    best = {u, v, face.evalUV(u, v)};
                }
                if (f <= opt.geom_tolerance) ++plane_in;
            }
        }
        const BBox bb = face.bbox();
        std::ostringstream os;
        os << "trim_empty face=" << iface.face_id << " nhit=" << seeds.nhit
           << " boundary=" << boundary3d.size() << " dom=[" << dom.umin << "," << dom.umax << "]x["
           << dom.vmin << "," << dom.vmax << "] bbox=[" << bb.xmin << "," << bb.xmax << "]x["
           << bb.ymin << "," << bb.ymax << "]x[" << bb.zmin << "," << bb.zmax
           << "] uv_grid In=" << in_trim << " On=" << on_trim << " Out=" << out_trim
           << " plane_hits_in_trim=" << plane_in << " best_|f|=" << best_abs_f << " at uv=("
           << best.u << "," << best.v << ") p=(" << best.p.x << "," << best.p.y << "," << best.p.z
           << ")";
        opt.constraint_audit->push_back(os.str());
    }

    if (chains.empty() && opt.constraint_audit) {
        logSetup("after_filter", seeds.nhit, 0);
    }

    // 3. Fitting + 4. makeSeg — emit filter-kept chains only (no pre-assemble recovery).
    replaceSeedOutKeptForFace(opt.seed_out, iface.face_id, chains);
    stats.kept = 0;
    stats.kept_chains = static_cast<int>(chains.size());
    for (const auto& ch : chains) stats.kept += static_cast<int>(ch.size());
    if (opt.constraint_audit) {
        auditFaceConstraintCoverage(iface.face_id, boundary3d, chains, {}, opt.tolerance,
                                    opt.constraint_audit);
    }
    std::vector<RawSegment> segs;
    int kept_chain_idx = 0;
    for (size_t ci = 0; ci < chains.size(); ++ci) {
        FitResult fit = fitChain(iface, std::move(chains[ci]), opt, boundary3d, seeds.link3d);
        // Consecutive index among emitted segs — assembler uses chain_idx to relate
        // pocket wall/cap locations on one face (ROBOT_7 z≈603 F1 blind hole).
        RawSegment rs = makeSeg(iface, std::move(fit), kept_chain_idx);
        if (rs.degenerate) continue;
        // Skip near-zero fit garbage left by tip stubs.
        {
            const Vec3 a = segTip(rs.geom, true);
            const Vec3 b = segTip(rs.geom, false);
            if (dist(a, b) <= std::max(opt.tolerance, 1e-4) && !rs.closed_loop) continue;
        }
        segs.push_back(std::move(rs));
        ++kept_chain_idx;
    }
    stats.segments = static_cast<int>(segs.size());

    if (opt.constraint_audit) {
        auditFaceConstraintCoverage(iface.face_id, boundary3d, {}, segs, opt.tolerance,
                                    opt.constraint_audit);
    }

    if (opt.face_seed_stats) opt.face_seed_stats->push_back(stats);
    return segs;
}

}  // namespace brepslicer
