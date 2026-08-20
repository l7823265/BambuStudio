#include <intersect/UvMarch.h>
#include <intersect/BSplineFit.h>
#include <geom/GeomUtil.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace brepslicer {
namespace {

struct UvPt {
    double u = 0;
    double v = 0;
    Vec3 p;
};

struct Seed {
    double u = 0;
    double v = 0;
    bool from_boundary = false;
    bool used = false;
};

bool inOrOn(PointClass c) { return c == PointClass::In || c == PointClass::On; }

double fval(const IFace& face, const Plane& pln, double u, double v) {
    return signedPlaneDist(pln, face.evalUV(u, v));
}

bool clampDomain(const UVBox& d, double& u, double& v, bool& hit_border) {
    hit_border = false;
    if (d.periodic_u && d.period_u > 0) {
        wrapUV(d, u, v);
    } else {
        if (u < d.umin) {
            u = d.umin;
            hit_border = true;
        }
        if (u > d.umax) {
            u = d.umax;
            hit_border = true;
        }
    }
    if (d.periodic_v && d.period_v > 0) {
        wrapUV(d, u, v);
    } else {
        if (v < d.vmin) {
            v = d.vmin;
            hit_border = true;
        }
        if (v > d.vmax) {
            v = d.vmax;
            hit_border = true;
        }
    }
    return !hit_border;
}

bool newton(const IFace& face, const Plane& pln, const UVBox& dom, double& u, double& v,
            double geomTol, bool& critical) {
    critical = false;
    double f = 0;
    for (int it = 0; it < 12; ++it) {
        bool border = false;
        clampDomain(dom, u, v, border);
        f = fval(face, pln, u, v);
        if (std::abs(f) <= geomTol) return true;
        Vec3 Su, Sv;
        if (!face.derivUV(u, v, Su, Sv)) return false;
        const double fu = dot(pln.n, Su);
        const double fv = dot(pln.n, Sv);
        const double g2 = fu * fu + fv * fv;
        if (g2 < 1e-18) {
            critical = true;
            return std::abs(f) <= 10.0 * geomTol;
        }
        u -= f * fu / g2;
        v -= f * fv / g2;
    }
    return std::abs(f) <= 10.0 * geomTol;
}

bool uvTangent(const IFace& face, const Plane& pln, double u, double v, double& du, double& dv,
               Vec3& T3) {
    Vec3 Su, Sv;
    if (!face.derivUV(u, v, Su, Sv)) return false;
    const double fu = dot(pln.n, Su);
    const double fv = dot(pln.n, Sv);
    du = -fv;
    dv = fu;
    T3 = Su * du + Sv * dv;
    const double m = length(T3);
    if (m < 1e-18) return false;
    T3 = T3 / m;
    // Scale (du,dv) so |Su du + Sv dv| = 1.
    du /= m;
    dv /= m;
    return true;
}

double uvDist(const UVBox& d, double u0, double v0, double u1, double v1) {
    const double du = uvDelta(u0, u1, d.periodic_u, d.period_u);
    const double dv = uvDelta(v0, v1, d.periodic_v, d.period_v);
    return std::hypot(du, dv);
}

std::vector<UvPt> march(const IFace& face, const Plane& pln, const UVBox& dom, double u, double v,
                        int dir, const SliceOptions& opt, double diag, bool try_close) {
    std::vector<UvPt> pts;
    bool crit = false;
    if (!newton(face, pln, dom, u, v, opt.geom_tolerance, crit)) return pts;
    if (crit) {
        pts.push_back({u, v, face.evalUV(u, v)});
        return pts;
    }
    if (!inOrOn(face.classifyUV(u, v, opt.tolerance))) return pts;

    const double stitch = opt.tolerance;
    // Prefer larger 3D steps on large faces; FaceClassifier cost dominates tiny steps.
    const double stepMax = std::min(1.0, std::max(2e-2, 0.02 * diag));
    double step = stepMax;
    const double stepMin = std::max(opt.geom_tolerance * 10.0, 1e-3);
    const double u0 = u, v0 = v;
    const Vec3 p0 = face.evalUV(u, v);
    double arclen = 0;

    Vec3 Tprev{};
    bool haveT = false;

    for (int k = 0; k < 250; ++k) {
        const Vec3 p = face.evalUV(u, v);
        if (pts.empty() || dist(p, pts.back().p) > 0.25 * stepMin) {
            if (!pts.empty()) arclen += dist(p, pts.back().p);
            pts.push_back({u, v, p});
        }

        if (try_close && k > 12 && arclen > 4.0 * stitch) {
            const double gap = dist(p, p0);
            if (gap < 2.0 * step) step = std::max(stepMin, 0.35 * gap);
            if (gap <= stitch) {
                pts.push_back({u0, v0, p0});
                break;
            }
        }
        // Open (boundary) traces should not wander forever on huge UV domains.
        if (!try_close && arclen > 2.0 * std::max(diag, 1.0)) break;

        double du = 0, dv = 0;
        Vec3 T3;
        if (!uvTangent(face, pln, u, v, du, dv, T3)) break;
        du *= dir;
        dv *= dir;
        T3 = T3 * dir;
        if (haveT && dot(T3, Tprev) < 0) {
            du = -du;
            dv = -dv;
            T3 = T3 * -1.0;
        }

        bool advanced = false;
        for (int attempt = 0; attempt < 8; ++attempt) {
            double un = u + du * step;
            double vn = v + dv * step;
            bool border = false;
            clampDomain(dom, un, vn, border);
            bool crit2 = false;
            if (!newton(face, pln, dom, un, vn, opt.geom_tolerance, crit2)) {
                step *= 0.5;
                if (step < stepMin) break;
                continue;
            }
            if (crit2) {
                u = un;
                v = vn;
                advanced = true;
                break;
            }
            // Classify most steps; FClass2d is cached on OccFace.
            const bool checkTrim = (k & 1) == 0 || border || attempt > 0;
            if (checkTrim && !inOrOn(face.classifyUV(un, vn, opt.tolerance))) {
                double lo = 0, hi = 1;
                double ub = u, vb = v;
                for (int b = 0; b < 20; ++b) {
                    const double m = 0.5 * (lo + hi);
                    double ut = u + du * step * m;
                    double vt = v + dv * step * m;
                    clampDomain(dom, ut, vt, border);
                    newton(face, pln, dom, ut, vt, opt.geom_tolerance, crit2);
                    if (inOrOn(face.classifyUV(ut, vt, opt.tolerance))) {
                        lo = m;
                        ub = ut;
                        vb = vt;
                    } else {
                        hi = m;
                    }
                }
                u = ub;
                v = vb;
                pts.push_back({u, v, face.evalUV(u, v)});
                return pts;
            }

            Vec3 Tn;
            double du2, dv2;
            uvTangent(face, pln, un, vn, du2, dv2, Tn);
            const double cang = std::abs(dot(T3, Tn));
            if (cang < std::cos(25.0 * kPi / 180.0) && step > stepMin) {
                step *= 0.5;
                continue;
            }
            if (cang > std::cos(4.0 * kPi / 180.0)) step = std::min(stepMax, step * 1.15);
            u = un;
            v = vn;
            Tprev = Tn;
            haveT = true;
            advanced = true;
            break;
        }
        if (!advanced) break;
    }
    return pts;
}

bool nearTrace(const UVBox& d, const std::vector<UvPt>& tr, double u, double v, double tol) {
    for (const UvPt& q : tr) {
        if (uvDist(d, u, v, q.u, q.v) <= tol) return true;
    }
    return false;
}

RawSegment toSegment(const FaceRecord& iface, const std::vector<UvPt>& tr, const SliceOptions& opt) {
    RawSegment rs;
    rs.solid_id = iface.solid_id;
    rs.shell_id = iface.shell_id;
    rs.face_id = iface.face_id;
    if (tr.size() == 1) {
        rs.degenerate = true;
        rs.degen_event = "tangent_point";
        rs.degen_point = tr[0].p;
        return rs;
    }
    std::vector<Vec3> pts;
    pts.reserve(tr.size());
    for (const UvPt& q : tr) {
        if (pts.empty() || dist(q.p, pts.back()) > 0.1 * opt.geom_tolerance) pts.push_back(q.p);
    }
    if (pts.size() < 2) {
        rs.degenerate = true;
        rs.degen_event = "tangent_point";
        rs.degen_point = tr.front().p;
        return rs;
    }
    const bool closed = dist(pts.front(), pts.back()) <= opt.tolerance && pts.size() >= 4;
    rs.closed_loop = closed;

    rs.geom = fitCubicBSpline(pts, closed, opt.tolerance);
    if (rs.geom.ctrl_pts.size() < 4) {
        rs.geom.type = SegmentType::Line;
        rs.geom.start = pts.front();
        rs.geom.end = pts.back();
        rs.closed_loop = false;
        return rs;
    }
    // March points already lie on the surface; skip expensive invertUV true-error sampling.
    if (rs.geom.fit_error > opt.tolerance) {
        rs.geom = cubicFromPolyline(pts, closed);
    }
    return rs;
}

}  // namespace

std::vector<RawSegment> intersectNurbsFaceWithPlane(const FaceRecord& iface, const Plane& pln,
                                                    const SliceFrame& /*frame*/,
                                                    const SliceOptions& opt,
                                                    const std::vector<Vec3>& boundary3d) {
    std::vector<RawSegment> segs;
    if (!iface.face) return segs;
    const IFace& face = *iface.face;
    const UVBox dom = face.uvDomain();
    const double du = std::max(1e-16, dom.umax - dom.umin);
    const double dv = std::max(1e-16, dom.vmax - dom.vmin);
    const BBox bb = face.bbox();
    const double diag = bb.valid ? std::hypot(bb.xmax - bb.xmin, std::hypot(bb.ymax - bb.ymin,
                                                                            bb.zmax - bb.zmin))
                                 : 1.0;

    std::vector<Seed> seeds;
    int n_seeds_bound = 0;
    int n_seeds_grid = 0;

    auto addSeed = [&](double u, double v, bool boundary) {
        bool crit = false;
        if (!newton(face, pln, dom, u, v, opt.geom_tolerance, crit)) return;
        if (!inOrOn(face.classifyUV(u, v, opt.tolerance))) return;
        for (const Seed& s : seeds) {
            if (uvDist(dom, u, v, s.u, s.v) < 1e-4 * (du + dv)) return;
        }
        seeds.push_back({u, v, boundary, false});
        if (boundary) ++n_seeds_bound;
        else ++n_seeds_grid;
    };

    for (const Vec3& p : boundary3d) {
        double u = 0, v = 0;
        if (face.invertUV(p, u, v, opt.tolerance)) addSeed(u, v, true);
    }

    const int N = 16;
    // Coarse probe first: AABB false-positives used to pay for a full dense evalUV grid
    // even when f never changes sign.
    const int Nc = 8;
    std::vector<double> Fc((Nc + 1) * (Nc + 1));
    auto atc = [&](int i, int j) -> double& { return Fc[i * (Nc + 1) + j]; };
    bool coarseHit = false;
    int nearPlane = 0;
    const int nCoarse = (Nc + 1) * (Nc + 1);
    const double nearTol = std::max(10.0 * opt.tolerance, 1e-3);
    for (int i = 0; i <= Nc; ++i) {
        for (int j = 0; j <= Nc; ++j) {
            const double u = dom.umin + du * (static_cast<double>(i) / Nc);
            const double v = dom.vmin + dv * (static_cast<double>(j) / Nc);
            atc(i, j) = fval(face, pln, u, v);
            if (std::abs(atc(i, j)) <= nearTol) ++nearPlane;
            if (std::abs(atc(i, j)) <= opt.tolerance) coarseHit = true;
        }
    }
    // Nearly coplanar NURBS×plane: prefer boundary hits only.
    const bool nearlyCoplanar = nearPlane * 5 >= nCoarse * 3;  // >=60% samples near plane
    if (!nearlyCoplanar && !coarseHit) {
        for (int i = 0; i <= Nc && !coarseHit; ++i) {
            for (int j = 0; j <= Nc && !coarseHit; ++j) {
                if (i < Nc && atc(i, j) * atc(i + 1, j) < 0 &&
                    std::max(std::abs(atc(i, j)), std::abs(atc(i + 1, j))) > opt.tolerance)
                    coarseHit = true;
                if (j < Nc && atc(i, j) * atc(i, j + 1) < 0 &&
                    std::max(std::abs(atc(i, j)), std::abs(atc(i, j + 1))) > opt.tolerance)
                    coarseHit = true;
            }
        }
    } else if (nearlyCoplanar) {
        coarseHit = false;
    }

    // If trim-boundary already provides plane hits, skip the dense UV seed grid.
    if (coarseHit && n_seeds_bound >= 2) coarseHit = false;

    if (coarseHit) {
        std::vector<double> F((N + 1) * (N + 1));
        auto at = [&](int i, int j) -> double& { return F[i * (N + 1) + j]; };
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                const double u = dom.umin + du * (static_cast<double>(i) / N);
                const double v = dom.vmin + dv * (static_cast<double>(j) / N);
                at(i, j) = fval(face, pln, u, v);
            }
        }
        auto consider = [&](int i0, int j0, int i1, int j1) {
            if (n_seeds_grid >= 8) return;
            const double f0 = at(i0, j0);
            const double f1 = at(i1, j1);
            if (f0 * f1 > 0) return;
            if (std::max(std::abs(f0), std::abs(f1)) <= opt.tolerance) return;
            double a = 0, b = 1;
            for (int k = 0; k < 18; ++k) {
                const double m = 0.5 * (a + b);
                const double u = dom.umin + du * ((i0 + m * (i1 - i0)) / N);
                const double v = dom.vmin + dv * ((j0 + m * (j1 - j0)) / N);
                const double fm = fval(face, pln, u, v);
                if (f0 * fm <= 0) b = m;
                else a = m;
            }
            const double m = 0.5 * (a + b);
            const double u = dom.umin + du * ((i0 + m * (i1 - i0)) / N);
            const double v = dom.vmin + dv * ((j0 + m * (j1 - j0)) / N);
            addSeed(u, v, false);
        };
        for (int i = 0; i <= N; ++i) {
            for (int j = 0; j <= N; ++j) {
                if (i < N) consider(i, j, i + 1, j);
                if (j < N) consider(i, j, i, j + 1);
            }
        }
    }

    const double cover = 1.25 * std::hypot(du / N, dv / N);
    std::vector<std::vector<UvPt>> traces;

    auto consume = [&](const std::vector<UvPt>& tr) {
        if (tr.empty()) return;
        traces.push_back(tr);
        for (Seed& s : seeds) {
            if (!s.used && nearTrace(dom, tr, s.u, s.v, cover)) s.used = true;
        }
    };

    for (Seed& s : seeds) {
        if (s.used) continue;
        if (s.from_boundary) {
            auto a = march(face, pln, dom, s.u, s.v, +1, opt, diag, false);
            auto b = march(face, pln, dom, s.u, s.v, -1, opt, diag, false);
            std::vector<UvPt> tr;
            for (int i = static_cast<int>(a.size()) - 1; i >= 0; --i) tr.push_back(a[i]);
            for (size_t i = 1; i < b.size(); ++i) tr.push_back(b[i]);
            if (tr.size() >= 2) consume(tr);
            else if (tr.size() == 1) consume(tr);
            s.used = true;
        }
    }

    for (Seed& s : seeds) {
        if (s.used) continue;
        auto tr = march(face, pln, dom, s.u, s.v, +1, opt, diag, true);
        if (tr.size() < 2) {
            auto tr2 = march(face, pln, dom, s.u, s.v, -1, opt, diag, true);
            if (tr2.size() > tr.size()) tr = std::move(tr2);
        }
        consume(tr);
        s.used = true;
    }

    for (const auto& tr : traces) {
        segs.push_back(toSegment(iface, tr, opt));
    }
    return segs;
}

}  // namespace brepslicer
