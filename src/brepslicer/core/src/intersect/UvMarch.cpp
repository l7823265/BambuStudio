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
    // Cap 3D step so chord sagitta h^2/(8r) stays under ~1e-4 for r on the order of 1 mm.
    const double stepMax = std::min(0.012, std::max(5e-4, 0.001 * diag));
    double step = stepMax;
    const double stepMin = std::max(opt.geom_tolerance * 10.0, 1e-5);
    const double u0 = u, v0 = v;
    const Vec3 p0 = face.evalUV(u, v);
    double arclen = 0;

    Vec3 Tprev{};
    bool haveT = false;

    for (int k = 0; k < 8000; ++k) {
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
            if (!inOrOn(face.classifyUV(un, vn, opt.tolerance))) {
                // Backtrack to the trim / domain boundary.
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

RawSegment toSegment(const FaceRecord& iface, const IFace& face, const Plane& pln, const UVBox& dom,
                     const std::vector<UvPt>& tr, const SliceOptions& opt) {
    RawSegment rs;
    rs.solid_id = iface.solid_id;
    rs.shell_id = iface.shell_id;
    rs.face_id = iface.face_id;
    rs.chain_idx = 0;
    rs.geom.face_id = iface.face_id;
    rs.geom.chain_idx = 0;
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

    auto trueErr = [&](const Segment& geom) {
        std::vector<Vec3> spl;
        const int ns = std::max(96, static_cast<int>(geom.ctrl_pts.size()));
        bsplineSample(geom, ns, spl);
        double e = 0;
        int hits = 0;
        for (const Vec3& p : spl) {
            double u = 0, v = 0;
            if (!face.invertUV(p, u, v, opt.tolerance)) continue;
            bool crit = false;
            newton(face, pln, dom, u, v, opt.geom_tolerance, crit);
            e = std::max(e, dist(p, face.evalUV(u, v)));
            ++hits;
        }
        return (hits > 0) ? e : geom.fit_error;
    };

    rs.geom = fitCubicBSpline(pts, closed, opt.tolerance);
    if (rs.geom.ctrl_pts.size() < 4) {
        rs.geom.type = SegmentType::Line;
        rs.geom.start = pts.front();
        rs.geom.end = pts.back();
        rs.closed_loop = false;
        return rs;
    }
    double te = trueErr(rs.geom);
    if (te > opt.tolerance) {
        // 3D cubics through on-surface points leave the surface; stay on the march chords.
        rs.geom = cubicFromPolyline(pts, closed);
        te = trueErr(rs.geom);
    }
    rs.geom.fit_error = std::max(rs.geom.fit_error, te);
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
    double dmin = 0, dmax = 0;
    if (bb.valid) aabbRangeAlong(bb, pln.n, dmin, dmax);
    const double diag = bb.valid ? std::hypot(bb.xmax - bb.xmin, std::hypot(bb.ymax - bb.ymin,
                                                                            bb.zmax - bb.zmin))
                                 : 1.0;

    std::vector<Seed> seeds;

    auto addSeed = [&](double u, double v, bool boundary) {
        bool crit = false;
        if (!newton(face, pln, dom, u, v, opt.geom_tolerance, crit)) return;
        if (!inOrOn(face.classifyUV(u, v, opt.tolerance))) return;
        const Vec3 p = face.evalUV(u, v);
        if (!inOrOn(face.classify(p, opt.tolerance))) return;
        for (const Seed& s : seeds) {
            if (uvDist(dom, u, v, s.u, s.v) < 1e-4 * (du + dv)) return;
        }
        seeds.push_back({u, v, boundary, false});
    };

    for (const Vec3& p : boundary3d) {
        double u = 0, v = 0;
        if (face.invertUV(p, u, v, opt.tolerance)) addSeed(u, v, true);
    }

    // OCC face∩plane samples (same idea as BRepAlgoAPI_Section on the trimmed face).
    // Critical for narrow trim ribbons that a coarse UV grid misses.
    {
        const double spacing = std::max(0.25, 0.02 * diag);
        const std::vector<Vec3> occ_hits = face.samplePlaneSection(pln, spacing);
        for (const Vec3& p : occ_hits) {
            double u = 0, v = 0;
            if (!face.invertUV(p, u, v, opt.tolerance)) continue;
            addSeed(u, v, false);
        }
    }

    // Adaptive UV grid: denser when the face is a thin slab / few edge hits (narrow trim).
    int N = 32;
    if (boundary3d.size() <= 4) N = 48;
    if (bb.valid) {
        const double thin = std::min({bb.xmax - bb.xmin, bb.ymax - bb.ymin, bb.zmax - bb.zmin});
        const double thick = std::max({bb.xmax - bb.xmin, bb.ymax - bb.ymin, bb.zmax - bb.zmin});
        if (thin > 1e-6 && thin < 0.25 * thick) N = std::max(N, 64);
    }
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
        const double f0 = at(i0, j0);
        const double f1 = at(i1, j1);
        if (f0 * f1 > 0 && std::abs(f0) > opt.tolerance && std::abs(f1) > opt.tolerance) return;
        if (f0 * f1 > 0) return;
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
        const PointClass c0 = face.classifyUV(dom.umin + du * i0 / N, dom.vmin + dv * j0 / N,
                                              opt.tolerance);
        const PointClass c1 = face.classifyUV(dom.umin + du * i1 / N, dom.vmin + dv * j1 / N,
                                              opt.tolerance);
        if (!inOrOn(c0) && !inOrOn(c1)) return;
        const double um = dom.umin + du * ((i0 + m * (i1 - i0)) / N);
        const double vm = dom.vmin + dv * ((j0 + m * (j1 - j0)) / N);
        if (!inOrOn(face.classifyUV(um, vm, opt.tolerance))) return;
        if (!inOrOn(face.classify(face.evalUV(um, vm), opt.tolerance))) return;
        addSeed(u, v, false);
    };
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            if (i < N) consider(i, j, i + 1, j);
            if (j < N) consider(i, j, i, j + 1);
        }
    }

    // Extra zero-cross probes on isos through inverted boundary UVs (narrow ribbon).
    for (const Vec3& p : boundary3d) {
        double ub = 0, vb = 0;
        if (!face.invertUV(p, ub, vb, opt.tolerance)) continue;
        if (dom.periodic_u || dom.periodic_v) wrapUV(dom, ub, vb);
        constexpr int kProbe = 48;
        for (int axis = 0; axis < 2; ++axis) {
            for (int i = 0; i < kProbe; ++i) {
                const double t0 = static_cast<double>(i) / kProbe;
                const double t1 = static_cast<double>(i + 1) / kProbe;
                double u0, v0, u1, v1;
                if (axis == 0) {
                    u0 = dom.umin + du * t0;
                    u1 = dom.umin + du * t1;
                    v0 = v1 = vb;
                } else {
                    v0 = dom.vmin + dv * t0;
                    v1 = dom.vmin + dv * t1;
                    u0 = u1 = ub;
                }
                const double f0 = fval(face, pln, u0, v0);
                const double f1 = fval(face, pln, u1, v1);
                if (f0 * f1 > 0.0 && std::abs(f0) > opt.geom_tolerance &&
                    std::abs(f1) > opt.geom_tolerance)
                    continue;
                double lo = 0, hi = 1, flo = f0;
                for (int it = 0; it < 20; ++it) {
                    const double m = 0.5 * (lo + hi);
                    const double um = u0 + m * (u1 - u0);
                    const double vm = v0 + m * (v1 - v0);
                    const double fm = fval(face, pln, um, vm);
                    if (flo * fm <= 0.0) hi = m;
                    else {
                        lo = m;
                        flo = fm;
                    }
                }
                const double m = 0.5 * (lo + hi);
                addSeed(u0 + m * (u1 - u0), v0 + m * (v1 - v0), false);
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
        RawSegment rs = toSegment(iface, face, pln, dom, tr, opt);
        segs.push_back(std::move(rs));
    }
    (void)dmin;
    (void)dmax;
    return segs;
}

}  // namespace brepslicer
