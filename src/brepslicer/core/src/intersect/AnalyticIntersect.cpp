#include <intersect/AnalyticIntersect.h>
#include <intersect/BSplineFit.h>
#include <intersect/QuadricPlane.h>
#include <intersect/UvMarch.h>
#include <intersect/UvMatch.h>
#include <geom/GeomUtil.h>
#include <topo/SolidAdjacency.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <vector>

namespace brepslicer {
namespace {

void reverseSegment(Segment& s) {
    std::swap(s.start, s.end);
    if (s.type == SegmentType::Arc || s.type == SegmentType::Ellipse) {
        s.start_angle = s.start_angle + s.sweep;
        s.sweep = -s.sweep;
    }
    if (s.type == SegmentType::BSpline) reverseBSpline(s);
    else std::reverse(s.ctrl_pts.begin(), s.ctrl_pts.end());
}

Segment makeLineSeg(const Vec3& a, const Vec3& b) {
    Segment s;
    s.type = SegmentType::Line;
    s.start = a;
    s.end = b;
    return s;
}

Segment makeArc(const Circle3& c, double t0, double t1, const SliceFrame& frame) {
    const Vec3 p0 = c.at(t0);
    const Vec3 p1 = c.at(t1);

    // Honor the parametric span from the caller (trimCurveToFace may select the
    // major arc). Prefer the mid-point that matches 0.5*(t0+t1); do NOT fold
    // |dt|>π down to the short chord (that undid ROBOT z≈70 cone wrap trims).
    double dt = t1 - t0;
    while (dt > kTwoPi) dt -= kTwoPi;
    while (dt < -kTwoPi) dt += kTwoPi;
    if (std::abs(dt) < 1e-15) {
        // Full turn request (t1 = t0 + 2π): keep closed circle.
        if (std::abs((t1 - t0) - kTwoPi) <= 1e-9 || std::abs((t1 - t0) + kTwoPi) <= 1e-9)
            dt = (t1 >= t0) ? kTwoPi : -kTwoPi;
    }
    const double alt = dt > 0.0 ? dt - kTwoPi : dt + kTwoPi;
    if (std::abs(alt) > 1e-15 && std::abs(alt) < kTwoPi - 1e-15) {
        const Vec3 refMid = c.at(0.5 * (t0 + t1));
        auto midErr = [&](double d) { return length(c.at(t0 + 0.5 * d) - refMid); };
        if (midErr(alt) + 1e-10 < midErr(dt)) dt = alt;
    }

    const Vec3 mid = c.at(t0 + 0.5 * dt);
    Segment s;
    s.type = SegmentType::Arc;
    s.start = p0;
    s.end = p1;
    s.center = c.center;
    s.normal = frame.n;
    s.radius = c.radius;
    s.start_angle = frame.angle(c.center, p0);
    const double mag = std::abs(dt);
    const double orient = dot(cross(p0 - c.center, mid - c.center), frame.n);
    s.sweep = (orient >= 0.0 ? mag : -mag);
    return s;
}

Segment makeEllipseSeg(const Ellipse3& e, double t0, double t1, const SliceFrame& frame) {
    const Vec3 p0 = e.at(t0);
    const Vec3 p1 = e.at(t1);

    // Same span logic as makeArc: honor trimCurveToFace's parametric interval
    // (major wrap or minor chord). Do not fold |dt|>π to the short lobe.
    double dt = t1 - t0;
    while (dt > kTwoPi) dt -= kTwoPi;
    while (dt < -kTwoPi) dt += kTwoPi;
    if (std::abs(dt) < 1e-15) {
        if (std::abs((t1 - t0) - kTwoPi) <= 1e-9 || std::abs((t1 - t0) + kTwoPi) <= 1e-9)
            dt = (t1 >= t0) ? kTwoPi : -kTwoPi;
    }
    const double alt = dt > 0.0 ? dt - kTwoPi : dt + kTwoPi;
    if (std::abs(alt) > 1e-15 && std::abs(alt) < kTwoPi - 1e-15) {
        const Vec3 refMid = e.at(0.5 * (t0 + t1));
        auto midErr = [&](double d) { return length(e.at(t0 + 0.5 * d) - refMid); };
        if (midErr(alt) + 1e-10 < midErr(dt)) dt = alt;
    }

    Segment s;
    s.type = SegmentType::Ellipse;
    s.start = p0;
    s.end = p1;
    s.center = e.center;
    s.normal = frame.n;
    s.major_axis = e.major;
    s.radius = e.a;
    s.radius_b = e.b;
    s.start_angle = t0;
    s.sweep = dt;
    return s;
}

bool pointOnPlane(const Vec3& p, const Plane& pln, double tol) {
    return std::abs(signedPlaneDist(pln, p)) <= tol;
}

// UV face classifier with period shifts — OCC often returns Out at the seam even when
// the 3D point lies on the trimmed cylinder/torus (ROBOT z≈373 face 889).
bool uvInOrOnPeriodic(const IFace& face, double u, double v, double tol) {
    const UVBox dom = face.uvDomain();
    wrapUV(dom, u, v);
    auto ok = [&](double uu, double vv) {
        const PointClass c = face.classifyUV(uu, vv, tol);
        return c == PointClass::In || c == PointClass::On;
    };
    if (ok(u, v)) return true;
    if (dom.periodic_u && dom.period_u > 0) {
        if (ok(u + dom.period_u, v) || ok(u - dom.period_u, v)) return true;
        // Face UV box may be a proper subset of one period (trimmed cylinder).
        double uw = u;
        if (uw < dom.umin) {
            while (uw < dom.umin) uw += dom.period_u;
        } else if (uw > dom.umax) {
            while (uw > dom.umax) uw -= dom.period_u;
        }
        if (std::abs(uw - u) > 1e-14 && ok(uw, v)) return true;
    }
    if (dom.periodic_v && dom.period_v > 0) {
        if (ok(u, v + dom.period_v) || ok(u, v - dom.period_v)) return true;
        double vw = v;
        if (vw < dom.vmin) {
            while (vw < dom.vmin) vw += dom.period_v;
        } else if (vw > dom.vmax) {
            while (vw > dom.vmax) vw -= dom.period_v;
        }
        if (std::abs(vw - v) > 1e-14 && ok(u, vw)) return true;
    }
    return false;
}

bool classifyInOrOnFace(const IFace& face, const Vec3& p, double tol) {
    // Keep In|On for trim mid-tests: OCC often marks near-boundary samples as On.
    // Shared-edge single-owner is enforced at assemble (coplanar / lower face_id).
    const PointClass st = face.classify(p, tol);
    if (st != PointClass::In && st != PointClass::On) return false;
    // 3D classifier can keep points on the infinite cylinder/cone outside the wire.
    // Require the UV projection to lie in/on the trimmed face (ROBOT_4 z≈202 orphans).
    double u = 0, v = 0;
    if (!face.invertUV(p, u, v, tol)) return false;
    return uvInOrOnPeriodic(face, u, v, tol);
}

void collectEdgePlanePoints(const IEdge& edge, const Plane& pln, double angTol, double tol,
                            std::vector<Vec3>& pts) {
    const CurveData c = edge.curve();
    if (c.degenerated) return;
    const Vec3 pF = edge.eval(c.first);
    const Vec3 pL = edge.eval(c.last);
    const Vec3 pM = edge.eval(0.5 * (c.first + c.last));
    const bool inPlane = pointOnPlane(pF, pln, tol) && pointOnPlane(pL, pln, tol) &&
                         pointOnPlane(pM, pln, tol);
    if (inPlane) {
        pts.push_back(pF);
        pts.push_back(pL);
        return;
    }
    intersectCurvePlane(c, pln, angTol, tol, pts, [&](double t) { return edge.eval(t); });
}

std::vector<Vec3> faceBoundaryHits(const IFace& face, const Plane& pln, double angTol, double tol) {
    std::vector<Vec3> pts;
    for (const auto& e : face.edges()) {
        if (e) collectEdgePlanePoints(*e, pln, angTol, tol, pts);
    }
    // Vertex-shared edges often report the same plane hit twice. Only merge true
    // duplicates — a hard 0.5 mm floor collapsed distinct tips on narrow faces
    // (WithSphere z≈69.37 F44: two hits 0.139 mm apart → one constraint, open gap to F51).
    const double merge = std::max(tol * 20.0, 1e-3);
    std::vector<Vec3> uniq;
    for (const Vec3& p : pts) {
        bool dup = false;
        for (const Vec3& q : uniq) {
            if (dist(p, q) <= merge) {
                dup = true;
                break;
            }
        }
        if (!dup) uniq.push_back(p);
    }
    return uniq;
}

// Debug: log when edge∩plane merge drops hits (distinct tips collapsed).
void auditFaceBoundaryHits(int face_id, const IFace& face, const Plane& pln, double angTol,
                           double tol, std::vector<std::string>* audit) {
    if (!audit) return;
    std::vector<Vec3> raw;
    for (const auto& e : face.edges()) {
        if (e) collectEdgePlanePoints(*e, pln, angTol, tol, raw);
    }
    const double merge = std::max(tol * 20.0, 1e-3);
    std::vector<Vec3> uniq;
    for (const Vec3& p : raw) {
        bool dup = false;
        for (const Vec3& q : uniq) {
            if (dist(p, q) <= merge) {
                dup = true;
                break;
            }
        }
        if (!dup) uniq.push_back(p);
    }
    if (raw.size() <= uniq.size()) return;
    std::ostringstream os;
    os << std::setprecision(9) << "boundary_merge face=" << face_id << " raw=" << raw.size()
       << " uniq=" << uniq.size() << " merge_tol=" << merge;
    for (size_t i = 0; i < raw.size(); ++i) {
        for (size_t j = i + 1; j < raw.size(); ++j) {
            const double d = dist(raw[i], raw[j]);
            if (d <= 1.0) os << " d" << i << j << "=" << d;
        }
    }
    audit->push_back(os.str());
}

Vec3 planeXYPoint(const SliceFrame& frame, const Plane& pln, const Vec2& q) {
    return frame.x * q.x + frame.y * q.y + pln.n * pln.d;
}

// Analytic circle–circle intersection on the slice plane (constraint points for free-face trim).
std::vector<Vec3> circleCircleIntersect(const Circle3& a, const Circle3& b, const Plane& pln,
                                        const SliceFrame& frame, double tol) {
    std::vector<Vec3> out;
    if (!almostParallel(a.n, b.n, std::max(tol, 1e-8))) return out;
    if (dist(a.center, b.center) <= tol && std::abs(a.radius - b.radius) <= tol) return out;

    const Vec2 ca = frame.toXY(a.center);
    const Vec2 cb = frame.toXY(b.center);
    const double r0 = a.radius;
    const double r1 = b.radius;
    const double dx = cb.x - ca.x;
    const double dy = cb.y - ca.y;
    const double d = std::hypot(dx, dy);
    if (d <= tol) return out;

    const double reach = r0 + r1;
    const double reach_min = std::abs(r0 - r1);
    if (d > reach + tol || d < reach_min - tol) return out;

    const double a_len = (r0 * r0 - r1 * r1 + d * d) / (2.0 * d);
    const double h2 = r0 * r0 - a_len * a_len;
    if (h2 < -tol * (1.0 + d)) return out;

    const double inv_d = 1.0 / d;
    const double xm = ca.x + a_len * dx * inv_d;
    const double ym = ca.y + a_len * dy * inv_d;
    if (h2 <= tol * tol) {
        out.push_back(planeXYPoint(frame, pln, {xm, ym}));
        return out;
    }
    const double h = std::sqrt(h2);
    const double rx = -dy * h * inv_d;
    const double ry = dx * h * inv_d;
    out.push_back(planeXYPoint(frame, pln, {xm + rx, ym + ry}));
    out.push_back(planeXYPoint(frame, pln, {xm - rx, ym - ry}));
    return out;
}

struct AnalyticCurve {
    enum Kind { Line, Circle, Ellipse } kind = Line;
    Line3 line;
    Circle3 circle;
    Ellipse3 ellipse;
    bool periodic = false;
    double period = kTwoPi;
};

bool snapToCurve(const AnalyticCurve& c, const Vec3& p, double tol, double& t) {
    if (c.kind == AnalyticCurve::Line) {
        t = dot(p - c.line.p, c.line.d);
        return length(c.line.at(t) - p) <= tol;
    }
    if (c.kind == AnalyticCurve::Circle) {
        t = c.circle.param(p);
        return length(c.circle.at(t) - p) <= tol;
    }
    t = c.ellipse.param(p);
    return length(c.ellipse.at(t) - p) <= tol;
}

// Project onto the analytic curve; returns distance to the closest curve point.
double projectToCurve(const AnalyticCurve& c, const Vec3& p, double& t) {
    if (c.kind == AnalyticCurve::Line) {
        t = dot(p - c.line.p, c.line.d);
        return length(c.line.at(t) - p);
    }
    if (c.kind == AnalyticCurve::Circle) {
        t = c.circle.param(p);
        return length(c.circle.at(t) - p);
    }
    t = c.ellipse.param(p);
    return length(c.ellipse.at(t) - p);
}

Vec3 evalCurve(const AnalyticCurve& c, double t) {
    if (c.kind == AnalyticCurve::Line) return c.line.at(t);
    if (c.kind == AnalyticCurve::Circle) return c.circle.at(t);
    return c.ellipse.at(t);
}

// Fitting BREP edges between a free/NURBS face and an analytic wall drift off the
// true plane∩analytic tip (ROBOT_4 z≈278: F257 edge∩plane vs cyl generator ~0.27 mm).
// Replace edge∩plane constraints with the min-dist projection onto each shared-edge
// neighbor's analytic∩plane curve (Approach B; coincides with reintersect for lines).
void snapBoundaryHitsToAnalyticNeighbors(std::vector<Vec3>& hits, const FaceRecord& iface,
                                         const Plane& pln, const SliceOptions& opt) {
    if (hits.empty() || !opt.plane_faces || !iface.face) return;
    const double fit_lateral = std::max(0.05, 500.0 * opt.tolerance);
    const SolidAdjacency* adj =
        (opt.solid_adjacency && !opt.solid_adjacency->empty()) ? opt.solid_adjacency : nullptr;

    struct NbCurve {
        AnalyticCurve c;
        int face_id = -1;
    };
    std::vector<NbCurve> curves;
    curves.reserve(8);
    for (const FaceRecord* nb : *opt.plane_faces) {
        if (!nb || !nb->face || nb->solid_id != iface.solid_id) continue;
        if (nb->face_id == iface.face_id) continue;
        if (adj && !adj->facesShareEdge(iface.face_id, nb->face_id)) continue;
        const SurfaceData surf = nb->face->surface();
        if (surf.kind == SurfaceKind::Other) continue;
        const std::vector<QuadricCurve> qcurves = intersectAnalyticSurfacePlane(
            surf, pln, opt.angular_tolerance, opt.geom_tolerance);
        for (const QuadricCurve& q : qcurves) {
            if (q.degenerate) continue;
            NbCurve nc;
            nc.face_id = nb->face_id;
            if (q.kind == QuadricKind::Line) {
                nc.c.kind = AnalyticCurve::Line;
                nc.c.line = q.line;
            } else if (q.kind == QuadricKind::Circle) {
                nc.c.kind = AnalyticCurve::Circle;
                nc.c.circle = q.circle;
                nc.c.periodic = true;
                nc.c.period = kTwoPi;
            } else if (q.kind == QuadricKind::Ellipse) {
                nc.c.kind = AnalyticCurve::Ellipse;
                nc.c.ellipse = q.ellipse;
                nc.c.periodic = true;
                nc.c.period = kTwoPi;
            } else {
                continue;
            }
            curves.push_back(nc);
        }
    }
    if (curves.empty()) return;

    for (Vec3& p : hits) {
        double best_d = fit_lateral;
        Vec3 best_p = p;
        bool found = false;
        for (const NbCurve& nc : curves) {
            double t = 0;
            const double d = projectToCurve(nc.c, p, t);
            if (d > best_d) continue;
            best_d = d;
            best_p = evalCurve(nc.c, t);
            found = true;
        }
        if (found && dist(p, best_p) > opt.geom_tolerance) {
            if (opt.constraint_audit) {
                std::ostringstream os;
                os << std::setprecision(9) << "constraint_snap_analytic face=" << iface.face_id
                   << " from=(" << p.x << "," << p.y << "," << p.z << ") to=(" << best_p.x << ","
                   << best_p.y << "," << best_p.z << ") derr=" << best_d;
                opt.constraint_audit->push_back(os.str());
            }
            p = best_p;
        }
    }
}

double arcLength(const AnalyticCurve& c, double t0, double t1) {
    const double dt = std::abs(t1 - t0);
    if (c.kind == AnalyticCurve::Line) return dt;
    if (c.kind == AnalyticCurve::Circle) return c.circle.radius * dt;
    return 0.5 * (c.ellipse.a + c.ellipse.b) * dt;
}

Segment emitInterval(const AnalyticCurve& c, double t0, double t1, const SliceFrame& frame) {
    if (c.kind == AnalyticCurve::Line) return makeLineSeg(evalCurve(c, t0), evalCurve(c, t1));
    if (c.kind == AnalyticCurve::Circle) return makeArc(c.circle, t0, t1, frame);
    return makeEllipseSeg(c.ellipse, t0, t1, frame);
}

// Project face AABB onto an infinite line → parameter range that can hit the face.
bool lineParamRangeOnBBox(const Line3& line, const BBox& bb, double pad, double& t0, double& t1) {
    if (!bb.valid) return false;
    const Vec3 corners[8] = {
        {bb.xmin, bb.ymin, bb.zmin}, {bb.xmax, bb.ymin, bb.zmin}, {bb.xmin, bb.ymax, bb.zmin},
        {bb.xmax, bb.ymax, bb.zmin}, {bb.xmin, bb.ymin, bb.zmax}, {bb.xmax, bb.ymin, bb.zmax},
        {bb.xmin, bb.ymax, bb.zmax}, {bb.xmax, bb.ymax, bb.zmax},
    };
    t0 = 1e300;
    t1 = -1e300;
    for (const Vec3& c : corners) {
        const double t = dot(c - line.p, line.d);
        t0 = std::min(t0, t);
        t1 = std::max(t1, t);
    }
    t0 -= pad;
    t1 += pad;
    return t1 > t0;
}

// When edge∩plane hits are missing/sparse (common off the mid-plane for // cylinders),
// recover In/On intervals by classifying samples along the line over the face bbox.
void classifyLineIntervals(const AnalyticCurve& curve, const IFace& face, double tol,
                           const std::vector<double>& splits, std::vector<std::pair<double, double>>& out) {
    out.clear();
    if (curve.kind != AnalyticCurve::Line) return;
    const BBox bb = face.bbox();
    double tLo = 0, tHi = 0;
    if (!lineParamRangeOnBBox(curve.line, bb, std::max(1.0, 50.0 * tol), tLo, tHi)) return;
    // When edge hits exist, do not scan outside their span — AABB pads past trim on
    // cylinders and produces exterior generator orphans.
    if (splits.size() >= 2) {
        double s0 = splits.front(), s1 = splits.front();
        for (double s : splits) {
            s0 = std::min(s0, s);
            s1 = std::max(s1, s);
        }
        tLo = std::max(tLo, s0);
        tHi = std::min(tHi, s1);
        if (tHi - tLo <= tol) return;
    }

    std::vector<double> knots = splits;
    knots.push_back(tLo);
    knots.push_back(tHi);
    std::sort(knots.begin(), knots.end());
    knots.erase(std::unique(knots.begin(), knots.end(),
                            [](double a, double b) { return std::abs(a - b) < 1e-12; }),
                knots.end());

    constexpr int kSamples = 48;
    std::vector<double> ts;
    ts.reserve(static_cast<size_t>(kSamples) + knots.size() * 2);
    for (int i = 0; i <= kSamples; ++i) {
        ts.push_back(tLo + (tHi - tLo) * (static_cast<double>(i) / kSamples));
    }
    for (double t : knots) {
        if (t >= tLo - 1e-14 && t <= tHi + 1e-14) ts.push_back(t);
    }
    // Midpoints between consecutive wire hits: hole Out bands sit between On-classified
    // boundary knots; uniform samples alone can skip a narrow notch (ROBOT z≈506 F209).
    for (size_t i = 0; i + 1 < knots.size(); ++i) {
        if (knots[i + 1] - knots[i] > tol) ts.push_back(0.5 * (knots[i] + knots[i + 1]));
    }
    std::sort(ts.begin(), ts.end());
    ts.erase(std::unique(ts.begin(), ts.end(),
                         [](double a, double b) { return std::abs(a - b) < 1e-14; }),
             ts.end());
    if (ts.size() < 2) return;

    auto inside = [&](double t) {
        return classifyInOrOnFace(face, evalCurve(curve, t), tol);
    };

    bool run = false;
    double run0 = 0;
    for (size_t i = 0; i < ts.size(); ++i) {
        const bool in = inside(ts[i]);
        if (in && !run) {
            run = true;
            run0 = ts[i];
        } else if (!in && run) {
            // End just before this Out sample.
            const double run1 = (i > 0) ? ts[i - 1] : ts[i];
            if (run1 - run0 > tol) out.push_back({run0, run1});
            run = false;
        }
    }
    if (run && ts.back() - run0 > tol) out.push_back({run0, ts.back()});
}

void appendUniquePoints(std::vector<Vec3>& dst, const std::vector<Vec3>& src, double tol) {
    for (const Vec3& p : src) {
        bool dup = false;
        for (const Vec3& q : dst) {
            if (dist(p, q) <= tol) {
                dup = true;
                break;
            }
        }
        if (!dup) dst.push_back(p);
    }
}

Circle3 segmentCircle(const Segment& s, const SliceFrame& frame) {
    Circle3 c;
    c.center = s.center;
    c.n = s.normal;
    c.xdir = frame.x;
    c.radius = s.radius;
    return c;
}

bool angleInSweepInterior(double a, double a0, double sweep) {
    constexpr double eps = 1e-9;
    double rel = a - a0;
    while (rel <= -kPi) rel += kTwoPi;
    while (rel > kPi) rel -= kTwoPi;
    if (sweep >= 0.0) return rel > eps && rel < sweep - eps;
    return rel < -eps && rel > sweep + eps;
}

void snapArcEndpoint(Segment& s, const Vec3& pt, bool is_start, const SliceFrame& frame) {
    const Circle3 c = segmentCircle(s, frame);
    const double a_hit = frame.angle(s.center, pt);
    const double old_abs = std::abs(s.sweep);
    const bool was_major = old_abs > kPi + 1e-9 && old_abs < kTwoPi - 1e-9;
    if (is_start) {
        const double a_end = s.start_angle + s.sweep;
        s.start_angle = a_hit;
        s.sweep = a_end - a_hit;
    } else {
        s.sweep = a_hit - s.start_angle;
    }
    // Preserve major vs minor: default fold-to-short only when the arc was minor.
    while (s.sweep > kTwoPi) s.sweep -= kTwoPi;
    while (s.sweep < -kTwoPi) s.sweep += kTwoPi;
    if (!was_major) {
        if (s.sweep > kPi) s.sweep -= kTwoPi;
        if (s.sweep < -kPi) s.sweep += kTwoPi;
    } else if (std::abs(s.sweep) <= kPi) {
        // Snap flipped a major into a minor — restore the long way.
        s.sweep = s.sweep > 0.0 ? s.sweep - kTwoPi : s.sweep + kTwoPi;
    }
    s.start = c.at(s.start_angle);
    s.end = c.at(s.start_angle + s.sweep);
}

void snapArcEndpoints(Segment& s, const std::vector<Vec3>& constraints, const SliceFrame& frame,
                      double snap_tol) {
    for (const Vec3& p : constraints) {
        if (std::abs(dist(p, s.center) - s.radius) > snap_tol) continue;
        if (dist(p, s.start) <= snap_tol) snapArcEndpoint(s, p, true, frame);
        if (dist(p, s.end) <= snap_tol) snapArcEndpoint(s, p, false, frame);
    }
}

std::vector<RawSegment> splitArcAtConstraints(RawSegment rs, const std::vector<Vec3>& constraints,
                                              const SliceFrame& frame, double snap_tol) {
    if (rs.geom.type != SegmentType::Arc || constraints.empty()) {
        if (rs.geom.type == SegmentType::Arc && !rs.closed_loop)
            snapArcEndpoints(rs.geom, constraints, frame, snap_tol);
        return {rs};
    }

    Segment& s = rs.geom;
    const Circle3 c = segmentCircle(s, frame);

    // Closed full circles (e.g. sphere∩plane) must still split at neighbor hits — otherwise
    // a free-face full loop blocks open chains whose tips lie on that circle (huapingdun+sphere).
    if (rs.closed_loop || std::abs(s.sweep) >= kTwoPi - 1e-6) {
        std::vector<double> angles;
        for (const Vec3& p : constraints) {
            if (std::abs(dist(p, s.center) - s.radius) > snap_tol) continue;
            angles.push_back(wrapTwoPi(frame.angle(s.center, p)));
        }
        std::sort(angles.begin(), angles.end());
        angles.erase(std::unique(angles.begin(), angles.end(),
                                 [](double x, double y) { return std::abs(x - y) <= 1e-9; }),
                     angles.end());
        if (angles.size() < 2) {
            snapArcEndpoints(s, constraints, frame, snap_tol);
            return {rs};
        }
        std::vector<RawSegment> out;
        out.reserve(angles.size());
        for (size_t i = 0; i < angles.size(); ++i) {
            const double t0 = angles[i];
            const double t1 = (i + 1 < angles.size()) ? angles[i + 1] : (angles.front() + kTwoPi);
            RawSegment piece = rs;
            piece.closed_loop = false;
            piece.geom = makeArc(c, t0, t1, frame);
            if (std::abs(piece.geom.sweep) > 1e-9 &&
                dist(piece.geom.start, piece.geom.end) > snap_tol * 1e-3)
                out.push_back(piece);
        }
        for (RawSegment& piece : out) snapArcEndpoints(piece.geom, constraints, frame, snap_tol);
        return out.empty() ? std::vector<RawSegment>{rs} : out;
    }

    std::vector<double> inside;
    inside.reserve(constraints.size());
    for (const Vec3& p : constraints) {
        if (std::abs(dist(p, s.center) - s.radius) > snap_tol) continue;
        const double a = frame.angle(s.center, p);
        if (angleInSweepInterior(a, s.start_angle, s.sweep)) inside.push_back(a);
    }
    std::sort(inside.begin(), inside.end());
    inside.erase(std::unique(inside.begin(), inside.end(),
                             [](double x, double y) { return std::abs(x - y) <= 1e-9; }),
                 inside.end());

    if (inside.empty()) {
        snapArcEndpoints(s, constraints, frame, snap_tol);
        return {rs};
    }

    std::vector<RawSegment> out;
    out.reserve(inside.size() + 1);
    double t0 = s.start_angle;
    const double t_end = s.start_angle + s.sweep;
    for (double a : inside) {
        RawSegment piece = rs;
        piece.geom = makeArc(c, t0, a, frame);
        if (std::abs(piece.geom.sweep) > 1e-9 &&
            dist(piece.geom.start, piece.geom.end) > snap_tol * 1e-3)
            out.push_back(piece);
        t0 = a;
    }
    {
        RawSegment piece = rs;
        piece.geom = makeArc(c, t0, t_end, frame);
        if (std::abs(piece.geom.sweep) > 1e-9 &&
            dist(piece.geom.start, piece.geom.end) > snap_tol * 1e-3)
            out.push_back(piece);
    }
    for (RawSegment& piece : out) snapArcEndpoints(piece.geom, constraints, frame, snap_tol);
    return out;
}

void applyNeighborConstraints(std::vector<RawSegment>& segs, const std::vector<Vec3>& constraints,
                              const SliceFrame& frame, double snap_tol) {
    if (constraints.empty()) return;
    std::vector<RawSegment> split;
    split.reserve(segs.size() + 4);
    for (RawSegment& rs : segs) {
        auto pieces = splitArcAtConstraints(std::move(rs), constraints, frame, snap_tol);
        split.insert(split.end(), std::make_move_iterator(pieces.begin()),
                     std::make_move_iterator(pieces.end()));
    }
    segs = std::move(split);
}

std::vector<Vec3> circleConstraintHits(const Circle3& circle, const Plane& pln,
                                         const SliceFrame& frame, const SliceOptions& opt,
                                         int solid_id) {
    std::vector<Vec3> out;
    if (!opt.plane_faces) return out;
    const double on_circ = std::max(opt.tolerance * 100.0, 0.05);
    for (const FaceRecord* nb : *opt.plane_faces) {
        if (!nb || !nb->face || nb->solid_id != solid_id) continue;
        const SurfaceData surf = nb->face->surface();
        if (surf.kind == SurfaceKind::Other) {
            // NURBS neighbors: use edge∩plane hits that lie on this circle (sphere∩vase).
            const std::vector<Vec3> hits =
                faceBoundaryHits(*nb->face, pln, opt.angular_tolerance, opt.tolerance);
            for (const Vec3& p : hits) {
                if (std::abs(dist(p, circle.center) - circle.radius) <= on_circ) out.push_back(p);
            }
            continue;
        }
        const std::vector<QuadricCurve> curves = intersectAnalyticSurfacePlane(
            surf, pln, opt.angular_tolerance, opt.geom_tolerance);
        for (const QuadricCurve& q : curves) {
            if (q.kind != QuadricKind::Circle || q.degenerate) continue;
            const std::vector<Vec3> hits =
                circleCircleIntersect(circle, q.circle, pln, frame, opt.geom_tolerance);
            out.insert(out.end(), hits.begin(), hits.end());
        }
    }
    return out;
}

std::vector<Vec3> neighborCircleConstraintHits(const AnalyticCurve& self, const FaceRecord& iface,
                                               const Plane& pln, const SliceFrame& frame,
                                               const SliceOptions& opt) {
    if (!opt.plane_faces || self.kind != AnalyticCurve::Circle) return {};
    return circleConstraintHits(self.circle, pln, frame, opt, iface.solid_id);
}

bool sameArcCircle(const Segment& a, const Segment& b, double tol) {
    if (a.type != SegmentType::Arc || b.type != SegmentType::Arc) return false;
    return dist(a.center, b.center) <= tol && std::abs(a.radius - b.radius) <= tol &&
           almostParallel(a.normal, b.normal, std::max(tol, 1e-8));
}

bool coversAngle(const Segment& s, double a_test) {
    constexpr double eps = 1e-8;
    double d = a_test - s.start_angle;
    while (d <= -kPi) d += kTwoPi;
    while (d > kPi) d -= kTwoPi;
    if (s.sweep >= 0.0) return d >= -eps && d <= s.sweep + eps;
    return d <= eps && d >= s.sweep - eps;
}

bool isCoplanarPlanarFace(const SurfaceData& surf, const Plane& pln, double distTol) {
    if (surf.kind != SurfaceKind::Plane) return false;
    const double nd = std::abs(dot(surf.plane.n, pln.n));
    if (nd < 1.0 - 1e-8) return false;
    const Vec3 loc = surf.plane.n * surf.plane.d;
    return std::abs(signedPlaneDist(pln, loc)) <= distTol;
}

bool curveToSegment(const IEdge& e, const SliceFrame& frame, Segment& out) {
    const CurveData c = e.curve();
    if (c.degenerated) return false;
    switch (c.kind) {
    case CurveKind::Line:
        out = makeLineSeg(e.eval(c.first), e.eval(c.last));
        break;
    case CurveKind::Circle:
        out = makeArc(c.circle, c.first, c.last, frame);
        break;
    case CurveKind::Ellipse:
        out = makeEllipseSeg(c.ellipse, c.first, c.last, frame);
        break;
    default:
        out = makeLineSeg(e.eval(c.first), e.eval(c.last));
        break;
    }
    if (c.reversed) reverseSegment(out);
    return true;
}

std::vector<RawSegment> dumpCoplanarFace(const FaceRecord& iface, const SliceFrame& frame,
                                         const SliceOptions& opt) {
    std::vector<RawSegment> segs;
    if (!iface.face) return segs;
    for (const auto& wire : iface.face->wires()) {
        for (const auto& e : wire) {
            if (!e) continue;
            Segment s;
            if (!curveToSegment(*e, frame, s)) continue;
            if (dist(s.start, s.end) <= opt.geom_tolerance && s.type == SegmentType::Line) continue;
            RawSegment rs;
            rs.geom = s;
            rs.solid_id = iface.solid_id;
            rs.shell_id = iface.shell_id;
            rs.face_id = iface.face_id;
            rs.chain_idx = 0;
            rs.geom.face_id = iface.face_id;
            rs.geom.chain_idx = 0;
            rs.coplanar = true;
            segs.push_back(rs);
        }
    }
    return segs;
}

RawSegment degenSeg(const FaceRecord& iface, const std::string& event, const Vec3& p) {
    RawSegment rs;
    rs.solid_id = iface.solid_id;
    rs.shell_id = iface.shell_id;
    rs.face_id = iface.face_id;
    rs.chain_idx = 0;
    rs.geom.face_id = iface.face_id;
    rs.geom.chain_idx = 0;
    rs.degenerate = true;
    rs.degen_event = event;
    rs.degen_point = p;
    return rs;
}

std::vector<RawSegment> trimCurveToFace(const AnalyticCurve& curve, const FaceRecord& iface,
                                        const std::vector<Vec3>& hits, const SliceFrame& frame,
                                        const SliceOptions& opt) {
    std::vector<RawSegment> segs;
    if (!iface.face) return segs;
    const double matchTol = std::max(opt.tolerance, 10.0 * opt.geom_tolerance);

    std::vector<double> ts;
    // Looser candidates for tip extension when the analytic generator classifies Out
    // a fraction of a mm before the shared corner (ROBOT z≈373 face 889).
    struct HitT {
        double t;
        Vec3 p;
        double derr;
    };
    std::vector<HitT> pull_hits;
    ts.reserve(hits.size());
    const double snapLoose = std::max({matchTol, 0.5, 50.0 * opt.tolerance});
    // Oblique-cyl ellipses: boundary corners can sit slightly off the analytic
    // plane∩quadric; under-snapping drops real trim ends and leaves a 2-hit sliver
    // (ROBOT_7 layers 4–20). Circles stay on matchTol (exact).
    const double snapTrim =
        (curve.kind == AnalyticCurve::Ellipse)
            ? std::max(matchTol, std::min(snapLoose, 0.2))
            : matchTol;
    for (const Vec3& p : hits) {
        double t = 0;
        const double derr = projectToCurve(curve, p, t);
        if (curve.periodic) t = wrapTwoPi(t);
        if (derr <= snapLoose) pull_hits.push_back({t, p, derr});
        // Keep primary trim hits tight for lines/circles so //cylinder generators
        // do not pick up nearby unrelated edge hits (ROBOT_4 z≈209).
        if (derr <= snapTrim) ts.push_back(t);
    }
    std::sort(ts.begin(), ts.end());
    const double merge = 1e-9;
    std::vector<double> uniq;
    for (double t : ts) {
        if (uniq.empty() || std::abs(t - uniq.back()) > merge) {
            if (!uniq.empty() && curve.periodic &&
                std::abs((t - uniq.front()) - curve.period) < merge) {
                continue;
            }
            uniq.push_back(t);
        }
    }

    constexpr double kMinOpenArc = 0.05;  // radians; skip trim slivers

    // Same on-face test as circles: 3D In|On then UV wire (classifyInOrOnFace).
    // For slightly-off analytic ellipse samples, snap to the surface first when the
    // projection stays within a tight band — do not use bare 3D classify alone
    // (that marks infinite-cylinder points In outside the trim).
    auto sampleOnFace = [&](double t) -> PointClass {
        const Vec3 p = evalCurve(curve, t);
        const double ct = opt.tolerance;
        double u = 0, v = 0;
        if (iface.face->invertUV(p, u, v, std::max(ct * 10.0, 1e-3))) {
            const Vec3 q = iface.face->evalUV(u, v);
            if (dist(p, q) <= std::max(0.05, 100.0 * ct)) {
                if (uvInOrOnPeriodic(*iface.face, u, v, ct)) return PointClass::In;
                return iface.face->classifyUV(u, v, ct);
            }
        }
        return classifyInOrOnFace(*iface.face, p, ct) ? PointClass::In : PointClass::Out;
    };

    auto intervalOnTrimFace = [&](double t0, double t1) -> bool {
        int in_n = 0, on_n = 0, out_n = 0;
        constexpr int kSamp = 5;
        for (int k = 1; k <= kSamp; ++k) {
            const double t = t0 + (t1 - t0) * (static_cast<double>(k) / (kSamp + 1));
            const PointClass st = sampleOnFace(t);
            if (st == PointClass::Out) {
                ++out_n;
                continue;
            }
            if (st == PointClass::In)
                ++in_n;
            else
                ++on_n;
        }
        // Short //cylinder generator closeouts (ROBOT z≈290.6 F206 ~1 mm): seam UV can
        // flick a single sample Out while the span is still the on-face ribbon.
        // Only forgive tip-adjacent Out — an interior Out is a real hole/notch crossing
        // (ROBOT z≈506 F209 ~0.7 mm mouth must stay split).
        if (curve.kind == AnalyticCurve::Line) {
            const double len = arcLength(curve, t0, t1);
            if (len > opt.geom_tolerance && len <= 8.0 && (in_n + on_n) >= 2 && out_n <= 1) {
                bool out_interior = false;
                for (int k = 1; k <= kSamp; ++k) {
                    const double alpha = static_cast<double>(k) / (kSamp + 1);
                    if (alpha <= 0.2 || alpha >= 0.8) continue;
                    const double t = t0 + (t1 - t0) * alpha;
                    if (sampleOnFace(t) == PointClass::Out) {
                        out_interior = true;
                        break;
                    }
                }
                if (!out_interior) return true;
            }
        }
        if (out_n > 0) return false;
        return (in_n + on_n) > 0;
    };

    auto tryEmit = [&](double t0, double t1, bool closed) -> bool {
        if (!closed && arcLength(curve, t0, t1) <= opt.geom_tolerance) return false;
        if (!closed && curve.periodic) {
            double span = t1 - t0;
            while (span <= -curve.period) span += curve.period;
            while (span > curve.period) span -= curve.period;
            // Preserve major spans (circle ROBOT z≈70 / ellipse ROBOT_7). Only drop
            // true microscopic slivers — do not fold |span|>π before this check.
            const double minOpen =
                (curve.kind == AnalyticCurve::Ellipse) ? 1e-3 : kMinOpenArc;
            if (std::abs(span) < minOpen) return false;
        }
        if (!intervalOnTrimFace(t0, t1)) {
            if (opt.constraint_audit &&
                (curve.kind == AnalyticCurve::Circle || curve.kind == AnalyticCurve::Ellipse)) {
                const Vec3 mid = evalCurve(curve, 0.5 * (t0 + t1));
                std::ostringstream os;
                os << "arc_trim_reject face=" << iface.face_id
                   << " kind=" << (curve.kind == AnalyticCurve::Ellipse ? "ell" : "cir")
                   << " t0=" << t0 << " t1=" << t1 << " mid=(" << mid.x << "," << mid.y
                   << ") class=" << static_cast<int>(sampleOnFace(0.5 * (t0 + t1)));
                opt.constraint_audit->push_back(os.str());
            }
            return false;
        }
        RawSegment rs;
        rs.geom = emitInterval(curve, t0, t1, frame);
        rs.solid_id = iface.solid_id;
        rs.shell_id = iface.shell_id;
        rs.face_id = iface.face_id;
        rs.chain_idx = 0;
        rs.geom.face_id = iface.face_id;
        rs.geom.chain_idx = 0;
        rs.closed_loop = closed;
        segs.push_back(rs);
        if (opt.constraint_audit &&
            (curve.kind == AnalyticCurve::Circle || curve.kind == AnalyticCurve::Ellipse ||
             curve.kind == AnalyticCurve::Line)) {
            std::ostringstream os;
            os << "arc_trim_keep face=" << iface.face_id
               << " kind="
               << (curve.kind == AnalyticCurve::Ellipse
                       ? "ell"
                       : (curve.kind == AnalyticCurve::Circle ? "cir" : "line"))
               << " t0=" << t0 << " t1=" << t1 << " sweep=" << (t1 - t0)
               << " len=" << arcLength(curve, t0, t1)
               << " start=(" << rs.geom.start.x << "," << rs.geom.start.y << ")"
               << " end=(" << rs.geom.end.x << "," << rs.geom.end.y << ")";
            opt.constraint_audit->push_back(os.str());
        }
        return true;
    };

    auto finalize = [&]() -> std::vector<RawSegment> {
        // Extend open tips along the curve onto nearby edge∩plane hits so adjacent
        // faces share vertices. Exact on-curve hits use a tight band; fitting-edge
        // partners (free/NURBS wall) may sit ~0.01–0.05 mm off the generator
        // (ROBOT_4 z≈278 F257↔cyl) — allow that lateral band but always snap the
        // endpoint onto the analytic curve (evalCurve), never the raw edge hit.
        // Loose 0.5 mm lateral pulls still break ROBOT_4 z≈209.
        // ROBOT z≈373 face 889: tip y=-41.86 vs corner y=-42.72 (derr≈2e-4).
        if (!pull_hits.empty() && curve.kind == AnalyticCurve::Line) {
            const double max_pull = 2.0;
            const double fit_lateral = std::max(0.05, 500.0 * opt.tolerance);
            for (RawSegment& rs : segs) {
                if (rs.closed_loop) continue;
                double t0 = 0, t1 = 0;
                projectToCurve(curve, rs.geom.start, t0);
                projectToCurve(curve, rs.geom.end, t1);
                const bool start_is_lo = t0 <= t1;
                double tLo = std::min(t0, t1), tHi = std::max(t0, t1);
                for (const HitT& h : pull_hits) {
                    if (h.derr > fit_lateral) continue;
                    // Refuse pulls across a trimmed hole/notch: the far hit may be the
                    // other mouth of the same generator (ROBOT z≈506 F209, gap≈0.7 mm
                    // < max_pull). Midpoint must stay In|On the face.
                    auto pullOnFace = [&](double from_t, double to_t) {
                        const double mid = 0.5 * (from_t + to_t);
                        return classifyInOrOnFace(*iface.face, evalCurve(curve, mid),
                                                  opt.tolerance);
                    };
                    const Vec3 on_p = evalCurve(curve, h.t);
                    if (h.t < tLo - 1e-9 && (tLo - h.t) <= max_pull) {
                        if (!pullOnFace(h.t, tLo)) continue;
                        if (start_is_lo) rs.geom.start = on_p;
                        else rs.geom.end = on_p;
                        tLo = h.t;
                    } else if (h.t > tHi + 1e-9 && (h.t - tHi) <= max_pull) {
                        if (!pullOnFace(tHi, h.t)) continue;
                        if (start_is_lo) rs.geom.end = on_p;
                        else rs.geom.start = on_p;
                        tHi = h.t;
                    }
                }
            }
        }
        return segs;
    };

    if (uniq.empty()) {
        if (curve.periodic) {
            const Vec3 sample = evalCurve(curve, 0.0);
            if (classifyInOrOnFace(*iface.face, sample, opt.tolerance)) {
                tryEmit(0.0, curve.period, true);
            }
            return finalize();
        }
        // Non-periodic (generator lines from // cylinder/cone): no edge hits on this
        // plane — classify along the face bbox instead of dropping the fillet.
        if (curve.kind == AnalyticCurve::Line) {
            std::vector<std::pair<double, double>> ivals;
            classifyLineIntervals(curve, *iface.face, opt.tolerance, {}, ivals);
            for (const auto& iv : ivals) tryEmit(iv.first, iv.second, false);
        }
        return finalize();
    }

    if (curve.periodic) {
        if (opt.constraint_audit &&
            (curve.kind == AnalyticCurve::Circle || curve.kind == AnalyticCurve::Ellipse)) {
            std::ostringstream os;
            os << "arc_trim_hits face=" << iface.face_id
               << " kind=" << (curve.kind == AnalyticCurve::Ellipse ? "ell" : "cir")
               << " n=" << uniq.size();
            for (double t : uniq) os << " " << t;
            opt.constraint_audit->push_back(os.str());
        }
        if (uniq.size() == 1) {
            tryEmit(uniq[0], uniq[0] + curve.period, true);
            return finalize();
        }
        if (uniq.size() == 2) {
            const double t0 = uniq[0], t1 = uniq[1];
            double gapInt = t1 - t0;
            if (gapInt <= 0.0) gapInt += curve.period;
            const double wrapSpan = curve.period - gapInt;
            auto sideScore = [&](double a, double b) {
                // Same scoring for circle and ellipse (ROBOT z≈70 circle trim).
                int score = 0;
                for (int k = 1; k <= 3; ++k) {
                    const double t = a + (b - a) * (static_cast<double>(k) / 4.0);
                    const PointClass st = sampleOnFace(t);
                    if (st == PointClass::In) score += 2;
                    else if (st == PointClass::On) score += 0;
                    else score -= 3;
                }
                return score;
            };
            const int scoreInt = sideScore(t0, t1);
            const int scoreWrap = sideScore(t1, t1 + wrapSpan);
            const bool inInt = scoreInt > 0;
            const bool inWrap = scoreWrap > 0;
            if (opt.constraint_audit) {
                std::ostringstream os;
                os << "arc_trim_2hit face=" << iface.face_id
                   << " kind=" << (curve.kind == AnalyticCurve::Ellipse ? "ell" : "cir")
                   << " gapInt=" << gapInt << " wrapSpan=" << wrapSpan
                   << " scoreInt=" << scoreInt << " scoreWrap=" << scoreWrap
                   << " inInt=" << inInt << " inWrap=" << inWrap;
                if (inInt && inWrap)
                    os << " keep=" << (scoreInt > scoreWrap
                                           ? "int"
                                           : (scoreWrap > scoreInt
                                                  ? "wrap"
                                                  : (gapInt >= wrapSpan ? "long_int" : "long_wrap")));
                else
                    os << " keep=" << (inInt ? "int" : (inWrap ? "wrap" : "none"));
                opt.constraint_audit->push_back(os.str());
            }
            // Identical selection to circle faces: stronger In score, else longer span.
            if (inInt && inWrap) {
                if (scoreInt > scoreWrap || (scoreInt == scoreWrap && gapInt >= wrapSpan)) {
                    if (!tryEmit(t0, t1, false)) tryEmit(t1, t1 + wrapSpan, false);
                } else {
                    if (!tryEmit(t1, t1 + wrapSpan, false)) tryEmit(t0, t1, false);
                }
            } else if (inInt) {
                tryEmit(t0, t1, false);
            } else if (inWrap) {
                tryEmit(t1, t1 + wrapSpan, false);
            }
            return finalize();
        }
        for (size_t i = 0; i + 1 < uniq.size(); ++i) tryEmit(uniq[i], uniq[i + 1], false);
        // Always emit the periodic wrap last→first (forward in param). When that wrap is
        // the long lobe (>π), the old branch substituted the complementary short span —
        // already covered by consecutive tryEmit above — and skipped the real on-face
        // ribbon (ROBOT z≈545 F125 hole: hits at π/2, mouth, mouth; long right-side arc
        // never tried). intervalOnTrimFace still rejects Out wraps.
        {
            double tWrap0 = uniq.back();
            double tWrap1 = uniq.front();
            if (tWrap1 <= tWrap0) tWrap1 += curve.period;
            tryEmit(tWrap0, tWrap1, false);
        }
    } else if (curve.kind == AnalyticCurve::Line) {
        // Promote slightly-off shared-edge hits (NURBS fillet tips) into the trim set.
        // matchTol alone drops F165's tip on F206 (derr≈0.002) so consecutive emit only
        // keeps the top stub (y=32→30.56) and never sees the ~1 mm planar closeout
        // (ROBOT z≈290.6 open gap). Cap promotion so //cyl orphans (ROBOT_4 z≈209) stay out.
        std::vector<double> line_hits = uniq;
        {
            const double promote = std::max({0.35, 50.0 * opt.tolerance, matchTol});
            for (const HitT& h : pull_hits) {
                if (h.derr > promote) continue;
                bool dup = false;
                for (double t : line_hits) {
                    if (std::abs(t - h.t) <= std::max(opt.tolerance, 1e-6)) {
                        dup = true;
                        break;
                    }
                }
                if (!dup) line_hits.push_back(h.t);
            }
            std::sort(line_hits.begin(), line_hits.end());
        }
        // Prefer consecutive edge∩plane hits. Bbox classifyLineIntervals often extends
        // //cylinder generators past the trimmed wire (ROBOT_4 z≈202 orphans).
        if (line_hits.size() >= 2) {
            for (size_t i = 0; i + 1 < line_hits.size(); ++i) tryEmit(line_hits[i], line_hits[i + 1], false);
        }
        // Always classify for leftover In ribbons: a successful long-neighbor stub must
        // not suppress a disjoint short closeout on the same generator (F206 right).
        {
            std::vector<std::pair<double, double>> ivals;
            classifyLineIntervals(curve, *iface.face, opt.tolerance, line_hits, ivals);
            for (const auto& iv : ivals) {
                // Skip intervals already covered by a kept segment.
                bool covered = false;
                for (const RawSegment& s : segs) {
                    if (s.degenerate) continue;
                    double ta = 0, tb = 0;
                    projectToCurve(curve, s.geom.start, ta);
                    projectToCurve(curve, s.geom.end, tb);
                    const double lo = std::min(ta, tb), hi = std::max(ta, tb);
                    if (iv.first >= lo - opt.tolerance && iv.second <= hi + opt.tolerance) {
                        covered = true;
                        break;
                    }
                }
                if (!covered) tryEmit(iv.first, iv.second, false);
            }
        }
        if (opt.constraint_audit && segs.empty()) {
            std::ostringstream os;
            os << "line_trim_empty face=" << iface.face_id << " hits=" << uniq.size()
               << " pull=" << pull_hits.size() << " line_hits=" << line_hits.size();
            opt.constraint_audit->push_back(os.str());
        }
    } else {
        for (size_t i = 0; i + 1 < uniq.size(); ++i) tryEmit(uniq[i], uniq[i + 1], false);
    }

    return finalize();
}

}  // namespace

void dedupeAnalyticArcs(std::vector<RawSegment>& segs, double tol) {
    const double kMinSweep = 1e-6;
    const double match = std::max(tol, 1e-6);

    auto nearlySameArc = [&](const Segment& a, const Segment& b) {
        const bool fwd =
            dist(a.start, b.start) <= match && dist(a.end, b.end) <= match;
        const bool rev =
            dist(a.start, b.end) <= match && dist(a.end, b.start) <= match;
        return fwd || rev;
    };

    auto arcCovers = [&](const Segment& outer, const Segment& inner) {
        if (std::abs(inner.sweep) < kMinSweep) return true;
        if (std::abs(outer.sweep) >= kTwoPi - 1e-6) return true;
        const double a0 = inner.start_angle;
        const double a1 = inner.start_angle + 0.5 * inner.sweep;
        const double a2 = inner.start_angle + inner.sweep;
        return coversAngle(outer, a0) && coversAngle(outer, a1) && coversAngle(outer, a2);
    };

    auto preferKeep = [&](const RawSegment& a, const RawSegment& b) {
        if (a.coplanar != b.coplanar) return a.coplanar;
        if (a.closed_loop != b.closed_loop) return a.closed_loop;
        if (std::abs(a.geom.sweep) != std::abs(b.geom.sweep))
            return std::abs(a.geom.sweep) >= std::abs(b.geom.sweep);
        return a.face_id <= b.face_id;
    };

    for (size_t i = 0; i < segs.size(); ++i) {
        if (segs[i].degenerate || segs[i].geom.type != SegmentType::Arc) continue;
        if (std::abs(segs[i].geom.sweep) < kMinSweep) {
            segs[i].degenerate = true;
            continue;
        }
        for (size_t j = i + 1; j < segs.size(); ++j) {
            if (segs[j].degenerate || segs[j].geom.type != SegmentType::Arc) continue;
            if (segs[i].solid_id != segs[j].solid_id) continue;
            if (!sameArcCircle(segs[i].geom, segs[j].geom, match)) continue;

            // Exact twin drop (existing).
            if (nearlySameArc(segs[i].geom, segs[j].geom)) {
                if (preferKeep(segs[i], segs[j])) segs[j].degenerate = true;
                else {
                    segs[i].degenerate = true;
                    break;
                }
                continue;
            }

            // Partial / subset arc: coversAngle — keep superset (or longer).
            if (arcCovers(segs[i].geom, segs[j].geom)) {
                segs[i].coplanar = segs[i].coplanar || segs[j].coplanar;
                segs[j].degenerate = true;
                continue;
            }
            if (arcCovers(segs[j].geom, segs[i].geom)) {
                segs[j].coplanar = segs[j].coplanar || segs[i].coplanar;
                segs[i].degenerate = true;
                break;
            }
        }
    }

    segs.erase(std::remove_if(segs.begin(), segs.end(),
                              [kMinSweep](const RawSegment& rs) {
                                  if (rs.degenerate) return true;
                                  return rs.geom.type == SegmentType::Arc &&
                                         std::abs(rs.geom.sweep) < kMinSweep;
                              }),
               segs.end());
}

std::vector<RawSegment> intersectFaceWithPlane(const FaceRecord& iface, const Plane& pln,
                                               const SliceFrame& frame, const SliceOptions& opt) {
    if (!iface.face) return {};
    const SurfaceData surf = iface.face->surface();
    if (isCoplanarPlanarFace(surf, pln, opt.tolerance)) {
        return dumpCoplanarFace(iface, frame, opt);
    }

    std::vector<Vec3> boundary =
        faceBoundaryHits(*iface.face, pln, opt.angular_tolerance, opt.tolerance);
    if (opt.constraint_audit) {
        auditFaceBoundaryHits(iface.face_id, *iface.face, pln, opt.angular_tolerance, opt.tolerance,
                              opt.constraint_audit);
    }

    if (surf.kind == SurfaceKind::Other) {
        snapBoundaryHitsToAnalyticNeighbors(boundary, iface, pln, opt);
        if (opt.nurbs_method == NurbsMethod::UvMarch) {
            return intersectNurbsFaceWithPlane(iface, pln, frame, opt, boundary);
        }
        return intersectNurbsFaceWithPlaneUvMatch(iface, pln, frame, opt, boundary);
    }

    const std::vector<QuadricCurve> qcurves =
        intersectAnalyticSurfacePlane(surf, pln, opt.angular_tolerance, opt.geom_tolerance);

    std::vector<RawSegment> segs;
    const double matchTol = std::max(opt.tolerance, opt.geom_tolerance);
    (void)matchTol;

    for (const QuadricCurve& q : qcurves) {
        if (q.kind == QuadricKind::Coplanar) {
            auto part = dumpCoplanarFace(iface, frame, opt);
            segs.insert(segs.end(), part.begin(), part.end());
            continue;
        }
        if (q.kind == QuadricKind::Empty) continue;
        if (q.kind == QuadricKind::Point || q.degenerate) {
            std::string ev = "tangent_point";
            if (q.kind == QuadricKind::Line) ev = "tangent_line";
            else if (q.kind == QuadricKind::Circle) ev = "degenerate_ring";
            segs.push_back(degenSeg(iface, ev, q.point));
            continue;
        }

        AnalyticCurve c;
        if (q.kind == QuadricKind::Line) {
            c.kind = AnalyticCurve::Line;
            c.line = q.line;
            c.periodic = false;
        } else if (q.kind == QuadricKind::Circle) {
            c.kind = AnalyticCurve::Circle;
            c.circle = q.circle;
            c.periodic = true;
            c.period = kTwoPi;
        } else if (q.kind == QuadricKind::Ellipse) {
            c.kind = AnalyticCurve::Ellipse;
            c.ellipse = q.ellipse;
            c.periodic = true;
            c.period = kTwoPi;
        } else {
            continue;
        }
        auto part = trimCurveToFace(c, iface, boundary, frame, opt);
        if (opt.plane_faces && c.kind == AnalyticCurve::Circle) {
            const std::vector<Vec3> constraints =
                neighborCircleConstraintHits(c, iface, pln, frame, opt);
            if (!constraints.empty()) {
                const double snap_tol = std::max(opt.tolerance * 1000.0, 0.2);
                applyNeighborConstraints(part, constraints, frame, snap_tol);
            }
        }
        segs.insert(segs.end(), part.begin(), part.end());
    }

    // Analytic plane∩cone/cyl parallel to the axis can be a hyperbola (unsupported) or
    // a false tangent Point — OCC section still has real fillet curves. Fall back to
    // UVMatch so those faces still contribute.
    const bool any_real =
        std::any_of(segs.begin(), segs.end(), [](const RawSegment& s) { return !s.degenerate; });
    if (!any_real && (surf.kind == SurfaceKind::Cone || surf.kind == SurfaceKind::Cylinder ||
                      surf.kind == SurfaceKind::Torus)) {
        snapBoundaryHitsToAnalyticNeighbors(boundary, iface, pln, opt);
        if (opt.nurbs_method == NurbsMethod::UvMarch) {
            return intersectNurbsFaceWithPlane(iface, pln, frame, opt, boundary);
        }
        return intersectNurbsFaceWithPlaneUvMatch(iface, pln, frame, opt, boundary);
    }

    if (segs.empty() && opt.constraint_audit && surf.kind != SurfaceKind::Other) {
        std::ostringstream os;
        os << "analytic_empty face=" << iface.face_id
           << " surf=" << static_cast<int>(surf.kind) << " boundary=" << boundary.size()
           << " qn=" << qcurves.size();
        for (const QuadricCurve& q : qcurves) {
            os << " k=" << static_cast<int>(q.kind) << (q.degenerate ? "*" : "");
        }
        opt.constraint_audit->push_back(os.str());
    }
    return segs;
}

}  // namespace brepslicer
