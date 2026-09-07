#include <intersect/AnalyticIntersect.h>
#include <intersect/BSplineFit.h>
#include <intersect/QuadricPlane.h>
#include <intersect/UvMarch.h>
#include <intersect/UvMatch.h>
#include <geom/GeomUtil.h>

#include <algorithm>
#include <cmath>
#include <fstream>
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

    double dt = t1 - t0;
    const double alt = dt > 0.0 ? dt - kTwoPi : dt + kTwoPi;
    const Vec3 refMid = c.at(0.5 * (t0 + t1));
    auto midErr = [&](double d) { return length(c.at(t0 + 0.5 * d) - refMid); };
    const double errDt = midErr(dt);
    const double errAlt = midErr(alt);
    constexpr double midTol = 1e-10;
    if (errAlt + midTol < errDt) {
        dt = alt;
    } else if (errAlt <= errDt + midTol && std::abs(dt) > kPi && std::abs(alt) <= kPi) {
        dt = alt;
    }
    if (std::abs(dt) > kPi && std::abs(dt) < kTwoPi - 1e-9) {
        dt = dt > 0.0 ? dt - kTwoPi : dt + kTwoPi;
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
    s.sweep = t1 - t0;
    return s;
}

bool pointOnPlane(const Vec3& p, const Plane& pln, double tol) {
    return std::abs(signedPlaneDist(pln, p)) <= tol;
}

bool classifyInOrOnFace(const IFace& face, const Vec3& p, double tol) {
    const PointClass st = face.classify(p, tol);
    return st == PointClass::In || st == PointClass::On;
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
    return pts;
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

Vec3 evalCurve(const AnalyticCurve& c, double t) {
    if (c.kind == AnalyticCurve::Line) return c.line.at(t);
    if (c.kind == AnalyticCurve::Circle) return c.circle.at(t);
    return c.ellipse.at(t);
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
    if (is_start) {
        const double a_end = s.start_angle + s.sweep;
        s.start_angle = a_hit;
        s.sweep = a_end - a_hit;
        if (s.sweep > kPi) s.sweep -= kTwoPi;
        if (s.sweep < -kPi) s.sweep += kTwoPi;
    } else {
        double sw = a_hit - s.start_angle;
        if (sw > kPi) sw -= kTwoPi;
        if (sw < -kPi) sw += kTwoPi;
        s.sweep = sw;
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
    if (rs.geom.type != SegmentType::Arc || rs.closed_loop || constraints.empty()) {
        if (rs.geom.type == SegmentType::Arc && !rs.closed_loop)
            snapArcEndpoints(rs.geom, constraints, frame, snap_tol);
        return {rs};
    }

    Segment& s = rs.geom;
    const Circle3 c = segmentCircle(s, frame);
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
    for (const FaceRecord* nb : *opt.plane_faces) {
        if (!nb || !nb->face || nb->solid_id != solid_id) continue;
        const SurfaceData surf = nb->face->surface();
        if (surf.kind == SurfaceKind::Other) continue;
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
    ts.reserve(hits.size());
    for (const Vec3& p : hits) {
        double t = 0;
        if (snapToCurve(curve, p, matchTol, t)) {
            if (curve.periodic) t = wrapTwoPi(t);
            ts.push_back(t);
        }
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

    auto tryEmit = [&](double t0, double t1, bool closed) -> bool {
        if (!closed && arcLength(curve, t0, t1) <= opt.geom_tolerance) return false;
        if (!closed && curve.periodic) {
            double span = t1 - t0;
            while (span <= -curve.period) span += curve.period;
            while (span > curve.period) span -= curve.period;
            if (std::abs(span) > curve.period * 0.5) span = span > 0 ? span - curve.period : span + curve.period;
            if (std::abs(span) < kMinOpenArc) return false;
        }
        const Vec3 mid = evalCurve(curve, 0.5 * (t0 + t1));
        if (!classifyInOrOnFace(*iface.face, mid, opt.tolerance)) {
            if (curve.kind == AnalyticCurve::Circle && opt.constraint_audit) {
                std::ostringstream os;
                os << "arc_trim_reject face=" << iface.face_id << " t0=" << t0 << " t1=" << t1
                   << " mid=(" << mid.x << "," << mid.y << "," << mid.z << ")"
                   << " class=" << static_cast<int>(iface.face->classify(mid, opt.tolerance));
                opt.constraint_audit->push_back(os.str());
            }
            return false;
        }
        RawSegment rs;
        rs.geom = emitInterval(curve, t0, t1, frame);
        rs.solid_id = iface.solid_id;
        rs.shell_id = iface.shell_id;
        rs.face_id = iface.face_id;
        rs.closed_loop = closed;
        segs.push_back(rs);
        if (curve.kind == AnalyticCurve::Circle && opt.constraint_audit) {
            std::ostringstream os;
            os << "arc_trim_keep face=" << iface.face_id << " t0=" << t0 << " t1=" << t1
               << " sweep=" << (t1 - t0) << " start=(" << rs.geom.start.x << "," << rs.geom.start.y
               << ") end=(" << rs.geom.end.x << "," << rs.geom.end.y << ")";
            opt.constraint_audit->push_back(os.str());
        }
        return true;
    };

    if (uniq.empty()) {
        if (curve.periodic) {
            const Vec3 sample = evalCurve(curve, 0.0);
            if (classifyInOrOnFace(*iface.face, sample, opt.tolerance)) {
                tryEmit(0.0, curve.period, true);
            }
        }
        return segs;
    }

    if (curve.periodic) {
        if (curve.kind == AnalyticCurve::Circle && opt.constraint_audit) {
            std::ostringstream os;
            os << "arc_trim_hits face=" << iface.face_id << " n=" << uniq.size();
            for (double t : uniq) os << " " << t;
            opt.constraint_audit->push_back(os.str());
        }
        if (uniq.size() == 1) {
            tryEmit(uniq[0], uniq[0] + curve.period, true);
            return segs;
        }
        if (uniq.size() == 2) {
            const double t0 = uniq[0], t1 = uniq[1];
            double gapInt = t1 - t0;
            if (gapInt <= 0.0) gapInt += curve.period;
            auto midInside = [&](double a, double b) {
                const Vec3 mid = evalCurve(curve, 0.5 * (a + b));
                return classifyInOrOnFace(*iface.face, mid, opt.tolerance);
            };
            const bool inInt = midInside(t0, t1);
            const bool inWrap = midInside(t1, t1 + (curve.period - gapInt));
            if (opt.constraint_audit) {
                std::ostringstream os;
                os << "arc_trim_2hit face=" << iface.face_id << " gapInt=" << gapInt
                   << " inInt=" << inInt << " inWrap=" << inWrap
                   << " keep=" << ((inInt && inWrap)
                                       ? (gapInt <= 0.5 * curve.period ? "short_int" : "short_wrap")
                                       : (inInt ? "int" : (inWrap ? "wrap" : "none")));
                opt.constraint_audit->push_back(os.str());
            }
            if (inInt && inWrap) {
                if (gapInt <= 0.5 * curve.period) {
                    tryEmit(t0, t1, false);
                } else {
                    tryEmit(t1, t1 + (curve.period - gapInt), false);
                }
            } else if (inInt) {
                tryEmit(t0, t1, false);
            } else if (inWrap) {
                tryEmit(t1, t1 + (curve.period - gapInt), false);
            }
            return segs;
        }
        for (size_t i = 0; i + 1 < uniq.size(); ++i) tryEmit(uniq[i], uniq[i + 1], false);
        double gapFwd = uniq.front() - uniq.back();
        if (gapFwd <= 0.0) gapFwd += curve.period;
        if (gapFwd > 0.5 * curve.period) {
            tryEmit(uniq.front(), uniq.front() + (curve.period - gapFwd), false);
        } else {
            tryEmit(uniq.back(), uniq.front(), false);
        }
    } else {
        for (size_t i = 0; i + 1 < uniq.size(); ++i) tryEmit(uniq[i], uniq[i + 1], false);
    }
    return segs;
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
            if (nearlySameArc(segs[i].geom, segs[j].geom)) segs[j].degenerate = true;
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

    const std::vector<Vec3> boundary =
        faceBoundaryHits(*iface.face, pln, opt.angular_tolerance, opt.tolerance);

    if (surf.kind == SurfaceKind::Other) {
        if (opt.nurbs_method == NurbsMethod::UvMarch) {
            return intersectNurbsFaceWithPlane(iface, pln, frame, opt, boundary);
        }
        return intersectNurbsFaceWithPlaneUvMatch(iface, pln, frame, opt, boundary);
    }

    const std::vector<QuadricCurve> qcurves =
        intersectAnalyticSurfacePlane(surf, pln, opt.angular_tolerance, opt.geom_tolerance);

    std::vector<RawSegment> segs;
    const double matchTol = std::max(opt.tolerance, opt.geom_tolerance);

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
    return segs;
}

}  // namespace brepslicer
