#include <intersect/AnalyticIntersect.h>
#include <intersect/BSplineFit.h>
#include <intersect/QuadricPlane.h>
#include <intersect/UvMarch.h>
#include <geom/GeomUtil.h>

#include <algorithm>
#include <cmath>
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
    const Vec3 mid = c.at(0.5 * (t0 + t1));
    Segment s;
    s.type = SegmentType::Arc;
    s.start = p0;
    s.end = p1;
    s.center = c.center;
    s.normal = frame.n;
    s.radius = c.radius;
    s.start_angle = frame.angle(c.center, p0);
    const double mag = std::abs(t1 - t0);
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

    // A single midpoint is not enough for trimmed faces (e.g. booleaned spheres):
    // the complete surface circle may pass one In sample while most of the arc is Out.
    auto intervalOnFace = [&](double t0, double t1, bool closed) -> bool {
        const int nSamp = closed ? 16 : 7;
        for (int i = 1; i <= nSamp; ++i) {
            const double t = t0 + (t1 - t0) * (double(i) / double(nSamp + 1));
            if (!classifyInOrOnFace(*iface.face, evalCurve(curve, t), opt.tolerance))
                return false;
        }
        return true;
    };

    auto pushInterval = [&](double t0, double t1, bool closed) {
        if (!closed && arcLength(curve, t0, t1) <= opt.geom_tolerance) return;
        if (!intervalOnFace(t0, t1, closed)) return;
        RawSegment rs;
        rs.geom = emitInterval(curve, t0, t1, frame);
        rs.solid_id = iface.solid_id;
        rs.shell_id = iface.shell_id;
        rs.face_id = iface.face_id;
        rs.closed_loop = closed;
        segs.push_back(rs);
    };

    // When boundary hits are missing/incomplete, classify densely on the analytic curve
    // and keep only contiguous In/On runs (respects OCC face trimming).
    auto emitBySampling = [&]() {
        if (!curve.periodic) return;
        constexpr int N = 64;
        std::vector<char> on(N, 0);
        int onCount = 0;
        for (int i = 0; i < N; ++i) {
            const double t = curve.period * (double(i) / double(N));
            on[i] = classifyInOrOnFace(*iface.face, evalCurve(curve, t), opt.tolerance) ? 1 : 0;
            onCount += on[i];
        }
        if (onCount == 0) return;
        if (onCount == N) {
            pushInterval(0.0, curve.period, true);
            return;
        }
        auto at = [&](int i) -> bool { return on[(i % N + N) % N] != 0; };
        int start = -1;
        for (int i = 0; i < N; ++i) {
            if (at(i) && !at(i - 1)) {
                start = i;
                break;
            }
        }
        if (start < 0) return;
        int i = start;
        do {
            const int run0 = i;
            int j = (i + 1) % N;
            while (j != run0 && at(j))
                j = (j + 1) % N;
            double t0 = curve.period * (double(run0) / double(N));
            double t1 = curve.period * (double(j) / double(N));
            if (j <= run0)
                t1 += curve.period;
            pushInterval(t0, t1, false);
            i = j;
            while (i != start && !at(i))
                i = (i + 1) % N;
        } while (i != start);
    };

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

    // 0 or 1 snap is not enough to split a periodic intersection on a trimmed face.
    if (curve.periodic && uniq.size() < 2) {
        emitBySampling();
        return segs;
    }

    if (uniq.empty()) {
        if (curve.periodic) emitBySampling();
        return segs;
    }

    if (curve.periodic) {
        for (size_t i = 0; i + 1 < uniq.size(); ++i) pushInterval(uniq[i], uniq[i + 1], false);
        pushInterval(uniq.back(), uniq.front() + curve.period, false);
    } else {
        for (size_t i = 0; i + 1 < uniq.size(); ++i) pushInterval(uniq[i], uniq[i + 1], false);
    }
    return segs;
}

}  // namespace

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
        return intersectNurbsFaceWithPlane(iface, pln, frame, opt, boundary);
    }

    const std::vector<QuadricCurve> hits =
        intersectAnalyticSurfacePlane(surf, pln, opt.angular_tolerance, opt.geom_tolerance);

    std::vector<RawSegment> segs;

    for (const QuadricCurve& q : hits) {
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
        segs.insert(segs.end(), part.begin(), part.end());
    }
    return segs;
}

}  // namespace brepslicer
