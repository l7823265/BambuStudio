#include <intersect/QuadricPlane.h>

#include <algorithm>
#include <cmath>

namespace brepslicer {
namespace {

int solveQuad(double a, double b, double c, double (&roots)[2], double tol) {
    if (std::abs(a) <= tol) {
        if (std::abs(b) <= tol) return 0;
        roots[0] = -c / b;
        return 1;
    }
    const double disc = b * b - 4.0 * a * c;
    if (disc < -tol * tol) return 0;
    if (disc <= tol * tol) {
        roots[0] = -0.5 * b / a;
        return 1;
    }
    const double s = std::sqrt(disc);
    const double q = -0.5 * (b + (b >= 0 ? s : -s));
    roots[0] = q / a;
    roots[1] = c / q;
    if (roots[0] > roots[1]) std::swap(roots[0], roots[1]);
    return 2;
}

QuadricCurve qcPoint(const Vec3& p, bool degenerate) {
    QuadricCurve c;
    c.kind = QuadricKind::Point;
    c.degenerate = degenerate;
    c.point = p;
    return c;
}

QuadricCurve qcLine(const Vec3& p, const Vec3& v, bool degenerate) {
    QuadricCurve c;
    c.kind = QuadricKind::Line;
    c.degenerate = degenerate;
    c.line = {p, normalized(v)};
    c.point = p;
    return c;
}

QuadricCurve qcCircle(const Vec3& center, const Vec3& n, double r, bool degenerate) {
    QuadricCurve c;
    c.kind = QuadricKind::Circle;
    c.degenerate = degenerate;
    c.circle = makeCircle(center, n, r);
    c.point = center;
    return c;
}

QuadricCurve qcEllipse(const Vec3& center, const Vec3& n, const Vec3& major, double a, double b) {
    QuadricCurve c;
    c.kind = QuadricKind::Ellipse;
    c.ellipse = makeEllipse(center, n, major, a, b);
    c.point = center;
    return c;
}

std::vector<QuadricCurve> planePlane(const Plane& a, const Plane& b, double angTol, double geomTol) {
    Vec3 dir;
    if (!unitize(cross(a.n, b.n), dir, std::sin(angTol))) {
        const Vec3 loc = a.n * a.d;
        if (std::abs(signedPlaneDist(b, loc)) <= geomTol) {
            QuadricCurve c;
            c.kind = QuadricKind::Coplanar;
            return {c};
        }
        return {};
    }
    const double d1 = a.d;
    const double d2 = b.d;
    const double cdot = dot(a.n, b.n);
    const double det = 1.0 - cdot * cdot;
    const double s1 = (d1 - d2 * cdot) / det;
    const double s2 = (d2 - d1 * cdot) / det;
    const Vec3 p = a.n * s1 + b.n * s2;
    return {qcLine(p, dir, false)};
}

std::vector<QuadricCurve> planeSphere(const Plane& pln, const Sphere3& sph, double geomTol) {
    const Vec3 c = sph.center;
    const double r = sph.radius;
    const Vec3 n = pln.n;
    const double h = signedPlaneDist(pln, c);
    const Vec3 foot = c - n * h;
    const double disc = r * r - h * h;
    if (disc < -2.0 * r * geomTol - geomTol * geomTol) return {};
    if (disc <= 2.0 * r * geomTol + geomTol * geomTol) {
        return {qcPoint(foot, true)};
    }
    return {qcCircle(foot, n, std::sqrt(disc), false)};
}

std::vector<QuadricCurve> planeCylinder(const Plane& pln, const Cylinder3& cyl, double angTol,
                                        double geomTol) {
    const Vec3 n = pln.n;
    const Vec3 v = cyl.axis;
    const Vec3 a = cyl.origin;
    const double r = cyl.radius;
    const double nv = dot(n, v);
    const double delta = signedPlaneDist(pln, a);

    if (std::abs(nv) <= std::sin(angTol)) {
        Vec3 u;
        if (!unitize(cross(v, n), u, 1e-14)) return {};
        const double disc = r * r - delta * delta;
        const Vec3 base = a - n * delta;
        if (disc < -2.0 * r * geomTol) return {};
        if (disc <= 2.0 * r * geomTol) {
            return {qcLine(base, v, true)};
        }
        const double w = std::sqrt(disc);
        return {qcLine(base + u * w, v, false), qcLine(base - u * w, v, false)};
    }

    const double t = -delta / nv;
    const Vec3 P = a + v * t;
    const double absnv = std::abs(nv);
    if (absnv >= std::cos(angTol)) {
        return {qcCircle(P, n, r, false)};
    }
    const Vec3 v_on_plane = v - n * nv;
    Vec3 major;
    if (!unitize(v_on_plane, major, 1e-14)) return {qcCircle(P, n, r, false)};
    const double aa = r / absnv;
    const double bb = r;
    if (std::abs(aa - bb) <= geomTol) return {qcCircle(P, n, r, false)};
    return {qcEllipse(P, n, major, aa, bb)};
}

std::vector<QuadricCurve> planeConeApex(const Plane& pln, const Vec3& A, const Vec3& v, double alpha,
                                        double angTol, double geomTol) {
    const Vec3 n = pln.n;
    const double nv = std::abs(dot(n, v));
    const double sin_phi = nv;
    const double sin_a = std::sin(alpha);
    const SliceFrame ax = makeSliceFrame(n);
    const Vec3 e1 = ax.x;
    const Vec3 e2 = ax.y;
    const double ev1 = dot(e1, v);
    const double ev2 = dot(e2, v);
    const double target = std::cos(alpha);
    auto linesFor = [&](double tgt, std::vector<QuadricCurve>& out) {
        const double m2 = ev1 * ev1 + ev2 * ev2;
        if (m2 <= 1e-30) return;
        const double proj = tgt;
        if (proj * proj > m2 + 1e-14) return;
        const Vec3 mdir = e1 * ev1 + e2 * ev2;
        const Vec3 base = A;
        if (std::abs(proj * proj - m2) <= 1e-12) {
            Vec3 g;
            if (unitize(mdir, g, geomTol)) out.push_back(qcLine(base, g, true));
            return;
        }
        const Vec3 g0 = mdir * (tgt / m2);
        Vec3 perp;
        if (!unitize(cross(n, mdir), perp, geomTol)) return;
        const double h = std::sqrt(std::max(0.0, 1.0 - (tgt * tgt) / m2));
        const Vec3 gA = g0 + perp * h;
        const Vec3 gB = g0 - perp * h;
        Vec3 dA, dB;
        if (unitize(gA, dA, geomTol)) out.push_back(qcLine(base, dA, false));
        if (unitize(gB, dB, geomTol)) out.push_back(qcLine(base, dB, false));
    };
    std::vector<QuadricCurve> out;
    if (sin_phi + 1e-12 < sin_a) {
        linesFor(target, out);
        if (out.empty()) linesFor(-target, out);
        if (out.empty()) out.push_back(qcPoint(A, true));
        return out;
    }
    if (std::abs(sin_phi - sin_a) <= std::sin(angTol) + 1e-12) {
        linesFor(target, out);
        if (out.empty()) linesFor(-target, out);
        for (auto& c : out) c.degenerate = true;
        if (out.empty()) out.push_back(qcPoint(A, true));
        return out;
    }
    return {qcPoint(A, true)};
}

std::vector<QuadricCurve> planeCone(const Plane& pln, const Cone3& cone, double angTol,
                                    double geomTol) {
    const Vec3 A = cone.apex;
    const Vec3 v = cone.axis;
    const double alpha = cone.semi_angle;
    const Vec3 n = pln.n;
    const double h = signedPlaneDist(pln, A);
    if (std::abs(h) <= geomTol) {
        return planeConeApex(pln, A, v, alpha, angTol, geomTol);
    }

    const double nv = dot(n, v);
    if (std::abs(nv) >= std::cos(angTol)) {
        const double t = -h / nv;
        const Vec3 C = A + v * t;
        const double rad = length(A - C) * std::tan(alpha);
        if (rad <= geomTol) return {qcPoint(C, true)};
        return {qcCircle(C, n, rad, false)};
    }

    const SliceFrame ax = makeSliceFrame(n);
    const Vec3 e1 = ax.x;
    const Vec3 e2 = ax.y;
    const Vec3 Q = n * pln.d;
    const Vec3 P0 = Q - A;
    const double c2 = std::cos(alpha) * std::cos(alpha);
    const double pv = dot(P0, v);
    const double evu = dot(e1, v);
    const double evv = dot(e2, v);
    const double pex = dot(P0, e1);
    const double pey = dot(P0, e2);
    const double auu = evu * evu - c2;
    const double avv = evv * evv - c2;
    const double hh = evu * evv;
    const double gu = pv * evu - c2 * pex;
    const double fv = pv * evv - c2 * pey;
    const double cc = pv * pv - c2 * dot(P0, P0);

    const double det = auu * avv - hh * hh;
    if (det <= 1e-18) {
        return {qcPoint(Q, true)};
    }

    const double u0 = (-gu * avv + fv * hh) / det;
    const double v0 = (-fv * auu + gu * hh) / det;
    const Vec3 center = Q + e1 * u0 + e2 * v0;

    const double theta = 0.5 * std::atan2(2.0 * hh, auu - avv);
    const double ct = std::cos(theta);
    const double st = std::sin(theta);
    const double ap = auu * ct * ct + 2.0 * hh * ct * st + avv * st * st;
    const double bp = auu * st * st - 2.0 * hh * ct * st + avv * ct * ct;
    const double val = auu * u0 * u0 + 2.0 * hh * u0 * v0 + avv * v0 * v0 + 2.0 * gu * u0 +
                       2.0 * fv * v0 + cc;
    const double rhs = -val;
    if (ap * rhs <= 0 || bp * rhs <= 0) return {};
    const double ra = std::sqrt(rhs / ap);
    const double rb = std::sqrt(rhs / bp);
    if (!(ra > geomTol && rb > geomTol)) return {qcPoint(center, true)};
    const Vec3 major = normalized(e1 * ct + e2 * st);
    if (std::abs(ra - rb) <= geomTol) return {qcCircle(center, n, 0.5 * (ra + rb), false)};
    return {qcEllipse(center, n, major, ra, rb)};
}

std::vector<QuadricCurve> planeTorus(const Plane& pln, const Torus3& tor, double angTol,
                                     double geomTol) {
    const Vec3 C = tor.center;
    const Vec3 k = tor.axis;
    const Vec3 n = pln.n;
    const double R = tor.major_radius;
    const double r = tor.minor_radius;
    const double h = signedPlaneDist(pln, C);

    if (almostParallel(n, k, angTol)) {
        const Vec3 foot = C - n * h;
        const double disc = r * r - h * h;
        if (disc < -2.0 * r * geomTol) return {};
        if (disc <= 2.0 * r * geomTol) {
            return {qcCircle(foot, n, R, std::abs(R) <= geomTol)};
        }
        const double s = std::sqrt(disc);
        std::vector<QuadricCurve> out;
        const double rho1 = R + s;
        const double rho2 = R - s;
        if (rho1 > geomTol) out.push_back(qcCircle(foot, n, rho1, false));
        if (rho2 > geomTol) out.push_back(qcCircle(foot, n, rho2, false));
        else if (std::abs(rho2) <= geomTol) out.push_back(qcPoint(foot, true));
        return out;
    }

    if (std::abs(dot(n, k)) <= std::sin(angTol)) {
        const double dist = std::abs(h);
        if (dist > R + r + geomTol) return {};
        if (std::abs(dist - (R + r)) <= geomTol || std::abs(dist - std::abs(R - r)) <= geomTol) {
            Vec3 u;
            if (!unitize(cross(k, n), u, 1e-14)) return {qcPoint(C, true)};
            const double side = (h >= 0) ? 1.0 : -1.0;
            const Vec3 p = C + n * (side * dist);
            return {qcPoint(p, true)};
        }
        if (dist <= geomTol) {
            Vec3 u;
            if (!unitize(cross(n, k), u, 1e-14)) return {};
            return {qcCircle(C + u * R, n, r, false), qcCircle(C - u * R, n, r, false)};
        }
        return {};
    }
    return {};
}

void linePlaneHits(const Line3& lin, const Plane& pln, double first, double last, double angTol,
                   double geomTol, std::vector<Vec3>& pts) {
    const double nd = dot(pln.n, lin.d);
    if (std::abs(nd) <= std::sin(angTol)) return;
    const double t = -signedPlaneDist(pln, lin.p) / nd;
    if (t < first - geomTol || t > last + geomTol) return;
    pts.push_back(lin.at(t));
}

void lineCircleHits(const Circle3& c, const Line3& lin, double geomTol, std::vector<Vec3>& pts) {
    const Vec3 C = c.center;
    const Vec3 D = lin.d;
    const Vec3 P = lin.p;
    const double t0 = dot(C - P, D);
    const Vec3 F = P + D * t0;
    const double dist = length(F - C);
    const double r = c.radius;
    const double disc = r * r - dist * dist;
    if (disc < -2.0 * r * geomTol) return;
    if (disc <= 2.0 * r * geomTol) {
        pts.push_back(F);
        return;
    }
    const double w = std::sqrt(disc);
    pts.push_back(F + D * w);
    pts.push_back(F - D * w);
}

void circlePlaneHits(const Circle3& c, const Plane& pln, double first, double last, double angTol,
                     double geomTol, std::vector<Vec3>& pts) {
    const Vec3 nc = c.n;
    const Vec3 n = pln.n;
    if (almostParallel(n, nc, angTol)) {
        if (std::abs(signedPlaneDist(pln, c.center)) <= geomTol) return;
        return;
    }
    Vec3 dir;
    if (!unitize(cross(n, nc), dir, 1e-14)) return;
    const double d1 = pln.d;
    const double d2 = dot(nc, c.center);
    const double cdot = dot(n, nc);
    const double det = 1.0 - cdot * cdot;
    const double s1 = (d1 - d2 * cdot) / det;
    const double s2 = (d2 - d1 * cdot) / det;
    const Vec3 origin = n * s1 + nc * s2;
    std::vector<Vec3> raw;
    lineCircleHits(c, Line3{origin, dir}, geomTol, raw);
    for (const Vec3& p : raw) {
        double u = wrapTwoPi(c.param(p));
        auto inRange = [&](double t) {
            if (last - first >= kTwoPi - 1e-12) return true;
            double tt = t;
            while (tt < first) tt += kTwoPi;
            return tt <= last + geomTol / std::max(c.radius, 1.0);
        };
        if (inRange(u) || inRange(u + kTwoPi) || inRange(u - kTwoPi)) pts.push_back(p);
    }
}

void ellipsePlaneHits(const Ellipse3& e, const Plane& pln, double first, double last, double angTol,
                      double geomTol, std::vector<Vec3>& pts) {
    const Vec3 ne = e.n;
    const Vec3 n = pln.n;
    if (almostParallel(n, ne, angTol)) return;
    Vec3 dir;
    if (!unitize(cross(n, ne), dir, 1e-14)) return;
    const double d1 = pln.d;
    const double d2 = dot(ne, e.center);
    const double cdot = dot(n, ne);
    const double det = 1.0 - cdot * cdot;
    const double s1 = (d1 - d2 * cdot) / det;
    const double s2 = (d2 - d1 * cdot) / det;
    const Vec3 origin = n * s1 + ne * s2;
    const Line3 lin{origin, dir};
    const Vec3 x = e.major;
    const Vec3 y = e.minor();
    const double a = e.a;
    const double b = e.b;
    const Vec3 p0 = lin.p - e.center;
    const double px = dot(p0, x);
    const double py = dot(p0, y);
    const double dx = dot(dir, x);
    const double dy = dot(dir, y);
    const double A = (dx * dx) / (a * a) + (dy * dy) / (b * b);
    const double B = 2.0 * (px * dx / (a * a) + py * dy / (b * b));
    const double C = (px * px) / (a * a) + (py * py) / (b * b) - 1.0;
    double roots[2];
    const int nr = solveQuad(A, B, C, roots, 1e-18);
    for (int i = 0; i < nr; ++i) {
        const Vec3 p = lin.at(roots[i]);
        double t = e.param(p);
        if (t < 0) t += kTwoPi;
        if (last - first >= kTwoPi - 1e-12 || (t >= first - 1e-8 && t <= last + 1e-8)) pts.push_back(p);
    }
}

void sampleCurvePlaneHits(const std::function<Vec3(double)>& eval, const Plane& pln, double first,
                          double last, double geomTol, std::vector<Vec3>& pts) {
    const int N = 64;
    auto f = [&](double t) { return signedPlaneDist(pln, eval(t)); };
    double t0 = first;
    double s0 = f(t0);
    if (std::abs(s0) <= geomTol) pts.push_back(eval(t0));
    for (int i = 1; i <= N; ++i) {
        const double t1 = first + (last - first) * (static_cast<double>(i) / N);
        const double s1 = f(t1);
        if (std::abs(s1) <= geomTol) {
            pts.push_back(eval(t1));
        } else if (s0 * s1 < 0.0) {
            double a = t0, b = t1;
            for (int k = 0; k < 40; ++k) {
                const double m = 0.5 * (a + b);
                if (f(a) * f(m) <= 0.0) b = m;
                else a = m;
            }
            pts.push_back(eval(0.5 * (a + b)));
        }
        t0 = t1;
        s0 = s1;
    }
}

}  // namespace

std::vector<QuadricCurve> intersectAnalyticSurfacePlane(const SurfaceData& surf, const Plane& pln,
                                                        double angTol, double geomTol) {
    switch (surf.kind) {
    case SurfaceKind::Plane:
        return planePlane(pln, surf.plane, angTol, geomTol);
    case SurfaceKind::Cylinder:
        return planeCylinder(pln, surf.cylinder, angTol, geomTol);
    case SurfaceKind::Sphere:
        return planeSphere(pln, surf.sphere, geomTol);
    case SurfaceKind::Cone:
        return planeCone(pln, surf.cone, angTol, geomTol);
    case SurfaceKind::Torus:
        return planeTorus(pln, surf.torus, angTol, geomTol);
    default:
        return {};
    }
}

void intersectCurvePlane(const CurveData& curve, const Plane& pln, double angTol, double geomTol,
                         std::vector<Vec3>& pts, const std::function<Vec3(double)>& eval) {
    switch (curve.kind) {
    case CurveKind::Line:
        linePlaneHits(curve.line, pln, curve.first, curve.last, angTol, geomTol, pts);
        break;
    case CurveKind::Circle:
        circlePlaneHits(curve.circle, pln, curve.first, curve.last, angTol, geomTol, pts);
        break;
    case CurveKind::Ellipse:
        ellipsePlaneHits(curve.ellipse, pln, curve.first, curve.last, angTol, geomTol, pts);
        break;
    default:
        sampleCurvePlaneHits(eval, pln, curve.first, curve.last, geomTol, pts);
        break;
    }
}

}  // namespace brepslicer
