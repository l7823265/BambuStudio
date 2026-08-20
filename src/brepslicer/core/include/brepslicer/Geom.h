#pragma once

#include <brepslicer/Types.h>

#include <algorithm>
#include <cmath>

namespace brepslicer {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTwoPi = 2.0 * kPi;

struct Vec2 {
    double x = 0;
    double y = 0;
};

inline Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Vec3 operator-(const Vec3& a) { return {-a.x, -a.y, -a.z}; }
inline Vec3 operator*(const Vec3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline Vec3 operator*(double s, const Vec3& a) { return a * s; }
inline Vec3 operator/(const Vec3& a, double s) { return {a.x / s, a.y / s, a.z / s}; }

inline double dot(const Vec3& a, const Vec3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline Vec3 cross(const Vec3& a, const Vec3& b) {
    return {a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
}
inline double length2(const Vec3& a) { return dot(a, a); }
inline double length(const Vec3& a) { return std::sqrt(length2(a)); }

inline Vec3 normalized(const Vec3& a) {
    const double m = length(a);
    return (m > 0) ? a / m : Vec3{};
}

inline bool unitize(Vec3 v, Vec3& out, double tol) {
    const double m = length(v);
    if (m <= tol) return false;
    out = v / m;
    return true;
}

inline double wrapTwoPi(double a) {
    a = std::fmod(a, kTwoPi);
    if (a < 0) a += kTwoPi;
    return a;
}

inline bool almostParallel(const Vec3& a, const Vec3& b, double angTol) {
    const double c = std::abs(dot(normalized(a), normalized(b)));
    return c >= std::cos(angTol) || c >= 1.0 - 1e-15;
}

// Plane n·x = d, |n| = 1.
struct Plane {
    Vec3 n{0, 0, 1};
    double d = 0;
};

inline Plane makePlane(const Vec3& n, double offset) {
    const Vec3 u = normalized(n);
    return {u, offset};
}

inline double signedPlaneDist(const Plane& pln, const Vec3& p) { return dot(pln.n, p) - pln.d; }
inline double planeOffset(const Vec3& n, const Vec3& p) { return dot(normalized(n), p); }

struct Line3 {
    Vec3 p;
    Vec3 d;  // unit
    Vec3 at(double t) const { return p + d * t; }
};

struct Circle3 {
    Vec3 center;
    Vec3 n;
    Vec3 xdir;
    double radius = 0;
    Vec3 ydir() const { return cross(n, xdir); }
    Vec3 at(double t) const { return center + xdir * (radius * std::cos(t)) + ydir() * (radius * std::sin(t)); }
    double param(const Vec3& p) const {
        const Vec3 v = p - center;
        return std::atan2(dot(v, ydir()), dot(v, xdir));
    }
};

struct Ellipse3 {
    Vec3 center;
    Vec3 n;
    Vec3 major;  // unit
    double a = 0;
    double b = 0;
    Vec3 minor() const { return cross(n, major); }
    Vec3 at(double t) const { return center + major * (a * std::cos(t)) + minor() * (b * std::sin(t)); }
    double param(const Vec3& p) const {
        const Vec3 v = p - center;
        return std::atan2(dot(v, minor()) / std::max(b, 1e-30), dot(v, major) / std::max(a, 1e-30));
    }
};

struct Cylinder3 {
    Vec3 origin;
    Vec3 axis;
    double radius = 0;
};

struct Sphere3 {
    Vec3 center;
    double radius = 0;
};

struct Cone3 {
    Vec3 apex;
    Vec3 axis;
    double semi_angle = 0;
};

struct Torus3 {
    Vec3 center;
    Vec3 axis;
    double major_radius = 0;
    double minor_radius = 0;
};

enum class SurfaceKind { Plane, Cylinder, Sphere, Cone, Torus, Other };
enum class CurveKind { Line, Circle, Ellipse, Other };
enum class PointClass { In, On, Out };

struct SurfaceData {
    SurfaceKind kind = SurfaceKind::Other;
    Plane plane;
    Cylinder3 cylinder;
    Sphere3 sphere;
    Cone3 cone;
    Torus3 torus;
};

struct CurveData {
    CurveKind kind = CurveKind::Other;
    Line3 line;
    Circle3 circle;
    Ellipse3 ellipse;
    double first = 0;
    double last = 0;
    bool reversed = false;
    bool degenerated = false;
};

struct BBox {
    bool valid = false;
    double xmin = 0, ymin = 0, zmin = 0;
    double xmax = 0, ymax = 0, zmax = 0;
};

struct UVBox {
    double umin = 0, umax = 1, vmin = 0, vmax = 1;
    bool periodic_u = false;
    bool periodic_v = false;
    double period_u = 0;
    double period_v = 0;
};

inline void wrapUV(const UVBox& d, double& u, double& v) {
    if (d.periodic_u && d.period_u > 0) {
        while (u < d.umin) u += d.period_u;
        while (u >= d.umin + d.period_u) u -= d.period_u;
    }
    if (d.periodic_v && d.period_v > 0) {
        while (v < d.vmin) v += d.period_v;
        while (v >= d.vmin + d.period_v) v -= d.period_v;
    }
}

inline double uvDelta(double a, double b, bool periodic, double period) {
    double d = b - a;
    if (periodic && period > 0) {
        while (d > 0.5 * period) d -= period;
        while (d < -0.5 * period) d += period;
    }
    return d;
}

struct SliceFrame {
    Vec3 n{0, 0, 1};
    Vec3 x{1, 0, 0};
    Vec3 y{0, 1, 0};

    Vec2 toXY(const Vec3& p) const { return {dot(p, x), dot(p, y)}; }

    double angle(const Vec3& center, const Vec3& p) const {
        const Vec3 v = p - center;
        return std::atan2(dot(v, y), dot(v, x));
    }
};

inline SliceFrame makeSliceFrame(const Vec3& n_in) {
    SliceFrame f;
    f.n = normalized(n_in);
    const Vec3 tmp = (std::abs(f.n.z) < 0.9) ? Vec3{0, 0, 1} : Vec3{0, 1, 0};
    f.x = normalized(cross(tmp, f.n));
    f.y = cross(f.n, f.x);
    return f;
}

inline Circle3 makeCircle(const Vec3& center, const Vec3& n, double r) {
    const SliceFrame f = makeSliceFrame(n);
    return {center, f.n, f.x, r};
}

inline Ellipse3 makeEllipse(const Vec3& center, const Vec3& n, const Vec3& major_in, double a, double b) {
    Ellipse3 e;
    e.center = center;
    e.n = normalized(n);
    if (a >= b) {
        e.major = normalized(major_in);
        e.a = a;
        e.b = b;
    } else {
        e.major = normalized(cross(e.n, normalized(major_in)));
        e.a = b;
        e.b = a;
    }
    return e;
}

inline void aabbRangeAlong(const BBox& b, const Vec3& n, double& dmin, double& dmax) {
    const double xs[2] = {b.xmin, b.xmax};
    const double ys[2] = {b.ymin, b.ymax};
    const double zs[2] = {b.zmin, b.zmax};
    dmin = 1e300;
    dmax = -1e300;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {
                const double d = n.x * xs[i] + n.y * ys[j] + n.z * zs[k];
                dmin = std::min(dmin, d);
                dmax = std::max(dmax, d);
            }
}

}  // namespace brepslicer
