#include <brepslicer/Slicer.h>
#include <contour/ContourAssembler.h>
#include <geom/GeomUtil.h>
#include <verify/SliceVerify.h>

#include <algorithm>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using namespace brepslicer;

static int g_fail = 0;

static void expect(bool cond, const std::string& msg) {
    if (!cond) {
        std::cerr << "FAIL: " << msg << "\n";
        ++g_fail;
    } else {
        std::cout << "  ok  " << msg << "\n";
    }
}

static const Vec3 kAxis{0.3, 0.7, 1.0};
static const Vec3 kTrans{3.1, -2.4, 5.7};
static constexpr double kAng = 0.41;

static Vec3 rotateL3(const Vec3& v) {
    const Vec3 axis = normalized(kAxis);
    const double c = std::cos(kAng);
    const double s = std::sin(kAng);
    return v * c + cross(axis, v) * s + axis * (dot(axis, v) * (1.0 - c));
}

static std::shared_ptr<IShape> xf(const std::shared_ptr<IShape>& s) {
    return ShapeEngine::Builder().transform(s, {0, 0, 0}, kAxis, kAng, kTrans);
}

static SliceOptions worldSlice(const Vec3& n_local, double d_local) {
    const Vec3 n0 = normalized(n_local);
    const Vec3 n = rotateL3(n0);
    const Vec3 p = rotateL3(n0 * d_local) + kTrans;
    SliceOptions opt;
    opt.normal = n;
    opt.explicit_heights = {planeOffset(n, p)};
    opt.tolerance = 1e-4;
    opt.geom_tolerance = 1e-9;
    opt.angular_tolerance = 1e-10;
    return opt;
}

static void expectOcctRel(const SliceResult& r, const std::shared_ptr<IShape>& shape,
                          const std::string& tag) {
    const VerifyReport v = verifyAgainstRef(r, shape);
    for (const auto& e : v.errors) std::cerr << "    ERROR " << tag << ": " << e << "\n";
    for (const auto& n : v.notes) std::cout << "    " << n << "\n";
    expect(v.ok, tag + " verifySlice+OCCT rel<=1e-4");
}

static int countType(const Layer& L, SegmentType t) {
    int n = 0;
    for (const Contour& c : L.contours)
        for (const Segment& s : c.segments)
            if (s.type == t) ++n;
    return n;
}

static bool hasEvent(const SliceResult& r, const std::string& ev) {
    for (const auto& e : r.topology_events) {
        if (e.event == ev) return true;
    }
    return false;
}

static void testCube() {
    std::cout << "[cube, transformed]\n";
    const auto box = xf(ShapeEngine::Builder().makeBox(10.0, 10.0, 10.0));

    {
        const SliceResult r = sliceShape(box, worldSlice({0, 0, 1}, 5.0));
        expect(!r.layers.empty() && r.layers[0].contours.size() == 1, "generic: 1 contour");
        if (!r.layers.empty() && !r.layers[0].contours.empty()) {
            expect(r.layers[0].contours[0].closed, "generic closed");
            expect(countType(r.layers[0], SegmentType::Line) == 4, "4 lines");
            const SliceFrame f = makeSliceFrame(r.normal);
            expect(std::abs(std::abs(contourSignedArea(r.layers[0].contours[0], f)) - 100.0) < 1e-6,
                   "area 100");
        }
        expectOcctRel(r, box, "cube generic");
    }

    {
        const SliceResult r = sliceShape(box, worldSlice({-2, -2, 5}, 0.0));
        expect(!r.layers.empty(), "vertex layer");
        if (!r.layers.empty()) {
            expect(r.layers[0].contours.size() == 1, "through vertex: 1 polygon");
            if (!r.layers[0].contours.empty()) {
                expect(r.layers[0].contours[0].closed, "vertex contour closed");
                expect(countType(r.layers[0], SegmentType::Line) >= 3, ">=3 lines at vertex");
            }
        }
        expectOcctRel(r, box, "cube through vertex");
    }

    {
        const SliceResult r = sliceShape(box, worldSlice({0, -1, 1}, 0.0));
        expect(!r.layers.empty() && r.layers[0].contours.size() == 1, "through edge: 1 contour");
        if (!r.layers.empty() && !r.layers[0].contours.empty()) {
            expect(r.layers[0].contours[0].closed, "edge contour closed");
            expect(countType(r.layers[0], SegmentType::Line) >= 3, "polygon includes the shared edge");
        }
        expectOcctRel(r, box, "cube through edge");
    }

    {
        const SliceResult r = sliceShape(box, worldSlice({0, 0, 1}, 0.0));
        expect(!r.layers.empty() && r.layers[0].contours.size() == 1, "coplanar: 1 contour");
        if (!r.layers.empty() && !r.layers[0].contours.empty()) {
            expect(r.layers[0].contours[0].coplanar, "coplanar:true");
            expect(r.layers[0].contours[0].closed, "coplanar closed");
        }
        expectOcctRel(r, box, "cube coplanar");
    }
}

static void testCylinder() {
    std::cout << "[cylinder r=5, transformed]\n";
    const auto cyl = xf(ShapeEngine::Builder().makeCylinder(5.0, 10.0));

    {
        const SliceResult r = sliceShape(cyl, worldSlice({0, 0, 1}, 5.0));
        expect(!r.layers.empty() && r.layers[0].contours.size() == 1, "perp: 1 circle");
        bool arc = false;
        if (!r.layers.empty() && !r.layers[0].contours.empty()) {
            for (const Segment& s : r.layers[0].contours[0].segments) {
                if (s.type != SegmentType::Arc) continue;
                arc = true;
                expect(std::abs(s.radius - 5.0) <= 1e-9, "radius <= 1e-9");
            }
        }
        expect(arc, "analytic arc");
        expectOcctRel(r, cyl, "cyl perp");
    }

    {
        const SliceResult r = sliceShape(cyl, worldSlice({1, 0, 0}, 5.0));
        int nClosed = 0;
        if (!r.layers.empty()) {
            for (const Contour& c : r.layers[0].contours)
                if (c.closed && !c.segments.empty()) ++nClosed;
        }
        expect(nClosed == 0, "tangent: no fake closed contour");
        expect(hasEvent(r, "tangent_line") || hasEvent(r, "tangent_point"),
               "tangent marked in topology_events");
        expectOcctRel(r, cyl, "cyl tangent");
    }
}

static void testCone() {
    std::cout << "[cone R=5 H=10, transformed]\n";
    const auto cone = xf(ShapeEngine::Builder().makeCone(5.0, 0.0, 10.0));

    {
        const SliceResult r = sliceShape(cone, worldSlice({0, 0, 1}, 5.0));
        expect(!r.layers.empty() && r.layers[0].contours.size() == 1, "cone mid: 1 circle");
        bool arc = false;
        if (!r.layers.empty() && !r.layers[0].contours.empty()) {
            for (const Segment& s : r.layers[0].contours[0].segments) {
                if (s.type != SegmentType::Arc) continue;
                arc = true;
                expect(std::abs(s.radius - 2.5) <= 1e-9, "cone radius 2.5");
            }
        }
        expect(arc, "cone analytic arc");
        expectOcctRel(r, cone, "cone mid");
    }

    {
        const SliceResult r = sliceShape(cone, worldSlice({0, 0, 1}, 10.0));
        int nClosed = 0;
        if (!r.layers.empty()) {
            for (const Contour& c : r.layers[0].contours)
                if (c.closed && !c.segments.empty()) ++nClosed;
        }
        expect(nClosed == 0, "apex: no fake contour");
        expect(hasEvent(r, "tangent_point") || hasEvent(r, "tangent_line") ||
                   hasEvent(r, "degenerate_ring"),
               "apex marked degenerate");
        expectOcctRel(r, cone, "cone apex");
    }
}

static void testTorus() {
    std::cout << "[torus R=5 r=2, transformed]\n";
    const auto tor = xf(ShapeEngine::Builder().makeTorus(5.0, 2.0));

    {
        const SliceResult r = sliceShape(tor, worldSlice({0, 0, 1}, 0.0));
        expect(!r.layers.empty() && r.layers[0].contours.size() == 2, "equator: 2 circles");
        if (!r.layers.empty()) {
            std::vector<double> rs;
            for (const Contour& c : r.layers[0].contours) {
                for (const Segment& s : c.segments) {
                    if (s.type == SegmentType::Arc) rs.push_back(s.radius);
                }
            }
            expect(rs.size() >= 2, "two analytic arcs");
            if (rs.size() >= 2) {
                std::sort(rs.begin(), rs.end());
                expect(std::abs(rs.front() - 3.0) <= 1e-9, "inner r=3");
                expect(std::abs(rs.back() - 7.0) <= 1e-9, "outer r=7");
            }
            int nInner = 0;
            for (const Contour& c : r.layers[0].contours)
                if (c.orientation == "inner") ++nInner;
            expect(nInner == 1, "outer+inner nesting");
        }
        expectOcctRel(r, tor, "torus equator");
    }

    {
        const SliceResult r = sliceShape(tor, worldSlice({1, 0, 0}, 7.0));
        int nClosed = 0;
        if (!r.layers.empty()) {
            for (const Contour& c : r.layers[0].contours)
                if (c.closed && !c.segments.empty()) ++nClosed;
        }
        expect(nClosed == 0, "torus outer tangent: no fake contour");
        expect(hasEvent(r, "tangent_point") || hasEvent(r, "tangent_line") ||
                   hasEvent(r, "degenerate_ring"),
               "torus tangent marked");
        expectOcctRel(r, tor, "torus tangent");
    }
}

static void testThinWall() {
    std::cout << "[thin wall gap=5e-4, transformed]\n";
    IShapeBuilder& b = ShapeEngine::Builder();
    const auto a = b.makeBox(10.0, 10.0, 10.0);
    const auto shifted = b.transform(a, {0, 0, 0}, {0, 0, 1}, 0.0, {10.0 + 5e-4, 0, 0});
    const auto shape = xf(b.compound({a, shifted}));
    const SliceResult r = sliceShape(shape, worldSlice({0, 0, 1}, 5.0));
    expect(!r.layers.empty() && r.layers[0].contours.size() == 2, "two contours, not merged");
    expectOcctRel(r, shape, "thin wall");
}

static void testNurbsBump() {
    std::cout << "[nurbs bump, interior closed loop, transformed]\n";
    const auto face = xf(ShapeEngine::Builder().makeBumpFace(10.0, 10.0, 4.0));
    const SliceResult r = sliceShape(face, worldSlice({0, 0, 1}, 2.0));
    expect(!r.layers.empty(), "has layer");
    int nClosed = 0;
    int nBsp = 0;
    double maxFit = 0;
    if (!r.layers.empty()) {
        for (const Contour& c : r.layers[0].contours) {
            if (c.closed && !c.segments.empty()) ++nClosed;
            for (const Segment& s : c.segments) {
                if (s.type == SegmentType::BSpline) {
                    ++nBsp;
                    maxFit = std::max(maxFit, s.fit_error);
                }
            }
        }
    }
    expect(nClosed == 1, "interior loop: 1 closed contour (no edge hits)");
    expect(nBsp >= 1, "cubic bspline output");
    expect(maxFit <= 1e-4, "fit_error Hausdorff <= 1e-4, got " + std::to_string(maxFit));
    expectOcctRel(r, face, "nurbs bump");
}

int main() {
    try {
        testCube();
        testCylinder();
        testCone();
        testTorus();
        testThinWall();
        testNurbsBump();
    } catch (const std::exception& ex) {
        std::cerr << "exception: " << ex.what() << "\n";
        return 1;
    }
    if (g_fail) {
        std::cerr << g_fail << " check(s) failed\n";
        return 2;
    }
    std::cout << "ALL SELFTESTS PASSED\n";
    return 0;
}
