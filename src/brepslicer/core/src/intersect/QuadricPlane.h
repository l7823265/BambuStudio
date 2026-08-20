#pragma once

#include <brepslicer/Geom.h>
#include <brepslicer/Shape.h>

#include <functional>
#include <vector>

namespace brepslicer {

enum class QuadricKind {
    Empty,
    Coplanar,
    Point,
    Line,
    Circle,
    Ellipse
};

struct QuadricCurve {
    QuadricKind kind = QuadricKind::Empty;
    bool degenerate = false;  // tangent point / tangent generator: mark, do not form a contour
    Vec3 point;
    Line3 line;
    Circle3 circle;
    Ellipse3 ellipse;
};

std::vector<QuadricCurve> intersectAnalyticSurfacePlane(const SurfaceData& surf,
                                                        const Plane& pln,
                                                        double angTol,
                                                        double geomTol);

void intersectCurvePlane(const CurveData& curve,
                         const Plane& pln,
                         double angTol,
                         double geomTol,
                         std::vector<Vec3>& pts,
                         const std::function<Vec3(double)>& eval);

}  // namespace brepslicer
