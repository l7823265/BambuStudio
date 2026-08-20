#pragma once

#include <brepslicer/Types.h>
#include <brepslicer/Geom.h>

#include <vector>

namespace brepslicer {

Vec3 bsplineEval(const Segment& s, double t01);
void bsplineSample(const Segment& s, int n, std::vector<Vec3>& out);
double bsplineLength(const Segment& s);
double bsplineHausdorff(const Segment& s, const std::vector<Vec3>& poly);

// Fit an open or closed cubic B-spline. fit_error is polyline↔spline Hausdorff.
Segment fitCubicBSpline(const std::vector<Vec3>& pts, bool closed, double max_err);

// C0 cubic encoding of the polyline (linear cubics). Interpolates every vertex.
Segment cubicFromPolyline(const std::vector<Vec3>& pts, bool closed);

void reverseBSpline(Segment& s);

}  // namespace brepslicer
