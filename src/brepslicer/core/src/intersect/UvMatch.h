#pragma once

#include <intersect/AnalyticIntersect.h>

namespace brepslicer {

// Wei et al., Sci. Rep. 2025: plane ∩ NURBS by iso-curve bisection + B-spline fit.
// Pipeline: seeding → trim filter → fitting → makeSeg.
std::vector<RawSegment> intersectNurbsFaceWithPlaneUvMatch(const FaceRecord& iface,
                                                          const Plane& pln,
                                                          const SliceFrame& frame,
                                                          const SliceOptions& opt,
                                                          const std::vector<Vec3>& boundary3d);

}  // namespace brepslicer
