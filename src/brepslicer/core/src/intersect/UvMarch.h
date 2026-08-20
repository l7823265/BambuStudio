#pragma once

#include <intersect/AnalyticIntersect.h>

namespace brepslicer {

std::vector<RawSegment> intersectNurbsFaceWithPlane(const FaceRecord& iface,
                                                    const Plane& pln,
                                                    const SliceFrame& frame,
                                                    const SliceOptions& opt,
                                                    const std::vector<Vec3>& boundary3d);

}  // namespace brepslicer
