#pragma once

#include <brepslicer/Types.h>
#include <brepslicer/Shape.h>
#include <brepslicer/Geom.h>

#include <vector>

namespace brepslicer {

struct RawSegment {
    Segment geom;
    int solid_id = -1;
    int shell_id = -1;
    int face_id = -1;
    bool coplanar = false;
    bool closed_loop = false;
    bool degenerate = false;
    std::string degen_event;
    Vec3 degen_point;
};

std::vector<RawSegment> intersectFaceWithPlane(const FaceRecord& iface,
                                               const Plane& pln,
                                               const SliceFrame& frame,
                                               const SliceOptions& opt);

}  // namespace brepslicer
