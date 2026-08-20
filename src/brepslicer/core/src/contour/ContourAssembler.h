#pragma once

#include <brepslicer/Types.h>
#include <geom/GeomUtil.h>
#include <intersect/AnalyticIntersect.h>

namespace brepslicer {

struct AssembleStats {
    int bridged = 0;
    double max_bridge = 0;
    int open_leftover = 0;
};

Layer assembleLayer(double z,
                    std::vector<RawSegment> segs,
                    const SliceFrame& frame,
                    const SliceOptions& opt,
                    AssembleStats& stats);

double contourSignedArea(const Contour& c, const SliceFrame& frame);
Vec3 contourStart(const Contour& c);
Vec3 contourEnd(const Contour& c);

}  // namespace brepslicer
