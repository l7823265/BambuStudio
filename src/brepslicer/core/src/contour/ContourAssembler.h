#pragma once

#include <brepslicer/Types.h>
#include <geom/GeomUtil.h>
#include <intersect/AnalyticIntersect.h>

namespace brepslicer {

struct AssembleStats {
    int bridged = 0;
    double max_bridge = 0;
    int open_leftover = 0;
    int spurious_closed_removed = 0;
};

// After edge stitching: drop closed loops whose vertices coincide with an open contour's
// endpoints (self-loop pocket next to a real open cap arc).
void fixContours(std::vector<Contour>& contours, double gap_tol = 0.001);

Layer assembleLayer(double z,
                    std::vector<RawSegment> segs,
                    const SliceFrame& frame,
                    const SliceOptions& opt,
                    AssembleStats& stats);

double contourSignedArea(const Contour& c, const SliceFrame& frame);
Vec3 contourStart(const Contour& c);
Vec3 contourEnd(const Contour& c);

}  // namespace brepslicer
