#pragma once

#include <brepslicer/Types.h>
#include <geom/GeomUtil.h>

namespace brepslicer {

double contourLength(const Contour& c);
double layerArea(const Layer& layer, const SliceFrame& frame);
double layerPerimeter(const Layer& layer);

}  // namespace brepslicer
