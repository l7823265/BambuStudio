#pragma once

#include <brepslicer/Geom.h>

namespace brepslicer {

inline double dist(const Vec3& a, const Vec3& b) { return length(a - b); }

}  // namespace brepslicer
