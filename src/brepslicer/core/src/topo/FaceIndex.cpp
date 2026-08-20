#include <topo/FaceIndex.h>

namespace brepslicer {

void FaceIndex::build(const std::vector<FaceRecord>& faces, const BBox& shapeBox) {
    faces_ = faces;
    shape_box_ = shapeBox;
}

std::vector<const FaceRecord*> FaceIndex::query(const Vec3& n, double d, double tol) const {
    std::vector<const FaceRecord*> out;
    out.reserve(faces_.size());
    for (const FaceRecord& f : faces_) {
        if (!f.box.valid) {
            out.push_back(&f);
            continue;
        }
        double dmin = 0, dmax = 0;
        aabbRangeAlong(f.box, n, dmin, dmax);
        if (d + tol >= dmin && d - tol <= dmax) out.push_back(&f);
    }
    return out;
}

void FaceIndex::shapeRange(const Vec3& n, double& dmin, double& dmax) const {
    if (!shape_box_.valid) {
        dmin = 0;
        dmax = 0;
        return;
    }
    aabbRangeAlong(shape_box_, n, dmin, dmax);
}

}  // namespace brepslicer
