#pragma once

#include <brepslicer/Shape.h>

#include <vector>

namespace brepslicer {

class FaceIndex {
public:
    void build(const std::vector<FaceRecord>& faces, const BBox& shapeBox);

    const std::vector<FaceRecord>& faces() const { return faces_; }

    // Conservative overlap of face AABB with plane n·x = d.
    std::vector<const FaceRecord*> query(const Vec3& n, double d, double tol) const;

    void shapeRange(const Vec3& n, double& dmin, double& dmax) const;

private:
    std::vector<FaceRecord> faces_;
    BBox shape_box_;
};

}  // namespace brepslicer
