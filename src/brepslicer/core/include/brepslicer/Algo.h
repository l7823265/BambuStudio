#pragma once

#include <brepslicer/Shape.h>
#include <brepslicer/Geom.h>

#include <memory>
#include <string>
#include <vector>

namespace brepslicer {

struct RawSegment;

class IBrepKernel {
public:
    virtual ~IBrepKernel() = default;
    virtual std::shared_ptr<IShape> readStep(const std::string& path) = 0;
    virtual std::vector<FaceRecord> exploreFaces(const IShape& shape) = 0;
};

class IShapeBuilder {
public:
    virtual ~IShapeBuilder() = default;
    virtual std::shared_ptr<IShape> makeBox(double dx, double dy, double dz) = 0;
    virtual std::shared_ptr<IShape> makeCylinder(double radius, double height) = 0;
    virtual std::shared_ptr<IShape> makeCone(double r1, double r2, double height) = 0;
    virtual std::shared_ptr<IShape> makeTorus(double major_radius, double minor_radius) = 0;
    // Interpolating Gaussian bump patch (NURBS). Interior z-slice yields a closed loop
    // that does not meet any edge — the L2 seed-search case.
    virtual std::shared_ptr<IShape> makeBumpFace(double dx, double dy, double bump_height) = 0;
    virtual std::shared_ptr<IShape> transform(const std::shared_ptr<IShape>& shape,
                                              const Vec3& axis_origin,
                                              const Vec3& axis_dir,
                                              double angle_rad,
                                              const Vec3& translation) = 0;
    // Affine map p' = R * p + t. R is row-major 3x3 (r00,r01,r02, r10,...), t is (tx,ty,tz).
    virtual std::shared_ptr<IShape> transformAffine(const std::shared_ptr<IShape>& shape,
                                                    const double R[9],
                                                    const double t[3]) = 0;
    virtual std::shared_ptr<IShape> compound(const std::vector<std::shared_ptr<IShape>>& parts) = 0;
};

class ISectionRef {
public:
    virtual ~ISectionRef() = default;
    virtual double sectionLength(const IShape& shape, const Plane& plane) = 0;
    // Drop face-trim segments whose interior samples are not on the OCCT section wire.
    virtual void clipSegmentsToSection(const IShape& shape, const Plane& plane,
                                       std::vector<RawSegment>& segs, double tol) = 0;
};

class IEngineStrategy {
public:
    virtual ~IEngineStrategy() = default;
    virtual IBrepKernel* GetKernel() = 0;
    virtual IShapeBuilder* GetShapeBuilder() = 0;
    virtual ISectionRef* GetSectionRef() = 0;
};

}  // namespace brepslicer
