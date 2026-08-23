#pragma once

#include <brepslicer/Geom.h>

#include <memory>
#include <vector>

namespace brepslicer {

class IShape {
public:
    enum class ShapeType { Compound, Solid, Shell, Face, Wire, Edge, Vertex, Other };
    virtual ~IShape() = default;
    virtual ShapeType type() const = 0;
    virtual bool empty() const = 0;
    virtual BBox bbox() const = 0;
};

class IEdge : public IShape {
public:
    ShapeType type() const override { return ShapeType::Edge; }
    virtual CurveData curve() const = 0;
    virtual Vec3 eval(double t) const = 0;
};

class IFace : public IShape {
public:
    ShapeType type() const override { return ShapeType::Face; }
    virtual SurfaceData surface() const = 0;
    virtual PointClass classify(const Vec3& p, double tol) const = 0;
    virtual UVBox uvDomain() const = 0;
    virtual Vec3 evalUV(double u, double v) const = 0;
    virtual bool derivUV(double u, double v, Vec3& Su, Vec3& Sv) const = 0;
    virtual bool invertUV(const Vec3& p, double& u, double& v, double tol) const = 0;
    virtual PointClass classifyUV(double u, double v, double tol) const = 0;
    virtual std::vector<std::shared_ptr<IEdge>> edges() const = 0;
    // Wires in explorer order; each edge already carries CurveData.reversed.
    virtual std::vector<std::vector<std::shared_ptr<IEdge>>> wires() const = 0;
    // Preferred iso-parameter samples (e.g. unique NURBS knots). Empty → caller densifies.
    virtual void uvIsoSamples(std::vector<double>& /*u_samples*/,
                              std::vector<double>& /*v_samples*/) const {}
};

class ISolid : public IShape {
public:
    ShapeType type() const override { return ShapeType::Solid; }
};

struct FaceRecord {
    std::shared_ptr<IFace> face;
    int solid_id = -1;
    int shell_id = -1;
    int face_id = -1;
    BBox box;
};

}  // namespace brepslicer
