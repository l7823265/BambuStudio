#pragma once

#include <brepslicer/Shape.h>

#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

namespace brepslicer {

class OccShape : public IShape {
public:
    explicit OccShape(const TopoDS_Shape& s);
    ShapeType type() const override;
    bool empty() const override;
    BBox bbox() const override;
    const TopoDS_Shape& occ() const { return shape_; }

protected:
    TopoDS_Shape shape_;
};

class OccEdge : public IEdge {
public:
    explicit OccEdge(const TopoDS_Edge& e);
    bool empty() const override;
    BBox bbox() const override;
    CurveData curve() const override;
    Vec3 eval(double t) const override;
    const TopoDS_Edge& occEdge() const { return edge_; }

private:
    TopoDS_Edge edge_;
};

class OccFace : public IFace {
public:
    explicit OccFace(const TopoDS_Face& f);
    bool empty() const override;
    BBox bbox() const override;
    SurfaceData surface() const override;
    PointClass classify(const Vec3& p, double tol) const override;
    UVBox uvDomain() const override;
    Vec3 evalUV(double u, double v) const override;
    bool derivUV(double u, double v, Vec3& Su, Vec3& Sv) const override;
    bool invertUV(const Vec3& p, double& u, double& v, double tol) const override;
    PointClass classifyUV(double u, double v, double tol) const override;
    std::vector<std::shared_ptr<IEdge>> edges() const override;
    std::vector<std::vector<std::shared_ptr<IEdge>>> wires() const override;
    void uvIsoSamples(std::vector<double>& u_samples,
                      std::vector<double>& v_samples) const override;
    const TopoDS_Face& occFace() const { return face_; }

private:
    TopoDS_Face face_;
};

const TopoDS_Shape& occShape(const IShape& s);

}  // namespace brepslicer
