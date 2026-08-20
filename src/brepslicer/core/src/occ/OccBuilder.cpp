#include <brepslicer/Algo.h>
#include <occ/OccShape.h>

#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_Transform.hxx>
#include <BRepPrimAPI_MakeBox.hxx>
#include <BRepPrimAPI_MakeCone.hxx>
#include <BRepPrimAPI_MakeCylinder.hxx>
#include <BRepPrimAPI_MakeTorus.hxx>
#include <BRep_Builder.hxx>
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <Geom_BSplineSurface.hxx>
#include <Precision.hxx>
#include <TColgp_Array2OfPnt.hxx>
#include <TopoDS_Compound.hxx>
#include <gp_Ax1.hxx>
#include <gp_Dir.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>

#include <cmath>
#include <stdexcept>

namespace brepslicer {

class OccShapeBuilder : public IShapeBuilder {
public:
    std::shared_ptr<IShape> makeBox(double dx, double dy, double dz) override {
        return std::make_shared<OccShape>(BRepPrimAPI_MakeBox(dx, dy, dz).Shape());
    }
    std::shared_ptr<IShape> makeCylinder(double radius, double height) override {
        return std::make_shared<OccShape>(BRepPrimAPI_MakeCylinder(radius, height).Shape());
    }
    std::shared_ptr<IShape> makeCone(double r1, double r2, double height) override {
        return std::make_shared<OccShape>(BRepPrimAPI_MakeCone(r1, r2, height).Shape());
    }
    std::shared_ptr<IShape> makeTorus(double major_radius, double minor_radius) override {
        return std::make_shared<OccShape>(BRepPrimAPI_MakeTorus(major_radius, minor_radius).Shape());
    }
    std::shared_ptr<IShape> makeBumpFace(double dx, double dy, double bump_height) override {
        const int n = 12;
        TColgp_Array2OfPnt pts(1, n, 1, n);
        for (int i = 1; i <= n; ++i) {
            for (int j = 1; j <= n; ++j) {
                const double x = dx * (i - 1.0) / (n - 1);
                const double y = dy * (j - 1.0) / (n - 1);
                const double ui = (i - 1.0) / (n - 1) * 2.0 - 1.0;
                const double vj = (j - 1.0) / (n - 1) * 2.0 - 1.0;
                const double z = bump_height * std::exp(-3.5 * (ui * ui + vj * vj));
                pts.SetValue(i, j, gp_Pnt(x, y, z));
            }
        }
        GeomAPI_PointsToBSplineSurface fitter;
        fitter.Interpolate(pts);
        const Handle(Geom_BSplineSurface) surf = fitter.Surface();
        return std::make_shared<OccShape>(BRepBuilderAPI_MakeFace(surf, Precision::Confusion()).Face());
    }
    std::shared_ptr<IShape> transform(const std::shared_ptr<IShape>& shape, const Vec3& axis_origin,
                                      const Vec3& axis_dir, double angle_rad,
                                      const Vec3& translation) override {
        if (!shape) throw std::runtime_error("null shape");
        gp_Trsf R;
        R.SetRotation(gp_Ax1(gp_Pnt(axis_origin.x, axis_origin.y, axis_origin.z),
                             gp_Dir(axis_dir.x, axis_dir.y, axis_dir.z)),
                      angle_rad);
        gp_Trsf T;
        T.SetTranslation(gp_Vec(translation.x, translation.y, translation.z));
        T.Multiply(R);
        const TopoDS_Shape out =
            BRepBuilderAPI_Transform(occShape(*shape), T, Standard_True).Shape();
        return std::make_shared<OccShape>(out);
    }
    std::shared_ptr<IShape> compound(const std::vector<std::shared_ptr<IShape>>& parts) override {
        TopoDS_Compound comp;
        BRep_Builder builder;
        builder.MakeCompound(comp);
        for (const auto& p : parts) {
            if (p) builder.Add(comp, occShape(*p));
        }
        return std::make_shared<OccShape>(comp);
    }
};

IShapeBuilder* makeOccBuilder() { return new OccShapeBuilder(); }

}  // namespace brepslicer
