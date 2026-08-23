#include <brepslicer/Algo.h>
#include <occ/OccShape.h>

#include <BRepBuilderAPI_GTransform.hxx>
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
#include <gp_GTrsf.hxx>
#include <gp_Mat.hxx>
#include <gp_Pnt.hxx>
#include <gp_Trsf.hxx>
#include <gp_Vec.hxx>
#include <gp_XYZ.hxx>

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
    std::shared_ptr<IShape> transformAffine(const std::shared_ptr<IShape>& shape, const double R[9],
                                            const double t[3]) override {
        if (!shape) throw std::runtime_error("null shape");
        // Prefer gp_Trsf when the map is a similarity (uniform scale + orthogonal).
        // That preserves Plane/Cylinder/... so AnalyticIntersect stays accurate.
        const gp_XYZ c0(R[0], R[3], R[6]);
        const gp_XYZ c1(R[1], R[4], R[7]);
        const gp_XYZ c2(R[2], R[5], R[8]);
        const double s0 = c0.Modulus();
        const double s1 = c1.Modulus();
        const double s2 = c2.Modulus();
        const double smean = (s0 + s1 + s2) / 3.0;
        const bool   scale_ok =
            smean > 1e-12 && std::abs(s0 - smean) <= 1e-6 * smean &&
            std::abs(s1 - smean) <= 1e-6 * smean && std::abs(s2 - smean) <= 1e-6 * smean;
        const double n01 = c0.Crossed(c1).Modulus();
        const double n12 = c1.Crossed(c2).Modulus();
        const double n20 = c2.Crossed(c0).Modulus();
        const bool   ortho_ok =
            scale_ok && std::abs(n01 - s0 * s1) <= 1e-5 * smean * smean &&
            std::abs(n12 - s1 * s2) <= 1e-5 * smean * smean &&
            std::abs(n20 - s2 * s0) <= 1e-5 * smean * smean;
        if (ortho_ok) {
            try {
                gp_Trsf tr;
                tr.SetValues(R[0], R[1], R[2], t[0], R[3], R[4], R[5], t[1], R[6], R[7], R[8], t[2]);
                const TopoDS_Shape out =
                    BRepBuilderAPI_Transform(occShape(*shape), tr, Standard_True).Shape();
                return std::make_shared<OccShape>(out);
            } catch (...) {
                // Fall through to GTrsf.
            }
        }
        gp_Mat mat(R[0], R[1], R[2], R[3], R[4], R[5], R[6], R[7], R[8]);
        gp_GTrsf gtr;
        gtr.SetVectorialPart(mat);
        gtr.SetTranslationPart(gp_XYZ(t[0], t[1], t[2]));
        const TopoDS_Shape out =
            BRepBuilderAPI_GTransform(occShape(*shape), gtr, Standard_True).Shape();
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
