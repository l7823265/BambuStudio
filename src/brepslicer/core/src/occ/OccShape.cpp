#include <occ/OccShape.h>

#include <BRepAdaptor_Curve.hxx>
#include <BRepAdaptor_Surface.hxx>
#include <BRepBndLib.hxx>
#include <BRepClass_FaceClassifier.hxx>
#include <BRepTools_WireExplorer.hxx>
#include <BRep_Tool.hxx>
#include <Bnd_Box.hxx>
#include <Extrema_ExtPS.hxx>
#include <GeomAbs_CurveType.hxx>
#include <GeomAbs_SurfaceType.hxx>
#include <Standard_Failure.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs.hxx>
#include <TopoDS.hxx>
#include <gp_Circ.hxx>
#include <gp_Cone.hxx>
#include <gp_Cylinder.hxx>
#include <gp_Dir.hxx>
#include <gp_Elips.hxx>
#include <gp_Lin.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>
#include <gp_Pnt2d.hxx>
#include <gp_Sphere.hxx>
#include <gp_Torus.hxx>
#include <gp_Vec.hxx>
#include <TopoDS_Wire.hxx>

#include <cmath>
#include <memory>
#include <stdexcept>

namespace brepslicer {
namespace {

Vec3 v3(const gp_Pnt& p) { return {p.X(), p.Y(), p.Z()}; }
Vec3 v3(const gp_Dir& d) { return {d.X(), d.Y(), d.Z()}; }

BBox boxOf(const TopoDS_Shape& s) {
    BBox b;
    Bnd_Box occ;
    BRepBndLib::Add(s, occ);
    if (occ.IsVoid()) return b;
    occ.Get(b.xmin, b.ymin, b.zmin, b.xmax, b.ymax, b.zmax);
    b.valid = true;
    return b;
}

}  // namespace

OccShape::OccShape(const TopoDS_Shape& s) : shape_(s) {}

IShape::ShapeType OccShape::type() const {
    switch (shape_.ShapeType()) {
    case TopAbs_COMPOUND:
    case TopAbs_COMPSOLID:
        return ShapeType::Compound;
    case TopAbs_SOLID:
        return ShapeType::Solid;
    case TopAbs_SHELL:
        return ShapeType::Shell;
    case TopAbs_FACE:
        return ShapeType::Face;
    case TopAbs_WIRE:
        return ShapeType::Wire;
    case TopAbs_EDGE:
        return ShapeType::Edge;
    case TopAbs_VERTEX:
        return ShapeType::Vertex;
    default:
        return ShapeType::Other;
    }
}

bool OccShape::empty() const { return shape_.IsNull(); }
BBox OccShape::bbox() const { return boxOf(shape_); }

OccEdge::OccEdge(const TopoDS_Edge& e) : edge_(e) {}
bool OccEdge::empty() const { return edge_.IsNull(); }
BBox OccEdge::bbox() const { return boxOf(edge_); }

CurveData OccEdge::curve() const {
    CurveData c;
    c.degenerated = BRep_Tool::Degenerated(edge_);
    c.reversed = (edge_.Orientation() == TopAbs_REVERSED);
    if (c.degenerated) return c;
    BRepAdaptor_Curve ac(edge_);
    c.first = ac.FirstParameter();
    c.last = ac.LastParameter();
    switch (ac.GetType()) {
    case GeomAbs_Line: {
        c.kind = CurveKind::Line;
        const gp_Lin ln = ac.Line();
        c.line = {v3(ln.Location()), v3(ln.Direction())};
        break;
    }
    case GeomAbs_Circle: {
        c.kind = CurveKind::Circle;
        const gp_Circ ci = ac.Circle();
        c.circle = {v3(ci.Location()), v3(ci.Axis().Direction()), v3(ci.XAxis().Direction()),
                    ci.Radius()};
        break;
    }
    case GeomAbs_Ellipse: {
        c.kind = CurveKind::Ellipse;
        const gp_Elips el = ac.Ellipse();
        c.ellipse = {v3(el.Location()), v3(el.Axis().Direction()), v3(el.XAxis().Direction()),
                     el.MajorRadius(), el.MinorRadius()};
        break;
    }
    default:
        c.kind = CurveKind::Other;
        break;
    }
    return c;
}

Vec3 OccEdge::eval(double t) const {
    BRepAdaptor_Curve ac(edge_);
    return v3(ac.Value(t));
}

OccFace::OccFace(const TopoDS_Face& f) : face_(f) {}
bool OccFace::empty() const { return face_.IsNull(); }
BBox OccFace::bbox() const { return boxOf(face_); }

SurfaceData OccFace::surface() const {
    SurfaceData s;
    BRepAdaptor_Surface ads(face_, Standard_True);
    switch (ads.GetType()) {
    case GeomAbs_Plane: {
        s.kind = SurfaceKind::Plane;
        const gp_Pln pl = ads.Plane();
        const Vec3 n = v3(pl.Axis().Direction());
        s.plane = {n, dot(n, v3(pl.Location()))};
        break;
    }
    case GeomAbs_Cylinder: {
        s.kind = SurfaceKind::Cylinder;
        const gp_Cylinder cy = ads.Cylinder();
        s.cylinder = {v3(cy.Location()), v3(cy.Axis().Direction()), cy.Radius()};
        break;
    }
    case GeomAbs_Sphere: {
        s.kind = SurfaceKind::Sphere;
        const gp_Sphere sp = ads.Sphere();
        s.sphere = {v3(sp.Location()), sp.Radius()};
        break;
    }
    case GeomAbs_Cone: {
        s.kind = SurfaceKind::Cone;
        const gp_Cone co = ads.Cone();
        s.cone = {v3(co.Apex()), v3(co.Axis().Direction()), std::abs(co.SemiAngle())};
        break;
    }
    case GeomAbs_Torus: {
        s.kind = SurfaceKind::Torus;
        const gp_Torus to = ads.Torus();
        s.torus = {v3(to.Location()), v3(to.Axis().Direction()), to.MajorRadius(), to.MinorRadius()};
        break;
    }
    default:
        s.kind = SurfaceKind::Other;
        break;
    }
    return s;
}

PointClass OccFace::classify(const Vec3& p, double tol) const {
    try {
        BRepClass_FaceClassifier fc(face_, gp_Pnt(p.x, p.y, p.z), tol);
        switch (fc.State()) {
        case TopAbs_IN:
            return PointClass::In;
        case TopAbs_ON:
            return PointClass::On;
        default:
            return PointClass::Out;
        }
    } catch (const Standard_Failure&) {
        return PointClass::Out;
    }
}

UVBox OccFace::uvDomain() const {
    UVBox b;
    BRepAdaptor_Surface ads(face_, Standard_True);
    b.umin = ads.FirstUParameter();
    b.umax = ads.LastUParameter();
    b.vmin = ads.FirstVParameter();
    b.vmax = ads.LastVParameter();
    b.periodic_u = ads.IsUPeriodic() == Standard_True;
    b.periodic_v = ads.IsVPeriodic() == Standard_True;
    b.period_u = b.periodic_u ? ads.UPeriod() : 0;
    b.period_v = b.periodic_v ? ads.VPeriod() : 0;
    return b;
}

Vec3 OccFace::evalUV(double u, double v) const {
    BRepAdaptor_Surface ads(face_, Standard_True);
    return v3(ads.Value(u, v));
}

bool OccFace::derivUV(double u, double v, Vec3& Su, Vec3& Sv) const {
    try {
        BRepAdaptor_Surface ads(face_, Standard_True);
        gp_Pnt p;
        gp_Vec du, dv;
        ads.D1(u, v, p, du, dv);
        Su = {du.X(), du.Y(), du.Z()};
        Sv = {dv.X(), dv.Y(), dv.Z()};
        return true;
    } catch (const Standard_Failure&) {
        return false;
    }
}

bool OccFace::invertUV(const Vec3& p, double& u, double& v, double tol) const {
    try {
        BRepAdaptor_Surface ads(face_, Standard_True);
        Extrema_ExtPS ext(gp_Pnt(p.x, p.y, p.z), ads, ads.FirstUParameter(), ads.LastUParameter(),
                          ads.FirstVParameter(), ads.LastVParameter(), 1e-12, 1e-12);
        if (!ext.IsDone() || ext.NbExt() < 1) return false;
        double best = 1e300;
        int ibest = 1;
        for (int i = 1; i <= ext.NbExt(); ++i) {
            if (ext.SquareDistance(i) < best) {
                best = ext.SquareDistance(i);
                ibest = i;
            }
        }
        ext.Point(ibest).Parameter(u, v);
        return std::sqrt(best) <= std::max(10.0 * tol, 1e-6);
    } catch (const Standard_Failure&) {
        return false;
    }
}

PointClass OccFace::classifyUV(double u, double v, double tol) const {
    try {
        BRepClass_FaceClassifier fc(face_, gp_Pnt2d(u, v), tol);
        switch (fc.State()) {
        case TopAbs_IN:
            return PointClass::In;
        case TopAbs_ON:
            return PointClass::On;
        default:
            return PointClass::Out;
        }
    } catch (const Standard_Failure&) {
        return PointClass::Out;
    }
}

std::vector<std::shared_ptr<IEdge>> OccFace::edges() const {
    std::vector<std::shared_ptr<IEdge>> out;
    for (TopExp_Explorer ex(face_, TopAbs_EDGE); ex.More(); ex.Next()) {
        out.push_back(std::make_shared<OccEdge>(TopoDS::Edge(ex.Current())));
    }
    return out;
}

std::vector<std::vector<std::shared_ptr<IEdge>>> OccFace::wires() const {
    std::vector<std::vector<std::shared_ptr<IEdge>>> out;
    for (TopExp_Explorer wr(face_, TopAbs_WIRE); wr.More(); wr.Next()) {
        std::vector<std::shared_ptr<IEdge>> wire;
        const TopoDS_Wire w = TopoDS::Wire(wr.Current());
        for (BRepTools_WireExplorer we(w, face_); we.More(); we.Next()) {
            wire.push_back(std::make_shared<OccEdge>(we.Current()));
        }
        out.push_back(std::move(wire));
    }
    return out;
}

const TopoDS_Shape& occShape(const IShape& s) {
    if (const auto* o = dynamic_cast<const OccShape*>(&s)) return o->occ();
    if (const auto* f = dynamic_cast<const OccFace*>(&s)) return f->occFace();
    if (const auto* e = dynamic_cast<const OccEdge*>(&s)) return e->occEdge();
    throw std::runtime_error("not an OCC shape");
}

}  // namespace brepslicer
