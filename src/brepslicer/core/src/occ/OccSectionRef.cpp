#include <brepslicer/Algo.h>
#include <geom/GeomUtil.h>
#include <intersect/AnalyticIntersect.h>
#include <intersect/BSplineFit.h>
#include <occ/OccShape.h>

#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <BRepExtrema_DistShapeShape.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>

#include <algorithm>
#include <cmath>
#include <vector>

namespace brepslicer {
namespace {

TopoDS_Shape sectionShape(const IShape& shape, const Plane& plane) {
    const gp_Pln pln(gp_Pnt(plane.n.x * plane.d, plane.n.y * plane.d, plane.n.z * plane.d),
                     gp_Dir(plane.n.x, plane.n.y, plane.n.z));
    const TopoDS_Face pface = BRepBuilderAPI_MakeFace(pln).Face();
    BRepAlgoAPI_Section section(occShape(shape), pface, Standard_False);
    section.Approximation(Standard_False);
    section.Build();
    if (!section.IsDone()) return {};
    return section.Shape();
}

Vec3 arcAt(const Segment& s, double t) {
    const double a = s.start_angle + t * s.sweep;
    Vec3 x = cross({0, 0, 1}, s.normal);
    if (length(x) < 1e-12) x = cross({0, 1, 0}, s.normal);
    x = normalized(x);
    const Vec3 y = cross(s.normal, x);
    return s.center + x * (s.radius * std::cos(a)) + y * (s.radius * std::sin(a));
}

void sampleSegment(const Segment& s, std::vector<Vec3>& pts) {
    pts.clear();
    if (s.type == SegmentType::BSpline) {
        for (int i = 0; i <= 4; ++i) pts.push_back(bsplineEval(s, static_cast<double>(i) / 4.0));
        return;
    }
    if (s.type == SegmentType::Arc || s.type == SegmentType::Ellipse) {
        for (int i = 0; i <= 4; ++i) {
            if (s.type == SegmentType::Arc) {
                pts.push_back(arcAt(s, static_cast<double>(i) / 4.0));
            } else {
                const Vec3 maj = normalized(s.major_axis);
                const Vec3 minv = cross(s.normal, maj);
                const double ang = s.start_angle + (static_cast<double>(i) / 4.0) * s.sweep;
                pts.push_back(s.center + maj * (s.radius * std::cos(ang)) +
                              minv * (s.radius_b * std::sin(ang)));
            }
        }
        return;
    }
    pts.push_back(s.start);
    pts.push_back(s.end);
}

double distToSection(const TopoDS_Shape& section, const Vec3& p) {
    BRepExtrema_DistShapeShape dist(BRepBuilderAPI_MakeVertex(gp_Pnt(p.x, p.y, p.z)).Vertex(),
                                    section);
    dist.Perform();
    return dist.IsDone() ? dist.Value() : 1e100;
}

}  // namespace

class OccSectionRef : public ISectionRef {
public:
    double sectionLength(const IShape& shape, const Plane& plane) override {
        const TopoDS_Shape section = sectionShape(shape, plane);
        double len = 0;
        for (TopExp_Explorer ex(section, TopAbs_EDGE); ex.More(); ex.Next()) {
            GProp_GProps props;
            BRepGProp::LinearProperties(ex.Current(), props);
            len += props.Mass();
        }
        return len;
    }

    void clipSegmentsToSection(const IShape& shape, const Plane& plane,
                               std::vector<RawSegment>& segs, double tol) override {
        const TopoDS_Shape section = sectionShape(shape, plane);
        if (section.IsNull()) return;
        const double snap = std::max(tol * 100.0, 0.01);
        std::vector<Vec3> samples;
        segs.erase(std::remove_if(segs.begin(), segs.end(),
                                  [&](const RawSegment& rs) {
                                      if (rs.degenerate) return false;
                                      sampleSegment(rs.geom, samples);
                                      for (const Vec3& p : samples) {
                                          if (distToSection(section, p) > snap) return true;
                                      }
                                      return false;
                                  }),
                   segs.end());
    }
};

ISectionRef* makeOccSectionRef() { return new OccSectionRef(); }

}  // namespace brepslicer
