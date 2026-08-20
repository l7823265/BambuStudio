#include <brepslicer/Algo.h>
#include <occ/OccShape.h>

#include <BRepAlgoAPI_Section.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepGProp.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS_Face.hxx>
#include <gp_Dir.hxx>
#include <gp_Pln.hxx>
#include <gp_Pnt.hxx>

namespace brepslicer {

class OccSectionRef : public ISectionRef {
public:
    double sectionLength(const IShape& shape, const Plane& plane) override {
        const gp_Pln pln(gp_Pnt(plane.n.x * plane.d, plane.n.y * plane.d, plane.n.z * plane.d),
                         gp_Dir(plane.n.x, plane.n.y, plane.n.z));
        const TopoDS_Face pface = BRepBuilderAPI_MakeFace(pln).Face();
        BRepAlgoAPI_Section section(occShape(shape), pface, Standard_False);
        section.Approximation(Standard_False);
        section.Build();
        double len = 0;
        if (!section.IsDone()) return 0;
        for (TopExp_Explorer ex(section.Shape(), TopAbs_EDGE); ex.More(); ex.Next()) {
            GProp_GProps props;
            BRepGProp::LinearProperties(ex.Current(), props);
            len += props.Mass();
        }
        return len;
    }
};

ISectionRef* makeOccSectionRef() { return new OccSectionRef(); }

}  // namespace brepslicer
