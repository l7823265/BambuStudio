#include <occ/OccTopology.h>
#include <occ/OccShape.h>

#include <BRep_Builder.hxx>
#include <TopExp.hxx>
#include <TopExp_Explorer.hxx>
#include <TopTools_IndexedDataMapOfShapeListOfShape.hxx>
#include <TopTools_IndexedMapOfShape.hxx>
#include <TopTools_ListIteratorOfListOfShape.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Compound.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shape.hxx>

#include <algorithm>
#include <map>

namespace brepslicer {

SolidAdjacency buildSolidAdjacency(const std::vector<FaceRecord>& faces) {
    SolidAdjacency adj;
    if (faces.empty()) return adj;

    TopTools_IndexedMapOfShape face_map;
    std::map<int, int> occ_face_index_to_id;  // 1-based OCC index -> FaceRecord.face_id
    int max_face_id = -1;
    for (const FaceRecord& fr : faces) {
        if (!fr.face) continue;
        const auto* of = dynamic_cast<const OccFace*>(fr.face.get());
        if (!of) continue;
        const int idx = face_map.Add(of->occFace());
        occ_face_index_to_id[idx] = fr.face_id;
        max_face_id = std::max(max_face_id, fr.face_id);
    }
    if (face_map.Extent() == 0) return adj;

    // Compound of all faces so EDGE→FACE ancestors cover the whole model.
    TopoDS_Compound comp;
    BRep_Builder builder;
    builder.MakeCompound(comp);
    for (int i = 1; i <= face_map.Extent(); ++i) builder.Add(comp, face_map(i));

    TopTools_IndexedDataMapOfShapeListOfShape edge_to_faces;
    TopExp::MapShapesAndAncestors(comp, TopAbs_EDGE, TopAbs_FACE, edge_to_faces);

    adj.edge_faces.resize(static_cast<size_t>(edge_to_faces.Extent()));
    adj.face_links.assign(static_cast<size_t>(std::max(0, max_face_id + 1)), {});

    for (int ei = 1; ei <= edge_to_faces.Extent(); ++ei) {
        const int edge_id = ei - 1;
        std::vector<int> fids;
        for (TopTools_ListIteratorOfListOfShape it(edge_to_faces(ei)); it.More(); it.Next()) {
            const TopoDS_Face f = TopoDS::Face(it.Value());
            const int loc = face_map.FindIndex(f);
            if (loc <= 0) continue;
            const auto jt = occ_face_index_to_id.find(loc);
            if (jt == occ_face_index_to_id.end()) continue;
            fids.push_back(jt->second);
        }
        // Unique face ids on this edge.
        std::sort(fids.begin(), fids.end());
        fids.erase(std::unique(fids.begin(), fids.end()), fids.end());
        adj.edge_faces[static_cast<size_t>(edge_id)] = fids;

        if (fids.size() == 2) {
            adj.face_links[static_cast<size_t>(fids[0])].push_back({edge_id, fids[1]});
            adj.face_links[static_cast<size_t>(fids[1])].push_back({edge_id, fids[0]});
        } else if (fids.size() == 1) {
            adj.face_links[static_cast<size_t>(fids[0])].push_back({edge_id, -1});
        } else {
            // Non-manifold: link every pair.
            for (size_t a = 0; a < fids.size(); ++a) {
                for (size_t b = 0; b < fids.size(); ++b) {
                    if (a == b) continue;
                    adj.face_links[static_cast<size_t>(fids[a])].push_back(
                        {edge_id, fids[b]});
                }
            }
        }
    }
    return adj;
}

}  // namespace brepslicer
