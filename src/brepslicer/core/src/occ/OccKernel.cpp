#include <brepslicer/Algo.h>
#include <occ/OccShape.h>

#include <IFSelect_ReturnStatus.hxx>
#include <STEPControl_Reader.hxx>
#include <TopExp_Explorer.hxx>
#include <TopAbs.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Shell.hxx>
#include <TopoDS_Solid.hxx>

#include <stdexcept>

namespace brepslicer {

class OccKernel : public IBrepKernel {
public:
    std::shared_ptr<IShape> readStep(const std::string& path) override {
        STEPControl_Reader reader;
        if (reader.ReadFile(path.c_str()) != IFSelect_RetDone) {
            throw std::runtime_error("Failed to read STEP file: " + path);
        }
        reader.TransferRoots();
        TopoDS_Shape shape = reader.OneShape();
        if (shape.IsNull()) throw std::runtime_error("STEP file produced an empty shape: " + path);
        return std::make_shared<OccShape>(shape);
    }

    std::vector<FaceRecord> exploreFaces(const IShape& shape) override {
        const TopoDS_Shape& root = occShape(shape);
        std::vector<FaceRecord> faces;
        int face_id = 0;
        int solid_count = 0;
        int shell_count = 0;
        bool has_solid = false;

        auto addFace = [&](const TopoDS_Face& f, int sid, int shid) {
            auto face = std::make_shared<OccFace>(f);
            FaceRecord rec;
            rec.face = face;
            rec.solid_id = sid;
            rec.shell_id = shid;
            rec.face_id = face_id++;
            rec.box = face->bbox();
            faces.push_back(std::move(rec));
        };

        for (TopExp_Explorer sol(root, TopAbs_SOLID); sol.More(); sol.Next()) {
            has_solid = true;
            const int sid = solid_count++;
            const TopoDS_Solid solid = TopoDS::Solid(sol.Current());
            for (TopExp_Explorer sh(solid, TopAbs_SHELL); sh.More(); sh.Next()) {
                const int shid = shell_count++;
                const TopoDS_Shell shell = TopoDS::Shell(sh.Current());
                for (TopExp_Explorer fa(shell, TopAbs_FACE); fa.More(); fa.Next()) {
                    addFace(TopoDS::Face(fa.Current()), sid, shid);
                }
            }
        }
        if (has_solid) return faces;

        for (TopExp_Explorer sh(root, TopAbs_SHELL, TopAbs_SOLID); sh.More(); sh.Next()) {
            const int shid = shell_count++;
            const TopoDS_Shell shell = TopoDS::Shell(sh.Current());
            for (TopExp_Explorer fa(shell, TopAbs_FACE); fa.More(); fa.Next()) {
                addFace(TopoDS::Face(fa.Current()), 0, shid);
            }
        }
        if (!faces.empty()) return faces;

        for (TopExp_Explorer fa(root, TopAbs_FACE); fa.More(); fa.Next()) {
            addFace(TopoDS::Face(fa.Current()), 0, 0);
        }
        return faces;
    }
};

IBrepKernel* makeOccKernel() { return new OccKernel(); }

}  // namespace brepslicer
