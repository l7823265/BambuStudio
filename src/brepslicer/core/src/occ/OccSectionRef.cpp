#include <brepslicer/Algo.h>
#include <geom/GeomUtil.h>
#include <intersect/AnalyticIntersect.h>
#include <intersect/BSplineFit.h>
#include <occ/OccShape.h>

#include <BRepAlgoAPI_Section.hxx>
#include <BRepAdaptor_Curve.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <BRepGProp.hxx>
#include <GCPnts_UniformAbscissa.hxx>
#include <GProp_GProps.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
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

// Sample section wires once; point queries use this cloud (avoids DistShapeShape per sample).
std::vector<Vec3> sampleSectionPoints(const TopoDS_Shape& section, double spacing) {
    std::vector<Vec3> pts;
    spacing = std::max(spacing, 1e-3);
    for (TopExp_Explorer ex(section, TopAbs_EDGE); ex.More(); ex.Next()) {
        try {
            BRepAdaptor_Curve c(TopoDS::Edge(ex.Current()));
            const double f = c.FirstParameter();
            const double l = c.LastParameter();
            if (l <= f) continue;
            GCPnts_UniformAbscissa discret(c, spacing, f, l);
            if (!discret.IsDone() || discret.NbPoints() < 2) {
                const gp_Pnt a = c.Value(f);
                const gp_Pnt b = c.Value(l);
                pts.push_back({a.X(), a.Y(), a.Z()});
                pts.push_back({b.X(), b.Y(), b.Z()});
                continue;
            }
            for (int i = 1; i <= discret.NbPoints(); ++i) {
                const gp_Pnt p = c.Value(discret.Parameter(i));
                pts.push_back({p.X(), p.Y(), p.Z()});
            }
        } catch (const Standard_Failure&) {
            continue;
        }
    }
    return pts;
}

double distToPointCloud(const std::vector<Vec3>& cloud, const Vec3& p) {
    double best = 1e100;
    for (const Vec3& q : cloud) {
        const double d = dist(p, q);
        if (d < best) best = d;
    }
    return best;
}

bool isAnalyticSegment(const Segment& s) {
    return s.type == SegmentType::Line || s.type == SegmentType::Arc ||
           s.type == SegmentType::Ellipse;
}

Segment subBSpline(const Segment& s, double t0, double t1) {
    t0 = std::min(1.0, std::max(0.0, t0));
    t1 = std::min(1.0, std::max(0.0, t1));
    if (t1 < t0) std::swap(t0, t1);
    const int n = std::max(4, static_cast<int>(std::ceil((t1 - t0) * 16)) + 1);
    std::vector<Vec3> pts;
    pts.reserve(static_cast<size_t>(n));
    for (int i = 0; i < n; ++i) {
        const double u = t0 + (t1 - t0) * (static_cast<double>(i) / (n - 1));
        pts.push_back(bsplineEval(s, u));
    }
    return fitCubicBSpline(pts, false, std::max(1e-4, s.fit_error));
}

// L2 B-splines only: keep parameter intervals near the OCC section wire.
void clipBSplineBySplitting(const RawSegment& rs, const std::vector<Vec3>& section_pts, double snap,
                            double minLen, std::vector<RawSegment>& out) {
    if (rs.degenerate) {
        out.push_back(rs);
        return;
    }
    const Segment& s = rs.geom;
    constexpr int n = 33;
    std::vector<char> inlier(static_cast<size_t>(n), 0);
    int n_in = 0;
    for (int i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / (n - 1);
        if (distToPointCloud(section_pts, bsplineEval(s, t)) <= snap) {
            inlier[static_cast<size_t>(i)] = 1;
            ++n_in;
        }
    }

    if (n_in == n) {
        out.push_back(rs);
        return;
    }

    for (int i = 1; i + 1 < n; ++i) {
        if (!inlier[static_cast<size_t>(i)] && inlier[static_cast<size_t>(i - 1)] &&
            inlier[static_cast<size_t>(i + 1)]) {
            inlier[static_cast<size_t>(i)] = 1;
            ++n_in;
        }
    }
    if (n_in == 0) return;

    int i = 0;
    while (i < n) {
        while (i < n && !inlier[static_cast<size_t>(i)]) ++i;
        if (i >= n) break;
        const int i0 = i;
        while (i < n && inlier[static_cast<size_t>(i)]) ++i;
        const int i1 = i - 1;
        if (i1 <= i0) continue;

        const double t0 = static_cast<double>(i0) / (n - 1);
        const double t1 = static_cast<double>(i1) / (n - 1);
        RawSegment piece = rs;
        piece.closed_loop = false;
        piece.geom = subBSpline(s, t0, t1);
        if (bsplineLength(piece.geom) < minLen) continue;
        out.push_back(std::move(piece));
    }
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
        // L1 line / arc / ellipse: trust face trim; never filter against OCC Section.
        bool need_section = false;
        for (const RawSegment& rs : segs) {
            if (!rs.degenerate && rs.geom.type == SegmentType::BSpline) {
                need_section = true;
                break;
            }
        }
        if (!need_section) return;

        const TopoDS_Shape section = sectionShape(shape, plane);
        if (section.IsNull()) return;
        const double snap = std::max(tol * 100.0, 0.01);
        const double minLen = std::max(tol, 1e-6);
        // Sample denser than snap so polyline distance ≈ true section distance.
        const double spacing = std::max(0.05, 0.5 * snap);
        const std::vector<Vec3> section_pts = sampleSectionPoints(section, spacing);
        if (section_pts.empty()) return;
        const double cloud_snap = snap + 0.5 * spacing;

        std::vector<RawSegment> clipped;
        clipped.reserve(segs.size() + 4);
        for (const RawSegment& rs : segs) {
            if (rs.degenerate || isAnalyticSegment(rs.geom)) {
                clipped.push_back(rs);
                continue;
            }
            clipBSplineBySplitting(rs, section_pts, cloud_snap, minLen, clipped);
        }
        segs = std::move(clipped);
    }
};

ISectionRef* makeOccSectionRef() { return new OccSectionRef(); }

}  // namespace brepslicer
