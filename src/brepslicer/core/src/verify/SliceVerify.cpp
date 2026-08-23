#include <brepslicer/Slicer.h>
#include <contour/ContourAssembler.h>
#include <geom/GeomUtil.h>
#include <intersect/BSplineFit.h>
#include <verify/SliceVerify.h>

#include <algorithm>
#include <cmath>
#include <sstream>

namespace brepslicer {

namespace {

Vec3 segStart(const Segment& s) {
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) return s.ctrl_pts.front();
    return s.start;
}
Vec3 segEnd(const Segment& s) {
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) return s.ctrl_pts.back();
    return s.end;
}

}  // namespace

double contourLength(const Contour& c) {
    double L = 0;
    for (const Segment& s : c.segments) {
        if (s.type == SegmentType::Arc) {
            L += std::abs(s.radius * s.sweep);
        } else if (s.type == SegmentType::Ellipse) {
            const double a = s.radius;
            const double b = s.radius_b;
            const double hh = ((a - b) / (a + b)) * ((a - b) / (a + b));
            const double full = kPi * (a + b) * (1.0 + 3.0 * hh / (10.0 + std::sqrt(4.0 - 3.0 * hh)));
            L += full * std::abs(s.sweep) / kTwoPi;
        } else if (s.type == SegmentType::BSpline) {
            L += bsplineLength(s);
        } else {
            L += dist(s.start, s.end);
        }
    }
    return L;
}

double layerPerimeter(const Layer& layer) {
    double L = 0;
    for (const Contour& c : layer.contours) L += contourLength(c);
    return L;
}

double layerArea(const Layer& layer, const SliceFrame& frame) {
    double a = 0;
    for (const Contour& c : layer.contours) a += contourSignedArea(c, frame);
    return a;
}

VerifyReport verifySlice(const SliceResult& result) {
    VerifyReport r;
    const SliceFrame frame = makeSliceFrame(result.normal);
    for (size_t li = 0; li < result.layers.size(); ++li) {
        const Layer& layer = result.layers[li];
        for (size_t ci = 0; ci < layer.contours.size(); ++ci) {
            const Contour& c = layer.contours[ci];
            if (c.segments.empty()) {
                r.errors.push_back("empty contour");
                r.ok = false;
                continue;
            }
            if (c.closed) {
                const double gap = dist(segStart(c.segments.front()), segEnd(c.segments.back()));
                const double tol = std::max(result.tolerance, result.stitch_tolerance);
                if (gap > tol) {
                    std::ostringstream os;
                    os << "layer " << li << " contour " << ci << " not closed, gap=" << gap;
                    r.errors.push_back(os.str());
                    r.ok = false;
                }
                for (size_t i = 0; i + 1 < c.segments.size(); ++i) {
                    const double g = dist(segEnd(c.segments[i]), segStart(c.segments[i + 1]));
                    if (g > tol) {
                        std::ostringstream os;
                        os << "layer " << li << " contour " << ci << " broken at seg " << i;
                        r.errors.push_back(os.str());
                        r.ok = false;
                    }
                }
            }
            const double a = contourSignedArea(c, frame);
            if (c.orientation == "outer" && a < -1e-12) {
                r.errors.push_back("outer contour is CW");
                r.ok = false;
            }
            if (c.orientation == "inner" && a > 1e-12) {
                r.errors.push_back("inner contour is CCW");
                r.ok = false;
            }
            if (c.parent) {
                if (*c.parent < 0 || *c.parent >= static_cast<int>(layer.contours.size())) {
                    r.errors.push_back("parent index out of range");
                    r.ok = false;
                }
            }
        }
    }

    if (result.layers.size() >= 2) {
        std::vector<double> areas;
        areas.reserve(result.layers.size());
        for (const Layer& layer : result.layers) areas.push_back(layerArea(layer, frame));
        for (size_t i = 1; i < areas.size(); ++i) {
            const double a0 = std::abs(areas[i - 1]);
            const double a1 = std::abs(areas[i]);
            const double denom = std::max(1e-12, 0.5 * (a0 + a1));
            const double rel = std::abs(a1 - a0) / denom;
            if (rel > 0.25) {
                std::ostringstream os;
                os << "area jump " << rel << " between z=" << result.layers[i - 1].z << " and z="
                   << result.layers[i].z << " (expect a topology_event)";
                if (result.topology_events.empty()) r.warnings.push_back(os.str());
                else r.notes.push_back(os.str());
            }
        }
    }
    return r;
}

VerifyReport verifyAgainstRef(const SliceResult& result, const std::shared_ptr<IShape>& shape) {
    VerifyReport r = verifySlice(result);
    if (!shape) {
        r.errors.push_back("null shape");
        r.ok = false;
        return r;
    }
    for (const Layer& layer : result.layers) {
        bool degenerate_layer = false;
        for (const TopologyEvent& e : result.topology_events) {
            if (std::abs(e.z - layer.z) <= result.tolerance) degenerate_layer = true;
        }
        bool coplanar_layer = false;
        for (const Contour& c : layer.contours) {
            if (c.coplanar) coplanar_layer = true;
        }
        try {
            const Plane pln = makePlane(result.normal, layer.z);
            const double refLen = ShapeEngine::SectionRef().sectionLength(*shape, pln);
            const double ourLen = layerPerimeter(layer);
            const double denom = std::max(1e-12, 0.5 * (ourLen + refLen));
            const double rel = std::abs(ourLen - refLen) / denom;
            std::ostringstream os;
            os << "z=" << layer.z << " perimeter our=" << ourLen << " occt=" << refLen
               << " rel=" << rel;
            r.notes.push_back(os.str());
            if (degenerate_layer) {
                r.notes.push_back("degenerate layer: skip relative-error fail vs OCCT");
                continue;
            }
            if (coplanar_layer && refLen < 1e-9) {
                r.notes.push_back("coplanar layer: OCCT section empty, skip rel");
                continue;
            }
            if (rel > 1e-4) {
                r.errors.push_back(os.str());
                r.ok = false;
            }
        } catch (...) {
            r.warnings.push_back("OCCT section failed on a layer");
        }
    }
    return r;
}

}  // namespace brepslicer
