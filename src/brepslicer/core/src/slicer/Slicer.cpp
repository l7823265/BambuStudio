#include <brepslicer/Slicer.h>

#include <contour/ContourAssembler.h>
#include <geom/GeomUtil.h>
#include <intersect/AnalyticIntersect.h>
#include <topo/FaceIndex.h>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>

namespace brepslicer {

std::shared_ptr<IShape> readStep(const std::string& path) {
    return ShapeEngine::Kernel().readStep(path);
}

SliceResult sliceShape(const std::shared_ptr<IShape>& shape, const SliceOptions& opt) {
    SliceResult result;
    result.tolerance = opt.tolerance;
    result.normal = opt.normal;

    if (!shape || shape->empty()) {
        result.logs.push_back("empty shape");
        return result;
    }

    const SliceFrame frame = makeSliceFrame(opt.normal);

    FaceIndex index;
    index.build(ShapeEngine::Kernel().exploreFaces(*shape), shape->bbox());

    double dmin = 0, dmax = 0;
    index.shapeRange(opt.normal, dmin, dmax);

    std::vector<double> heights = opt.explicit_heights;
    if (heights.empty()) {
        double start = opt.start_height;
        if (!opt.start_height_set) {
            start = dmin + 0.5 * opt.layer_height;
        }
        int count = opt.layer_count;
        if (count < 0) {
            if (opt.layer_height <= 0) {
                throw std::runtime_error("layer-height must be positive");
            }
            count = static_cast<int>(std::floor((dmax - start) / opt.layer_height)) + 1;
        }
        heights.reserve(static_cast<size_t>(std::max(count, 0)));
        for (int i = 0; i < count; ++i) {
            const double h = start + i * opt.layer_height;
            if (h > dmax + opt.tolerance) break;
            if (h < dmin - opt.tolerance) continue;
            heights.push_back(h);
        }
    }

    int total_bridge = 0;
    double max_bridge = 0;
    for (double h : heights) {
        const Plane pln = makePlane(opt.normal, h);
        auto candidates = index.query(opt.normal, h, opt.tolerance);
        std::vector<RawSegment> segs;
        for (const FaceRecord* f : candidates) {
            auto part = intersectFaceWithPlane(*f, pln, frame, opt);
            for (RawSegment& s : part) {
                if (s.degenerate) {
                    TopologyEvent ev;
                    ev.z = h;
                    ev.event = s.degen_event.empty() ? "tangent_point" : s.degen_event;
                    ev.point = s.degen_point;
                    ev.has_point = true;
                    result.topology_events.push_back(ev);
                } else {
                    segs.push_back(std::move(s));
                }
            }
        }
        AssembleStats stats;
        Layer layer = assembleLayer(h, std::move(segs), frame, opt, stats);
        if (stats.bridged > 0) {
            total_bridge += stats.bridged;
            max_bridge = std::max(max_bridge, stats.max_bridge);
            std::ostringstream os;
            os << "z=" << h << " bridged " << stats.bridged
               << " gap(s), max=" << stats.max_bridge;
            result.logs.push_back(os.str());
        }
        if (stats.open_leftover > 0) {
            std::ostringstream os;
            os << "z=" << h << " left " << stats.open_leftover << " open chain(s)";
            result.logs.push_back(os.str());
        }
        result.layers.push_back(std::move(layer));
    }

    if (total_bridge > 0) {
        std::ostringstream os;
        os << "stitched " << total_bridge << " open gap(s), max bridge distance " << max_bridge;
        result.logs.push_back(os.str());
    }

    if (!opt.svg_dir.empty()) {
        writeAllLayerSvg(result, opt.svg_dir);
    }
    if (!opt.dxf_path.empty()) {
        writeAllLayerDxf(result, opt.dxf_path);
    }
    return result;
}

SliceResult sliceFile(const std::string& stepPath, const SliceOptions& opt) {
    return sliceShape(readStep(stepPath), opt);
}

}  // namespace brepslicer
