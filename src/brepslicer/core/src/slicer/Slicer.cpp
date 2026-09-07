#include <brepslicer/Slicer.h>

#include <contour/ContourAssembler.h>
#include <geom/GeomUtil.h>
#include <intersect/AnalyticIntersect.h>
#include <topo/FaceIndex.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <set>
#include <sstream>
#include <stdexcept>

namespace brepslicer {

namespace {

using SteadyClock = std::chrono::steady_clock;

double elapsedMs(const SteadyClock::time_point& t0) {
    return std::chrono::duration<double, std::milli>(SteadyClock::now() - t0).count();
}

void appendTimingLog(std::vector<std::string>& logs, const char* label, double ms, int layers = 0) {
    std::ostringstream os;
    os << std::fixed << std::setprecision(1);
    os << "timing: " << label << "=" << ms << "ms";
    if (layers > 0) os << " (" << layers << " layers, " << (ms / layers) << "ms/layer)";
    logs.push_back(os.str());
}

void appendFaceSeedLogs(double z, const std::vector<FaceSeedStats>& stats,
                        std::vector<std::string>& logs) {
    if (stats.empty()) return;
    std::ostringstream os;
    os << "z=" << z << " uvmatch face seeds (id:raw/kept/kept_chains/segs):";
    int zero_raw = 0;
    int zero_kept = 0;
    int kept_no_seg = 0;
    for (const FaceSeedStats& s : stats) {
        os << " " << s.face_id << ":" << s.raw << "/" << s.kept << "/" << s.kept_chains << "/"
           << s.segments << "/n" << s.neighbor_splits;
        if (s.uvmarch_fallback) os << "*";
        if (s.raw == 0) ++zero_raw;
        if (s.raw > 0 && s.kept == 0) ++zero_kept;
        if (s.kept_chains > 0 && s.segments == 0) ++kept_no_seg;
    }
    logs.push_back(os.str());
    if (zero_raw > 0) {
        std::ostringstream warn;
        warn << "z=" << z << " seeding: " << zero_raw << " NURBS face(s) with 0 raw seeds";
        logs.push_back(warn.str());
    }
    if (zero_kept > 0) {
        std::ostringstream warn;
        warn << "z=" << z << " seeding: " << zero_kept
             << " NURBS face(s) with raw seeds but 0 kept after trim filter";
        logs.push_back(warn.str());
    }
    if (kept_no_seg > 0) {
        std::ostringstream warn;
        warn << "z=" << z << " seeding: " << kept_no_seg
             << " NURBS face(s) with kept chains but 0 segments";
        logs.push_back(warn.str());
    }
}

}  // namespace

std::shared_ptr<IShape> readStep(const std::string& path) {
    return ShapeEngine::Kernel().readStep(path);
}

SliceResult sliceShape(const std::shared_ptr<IShape>& shape, const SliceOptions& opt) {
    SliceResult result;
    result.tolerance = opt.tolerance;
    result.stitch_tolerance = opt.stitch_tolerance;
    result.normal = opt.normal;
    const auto t_slice = SteadyClock::now();

    if (!shape || shape->empty()) {
        result.logs.push_back("empty shape");
        return result;
    }

    const SliceFrame frame = makeSliceFrame(opt.normal);

    auto t0 = SteadyClock::now();
    FaceIndex index;
    index.build(ShapeEngine::Kernel().exploreFaces(*shape), shape->bbox());
    const double ms_index = elapsedMs(t0);
    const int total_faces = static_cast<int>(index.faces().size());
    {
        std::ostringstream os;
        os << "shape faces=" << total_faces;
        result.logs.push_back(os.str());
    }

    double dmin = 0, dmax = 0;
    index.shapeRange(opt.normal, dmin, dmax);

    std::vector<double> heights = opt.explicit_heights;
    if (heights.empty()) {
        const double span = dmax - dmin;
        double lh = opt.layer_height;
        int count = opt.layer_count;
        if (count > 0) {
            if (span <= opt.tolerance) {
                throw std::runtime_error("shape has no extent along slice normal");
            }
            lh = span / static_cast<double>(count);
        } else {
            if (lh <= 0) throw std::runtime_error("layer-height must be positive");
        }
        double start = opt.start_height;
        if (!opt.start_height_set) {
            start = dmin + 0.5 * lh;
        }
        if (count < 0) {
            count = static_cast<int>(std::floor((dmax - start) / lh)) + 1;
        }
        heights.reserve(static_cast<size_t>(std::max(count, 0)));
        for (int i = 0; i < count; ++i) {
            const double h = start + i * lh;
            if (h > dmax + opt.tolerance) break;
            if (h < dmin - opt.tolerance) continue;
            heights.push_back(h);
        }
    }

    int total_bridge = 0;
    double max_bridge = 0;
    int total_open_contours = 0;
    std::vector<std::pair<double, int>> unclosed_layers;
    double ms_query = 0;
    double ms_intersect = 0;
    double ms_dedupe = 0;
    double ms_clip = 0;
    double ms_assemble = 0;
    for (double h : heights) {
        if (opt.throw_on_cancel)
            opt.throw_on_cancel();
        const Plane pln = makePlane(opt.normal, h);
        t0 = SteadyClock::now();
        auto candidates = index.query(opt.normal, h, opt.tolerance);
        ms_query += elapsedMs(t0);
        std::vector<RawSegment> segs;
        SeedLayer seed_layer;
        seed_layer.z = h;
        std::vector<FaceSeedStats> face_seed_stats;
        if (!opt.seed_dxf_path.empty()) opt.seed_out = &seed_layer;
        opt.face_seed_stats = &face_seed_stats;
        opt.plane_faces = &candidates;
        std::vector<std::string> constraint_audit;
        if (!opt.seed_dxf_path.empty()) opt.constraint_audit = &constraint_audit;
        t0 = SteadyClock::now();
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
        ms_intersect += elapsedMs(t0);
        std::set<int> active_faces;
        for (const FaceSeedStats& s : face_seed_stats) {
            if (s.segments > 0 || s.raw > 0 || s.kept > 0) active_faces.insert(s.face_id);
        }
        for (const RawSegment& s : segs) active_faces.insert(s.face_id);
        {
            std::ostringstream os;
            os << "z=" << h << " plane_faces=" << candidates.size()
               << " active_faces=" << active_faces.size();
            result.logs.push_back(os.str());
        }
        t0 = SteadyClock::now();
        dedupeAnalyticArcs(segs, opt.tolerance);
        ms_dedupe += elapsedMs(t0);
        opt.seed_out = nullptr;
        opt.face_seed_stats = nullptr;
        opt.plane_faces = nullptr;
        for (const std::string& line : constraint_audit) result.logs.push_back(line);
        opt.constraint_audit = nullptr;
        appendFaceSeedLogs(h, face_seed_stats, result.logs);
        t0 = SteadyClock::now();
        ShapeEngine::SectionRef().clipSegmentsToSection(*shape, pln, segs, opt.tolerance);
        ms_clip += elapsedMs(t0);
        AssembleStats stats;
        t0 = SteadyClock::now();
        Layer layer = assembleLayer(h, std::move(segs), frame, opt, stats);
        ms_assemble += elapsedMs(t0);
        if (!opt.seed_dxf_path.empty()) {
            seed_layer.face_stats = std::move(face_seed_stats);
            result.seed_layers.push_back(std::move(seed_layer));
        }
        if (stats.bridged > 0) {
            total_bridge += stats.bridged;
            max_bridge = std::max(max_bridge, stats.max_bridge);
            std::ostringstream os;
            os << "z=" << h << " bridged " << stats.bridged
               << " gap(s), max=" << stats.max_bridge;
            result.logs.push_back(os.str());
        }
        if (stats.open_leftover > 0) {
            total_open_contours += stats.open_leftover;
            unclosed_layers.push_back({h, stats.open_leftover});
            std::ostringstream os;
            os << "z=" << h << " not closed: " << stats.open_leftover << " open contour(s)";
            result.logs.push_back(os.str());
        }
        result.layers.push_back(std::move(layer));
    }

    if (total_bridge > 0) {
        std::ostringstream os;
        os << "stitched " << total_bridge << " open gap(s), max bridge distance " << max_bridge;
        result.logs.push_back(os.str());
    }

    if (!unclosed_layers.empty()) {
        std::ostringstream os;
        os << "unclosed layers (" << unclosed_layers.size() << "/" << heights.size()
           << ", " << total_open_contours << " open contour(s)): ";
        for (size_t i = 0; i < unclosed_layers.size(); ++i) {
            if (i > 0) os << ", ";
            os << "z=" << unclosed_layers[i].first << " (" << unclosed_layers[i].second << ")";
        }
        result.logs.push_back(os.str());
    }

    const int layer_count = static_cast<int>(heights.size());
    appendTimingLog(result.logs, "index", ms_index);
    appendTimingLog(result.logs, "query", ms_query, layer_count);
    appendTimingLog(result.logs, "intersect", ms_intersect, layer_count);
    appendTimingLog(result.logs, "dedupe", ms_dedupe, layer_count);
    appendTimingLog(result.logs, "section_clip", ms_clip, layer_count);
    appendTimingLog(result.logs, "assemble", ms_assemble, layer_count);

    double ms_export = 0;
    if (!opt.svg_dir.empty()) {
        t0 = SteadyClock::now();
        writeAllLayerSvg(result, opt.svg_dir);
        ms_export += elapsedMs(t0);
    }
    if (!opt.dxf_path.empty()) {
        t0 = SteadyClock::now();
        writeAllLayerDxf(result, opt.dxf_path);
        ms_export += elapsedMs(t0);
    }
    if (!opt.open_dxf_dir.empty()) {
        t0 = SteadyClock::now();
        writeOpenLayerDxf(result, unclosed_layers, opt.open_dxf_dir);
        ms_export += elapsedMs(t0);
    }
    if (!opt.seed_dxf_path.empty()) {
        t0 = SteadyClock::now();
        writeSeedDxfFile(result, opt.seed_dxf_path);
        std::string summary = opt.seed_dxf_path;
        if (summary.size() >= 4 && summary.compare(summary.size() - 4, 4, ".dxf") == 0)
            summary.replace(summary.size() - 4, 4, "_summary.txt");
        else
            summary += "_summary.txt";
        writeSeedSummaryFile(result, summary);
        ms_export += elapsedMs(t0);
    }
    if (ms_export > 0) appendTimingLog(result.logs, "export", ms_export);
    appendTimingLog(result.logs, "slice_total", elapsedMs(t_slice));
    return result;
}

SliceResult sliceFile(const std::string& stepPath, const SliceOptions& opt) {
    return sliceShape(readStep(stepPath), opt);
}

}  // namespace brepslicer
