#include <brepslicer/Slicer.h>

#include <contour/ContourAssembler.h>
#include <geom/GeomUtil.h>
#include <intersect/AnalyticIntersect.h>
#include <topo/FaceIndex.h>
#include <topo/SolidAdjacency.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <vector>

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
    SolidAdjacency solid_adj = ShapeEngine::Kernel().buildSolidAdjacency(index.faces());
    opt.solid_adjacency = solid_adj.empty() ? nullptr : &solid_adj;
    {
        int kind_count[6] = {};
        int cyl_vert = 0, cyl_horiz = 0, cyl_oblique = 0;
        int cone_vert = 0, cone_horiz = 0;
        int plane_horiz = 0, plane_vert = 0;
        int tor = 0;
        const Vec3 n = normalized(opt.normal);
        for (const FaceRecord& fr : index.faces()) {
            if (!fr.face) continue;
            const SurfaceData sd = fr.face->surface();
            const int k = static_cast<int>(sd.kind);
            if (k >= 0 && k < 6) ++kind_count[k];
            if (sd.kind == SurfaceKind::Cylinder) {
                const double nv = std::abs(dot(n, sd.cylinder.axis));
                if (nv >= 0.707) ++cyl_vert;
                else if (nv <= 0.1) ++cyl_horiz;
                else ++cyl_oblique;
            } else if (sd.kind == SurfaceKind::Cone) {
                const double nv = std::abs(dot(n, sd.cone.axis));
                if (nv >= 0.707) ++cone_vert;
                else ++cone_horiz;
            } else if (sd.kind == SurfaceKind::Torus) {
                ++tor;
            } else if (sd.kind == SurfaceKind::Plane) {
                const double pn = std::abs(dot(n, sd.plane.n));
                if (pn >= 0.707) ++plane_horiz;  // face // slice → coplanar candidate
                else ++plane_vert;
            }
        }
        std::ostringstream os;
        os << "face_kinds Plane=" << kind_count[0] << " Cyl=" << kind_count[1]
           << " Sph=" << kind_count[2] << " Cone=" << kind_count[3]
           << " Tor=" << kind_count[4] << " Other=" << kind_count[5]
           << " | cyl: vert=" << cyl_vert << " horiz=" << cyl_horiz
           << " oblique=" << cyl_oblique << " cone_vert=" << cone_vert
           << " cone_horiz=" << cone_horiz << " plane_horiz=" << plane_horiz
           << " plane_vert=" << plane_vert;
        result.logs.push_back(os.str());
    }
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

        // Solid-first: intersect + assemble each solid independently, then cross-solid nest.
        std::map<int, std::vector<const FaceRecord*>> by_solid;
        for (const FaceRecord* f : candidates) {
            if (f) by_solid[f->solid_id].push_back(f);
        }

        SeedLayer seed_layer;
        seed_layer.z = h;
        std::vector<FaceSeedStats> face_seed_stats;
        const bool want_seeds =
            !opt.seed_dxf_path.empty() || !opt.seed_open_dxf_dir.empty();
        if (want_seeds) opt.seed_out = &seed_layer;
        opt.face_seed_stats = &face_seed_stats;
        std::vector<std::string> constraint_audit;
        if (want_seeds) opt.constraint_audit = &constraint_audit;

        std::vector<Contour> all_contours;
        AssembleStats stats;
        double ms_dedupe_layer = 0;
        double ms_clip_layer = 0;
        double ms_assemble_layer = 0;
        for (auto& kv : by_solid) {
            std::vector<const FaceRecord*>& solid_faces = kv.second;
            opt.plane_faces = &solid_faces;

            std::vector<RawSegment> segs;
            int emit_by_surf[6] = {};
            int emit_line = 0, emit_arc = 0, emit_ell = 0, emit_bsp = 0, emit_empty = 0;
            t0 = SteadyClock::now();
            for (const FaceRecord* f : solid_faces) {
                auto part = intersectFaceWithPlane(*f, pln, frame, opt);
                const SurfaceData sd =
                    (f && f->face) ? f->face->surface() : SurfaceData{};
                const int sk = static_cast<int>(sd.kind);
                int n_keep = 0;
                for (RawSegment& s : part) {
                    if (s.degenerate) {
                        TopologyEvent ev;
                        ev.z = h;
                        ev.event = s.degen_event.empty() ? "tangent_point" : s.degen_event;
                        ev.point = s.degen_point;
                        ev.has_point = true;
                        result.topology_events.push_back(ev);
                    } else {
                        ++n_keep;
                        switch (s.geom.type) {
                        case SegmentType::Line: ++emit_line; break;
                        case SegmentType::Arc: ++emit_arc; break;
                        case SegmentType::Ellipse: ++emit_ell; break;
                        case SegmentType::BSpline: ++emit_bsp; break;
                        }
                        segs.push_back(std::move(s));
                    }
                }
                if (n_keep == 0)
                    ++emit_empty;
                else if (sk >= 0 && sk < 6)
                    emit_by_surf[sk] += n_keep;
            }
            {
                std::ostringstream os;
                os << "z=" << h << " emit_segs line=" << emit_line << " arc=" << emit_arc
                   << " ell=" << emit_ell << " bsp=" << emit_bsp
                   << " empty_faces=" << emit_empty << " by_surf[P,Cyl,Sph,Cone,Tor,Oth]="
                   << emit_by_surf[0] << "," << emit_by_surf[1] << "," << emit_by_surf[2]
                   << "," << emit_by_surf[3] << "," << emit_by_surf[4] << ","
                   << emit_by_surf[5];
                result.logs.push_back(os.str());
            }
            ms_intersect += elapsedMs(t0);

            t0 = SteadyClock::now();
            dedupeAnalyticArcs(segs, opt.tolerance);
            ms_dedupe_layer += elapsedMs(t0);

            t0 = SteadyClock::now();
            if (opt.use_occ_section_clip) {
                ShapeEngine::SectionRef().clipSegmentsToSection(*shape, pln, segs, opt.tolerance);
            }
            ms_clip_layer += elapsedMs(t0);

            AssembleStats solid_stats;
            t0 = SteadyClock::now();
            Layer solid_layer = assembleLayer(h, std::move(segs), frame, opt, solid_stats);
            ms_assemble_layer += elapsedMs(t0);

            stats.bridged += solid_stats.bridged;
            stats.max_bridge = std::max(stats.max_bridge, solid_stats.max_bridge);
            stats.spurious_closed_removed += solid_stats.spurious_closed_removed;
            for (Contour& c : solid_layer.contours) all_contours.push_back(std::move(c));
        }
        ms_dedupe += ms_dedupe_layer;
        ms_clip += ms_clip_layer;
        ms_assemble += ms_assemble_layer;

        opt.seed_out = nullptr;
        opt.face_seed_stats = nullptr;
        opt.plane_faces = nullptr;
        for (const std::string& line : constraint_audit) result.logs.push_back(line);
        opt.constraint_audit = nullptr;

        std::set<int> active_faces;
        for (const FaceSeedStats& s : face_seed_stats) {
            if (s.segments > 0 || s.raw > 0 || s.kept > 0) active_faces.insert(s.face_id);
        }
        {
            std::ostringstream os;
            os << "z=" << h << " plane_faces=" << candidates.size()
               << " solids=" << by_solid.size() << " active_faces=" << active_faces.size();
            result.logs.push_back(os.str());
        }
        appendFaceSeedLogs(h, face_seed_stats, result.logs);

        t0 = SteadyClock::now();
        reassembleCrossSolidContours(all_contours, frame, opt, stats);
        ms_assemble += elapsedMs(t0);

        Layer layer;
        layer.z = h;
        layer.contours = std::move(all_contours);

        if (opt.compare_occ_section) {
            try {
                const double refLen = ShapeEngine::SectionRef().sectionLength(*shape, pln);
                double ourLen = 0;
                for (const Contour& c : layer.contours) {
                    for (const Segment& s : c.segments) {
                        if (s.type == SegmentType::Line) {
                            ourLen += dist(s.start, s.end);
                        } else if (s.type == SegmentType::Arc) {
                            ourLen += std::abs(s.sweep) * s.radius;
                        } else if (s.type == SegmentType::BSpline && s.ctrl_pts.size() >= 2) {
                            for (size_t i = 0; i + 1 < s.ctrl_pts.size(); ++i)
                                ourLen += dist(s.ctrl_pts[i], s.ctrl_pts[i + 1]);
                        }
                    }
                }
                std::ostringstream os;
                os << std::fixed << std::setprecision(4);
                os << "z=" << h << " perimeter our=" << ourLen << " occt=" << refLen
                   << " rel=" << (std::abs(ourLen - refLen) / std::max(1e-9, 0.5 * (ourLen + refLen)));
                result.logs.push_back(os.str());
            } catch (...) {
                result.logs.push_back("z=" + std::to_string(h) + " occt section failed");
            }
        }

        if (!opt.seed_dxf_path.empty() || !opt.seed_open_dxf_dir.empty()) {
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
            summary += "/seeds_summary.txt";
        writeSeedSummaryFile(result, summary);
        ms_export += elapsedMs(t0);
    }
    if (!opt.seed_open_dxf_dir.empty()) {
        t0 = SteadyClock::now();
        writeOpenLayerSeedDxf(result, unclosed_layers, opt.seed_open_dxf_dir);
        writeSeedSummaryFile(result, opt.seed_open_dxf_dir + "/seeds_summary.txt");
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
