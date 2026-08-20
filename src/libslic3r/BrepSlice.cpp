#include "BrepSlice.hpp"

#include "ClipperUtils.hpp"
#include "Model.hpp"
#include "Polygon.hpp"

#include <brepslicer/Geom.h>
#include <brepslicer/Slicer.h>
#include <intersect/BSplineFit.h>

#include <boost/algorithm/string/predicate.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/log/trivial.hpp>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <fstream>
#include <stdexcept>
#include <sstream>

namespace Slic3r {

// #region agent log
static void agent_dbg(const char *hypothesisId, const char *location, const char *message, const std::string &data_json)
{
    try {
        std::ofstream f("E:/learning/slicer/BambuStudio/debug-9ef780.log", std::ios::app);
        if (!f)
            return;
        const auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                            std::chrono::system_clock::now().time_since_epoch())
                            .count();
        f << "{\"sessionId\":\"9ef780\",\"runId\":\"brep-gcode\",\"hypothesisId\":\"" << hypothesisId
          << "\",\"location\":\"" << location << "\",\"message\":\"" << message << "\",\"data\":" << data_json
          << ",\"timestamp\":" << ms << "}\n";
    } catch (...) {}
}
// #endregion

static bool is_step_path(const std::string &path)
{
    return boost::algorithm::iends_with(path, ".step") || boost::algorithm::iends_with(path, ".stp");
}

static std::string step_path_of(const ModelObject &object)
{
    if (is_step_path(object.input_file))
        return object.input_file;
    for (const ModelVolume *v : object.volumes) {
        if (v && is_step_path(v->source.input_file))
            return v->source.input_file;
    }
    return {};
}

bool model_object_has_brep_step(const ModelObject &object)
{
    const std::string path = step_path_of(object);
    if (path.empty())
        return false;
    boost::system::error_code ec;
    return boost::filesystem::exists(boost::filesystem::path(path), ec);
}

static Point to_scaled_xy(const brepslicer::Vec3 &p, const Transform3d &T)
{
    const Vec3d q = T * Vec3d(p.x, p.y, p.z);
    return Point::new_scale(q.x(), q.y());
}

static int sample_count(double radius, double sweep, double chord)
{
    radius = std::max(std::abs(radius), 1e-9);
    chord  = std::max(chord, 1e-6);
    double x = 1.0 - chord / radius;
    x        = std::clamp(x, -1.0, 1.0);
    double dtheta = 2.0 * std::acos(x);
    dtheta        = std::clamp(dtheta, 1.0 * brepslicer::kPi / 180.0, 30.0 * brepslicer::kPi / 180.0);
    return std::max(2, int(std::ceil(std::abs(sweep) / dtheta)));
}

static void append_segment_open(Points                           &pts,
                                const brepslicer::Segment        &s,
                                const Transform3d                &T,
                                double                            chord)
{
    auto add = [&](const brepslicer::Vec3 &p) {
        Point q = to_scaled_xy(p, T);
        if (pts.empty() || pts.back() != q)
            pts.push_back(q);
    };

    if (s.type == brepslicer::SegmentType::Line) {
        add(s.start);
        return;
    }

    if (s.type == brepslicer::SegmentType::BSpline) {
        const double len = std::max(brepslicer::bsplineLength(s), chord);
        const int    n   = std::max(2, int(std::ceil(len / chord)) + 1);
        std::vector<brepslicer::Vec3> samples;
        brepslicer::bsplineSample(s, n, samples);
        for (size_t i = 0; i + 1 < samples.size(); ++i)
            add(samples[i]);
        return;
    }

    if (s.type == brepslicer::SegmentType::Arc) {
        const int n = sample_count(s.radius, s.sweep, chord);
        const brepslicer::SliceFrame frame = brepslicer::makeSliceFrame(s.normal);
        for (int i = 0; i < n; ++i) {
            const double a = s.start_angle + s.sweep * (double(i) / double(n));
            const brepslicer::Vec3 p = s.center
                + frame.x * (s.radius * std::cos(a))
                + frame.y * (s.radius * std::sin(a));
            add(p);
        }
        return;
    }

    if (s.type == brepslicer::SegmentType::Ellipse) {
        const double a = s.radius;
        const double b = s.radius_b > 0 ? s.radius_b : s.radius;
        const int    n = sample_count(std::max(a, b), s.sweep, chord);
        brepslicer::Vec3 major = s.major_axis;
        const double ml        = brepslicer::length(major);
        if (ml > 1e-15)
            major = major / ml;
        else
            major = brepslicer::makeSliceFrame(s.normal).x;
        const brepslicer::Vec3 minor = brepslicer::cross(s.normal, major);
        for (int i = 0; i < n; ++i) {
            const double t = s.start_angle + s.sweep * (double(i) / double(n));
            const brepslicer::Vec3 p = s.center + major * (a * std::cos(t)) + minor * (b * std::sin(t));
            add(p);
        }
        return;
    }

    add(s.start);
}

static Polygon contour_to_polygon(const brepslicer::Contour &c, const Transform3d &T, double chord)
{
    Points pts;
    pts.reserve(c.segments.size() * 4);
    for (const brepslicer::Segment &s : c.segments)
        append_segment_open(pts, s, T, chord);
    if (!c.segments.empty() && !c.closed)
        pts.push_back(to_scaled_xy(c.segments.back().end, T));
    if (pts.size() >= 3 && pts.front() == pts.back())
        pts.pop_back();
    return Polygon(std::move(pts));
}

// Near the solid tip / grazing planes, ContourAssembler may leave a tiny gap and
// mark the loop open. Treat those as closed so we do not discard valid area.
static bool contour_usable_as_polygon(const brepslicer::Contour &c, double gap_tol)
{
    if (c.segments.empty())
        return false;
    if (c.closed)
        return true;
    const brepslicer::Vec3 &a = c.segments.front().start;
    const brepslicer::Vec3 &b = c.segments.back().end;
    const double dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return (dx * dx + dy * dy + dz * dz) <= gap_tol * gap_tol;
}

static ExPolygons layer_to_expolygons(const brepslicer::Layer   &layer,
                                      const Transform3d         &T,
                                      double                     chord,
                                      int                        solid_id_filter)
{
    Polygons polys;
    polys.reserve(layer.contours.size());
    const double gap_tol = std::max(chord * 2.0, 1e-3);
    for (const brepslicer::Contour &c : layer.contours) {
        if (solid_id_filter >= 0 && c.solid_id != solid_id_filter)
            continue;
        if (!contour_usable_as_polygon(c, gap_tol))
            continue;
        Polygon p = contour_to_polygon(c, T, chord);
        if (p.size() >= 3)
            polys.emplace_back(std::move(p));
    }
    if (polys.empty())
        return {};
    if (T.linear().determinant() < 0.) {
        for (Polygon &p : polys)
            p.reverse();
    }
    return union_ex(polys);
}

// Match mesh slicing: trafo_centered * volume_matrix * unit_scale * (p_step - mesh_offset).
// mesh_offset is recorded before convert_from_inches/meters, so it stays in STEP units;
// the mesh itself was scaled to mm. Apply the same unit scale here.
static double volume_unit_scale_to_mm(const ModelVolume &volume)
{
    if (volume.source.is_converted_from_meters)
        return 1000.0;
    if (volume.source.is_converted_from_inches)
        return 25.4;
    return 1.0;
}

static Transform3d brep_volume_trafo(const Transform3d &trafo_centered, const ModelVolume &volume)
{
    const double us = volume_unit_scale_to_mm(volume);
    return trafo_centered * volume.get_matrix() * Eigen::Scaling(us, us, us) *
           Eigen::Translation3d(-volume.source.mesh_offset);
}

static bool layer_has_geometry(const std::vector<ExPolygons> &layers)
{
    for (const ExPolygons &ex : layers)
        if (!ex.empty())
            return true;
    return false;
}

static bool bottom_layers_empty(const std::vector<ExPolygons> &layers, size_t count)
{
    count = std::min(count, layers.size());
    for (size_t i = 0; i < count; ++i)
        if (!layers[i].empty())
            return false;
    return count > 0;
}

bool slice_model_object_brep(const ModelObject           &object,
                             const Transform3d           &trafo_centered,
                             const std::vector<float>    &zs,
                             double                       chord_error,
                             const std::function<void()> &throw_on_cancel,
                             std::vector<BrepVolumeSlices> &out_by_volume)
{
    out_by_volume.clear();
    if (zs.empty()) {
        // #region agent log
        agent_dbg("C", "BrepSlice.cpp:entry", "reject_empty_zs", "{\"reason\":\"zs_empty\"}");
        // #endregion
        return false;
    }

    const std::string path = step_path_of(object);
    if (path.empty()) {
        // #region agent log
        agent_dbg("C", "BrepSlice.cpp:entry", "reject_no_step_path", "{\"reason\":\"path_empty\"}");
        // #endregion
        return false;
    }

    for (const ModelVolume *v : object.volumes) {
        if (!v)
            continue;
        if (v->is_negative_volume() || v->is_modifier() || v->is_mm_painted() || v->is_fuzzy_skin_facets_painted()) {
            // #region agent log
            agent_dbg("C", "BrepSlice.cpp:gates", "reject_modifier_or_paint",
                      std::string("{\"volume\":\"") + v->name + "\"}");
            // #endregion
            return false;
        }
    }

    std::vector<const ModelVolume *> parts;
    parts.reserve(object.volumes.size());
    for (const ModelVolume *v : object.volumes) {
        if (v && v->is_model_part())
            parts.push_back(v);
    }
    if (parts.empty()) {
        // #region agent log
        agent_dbg("C", "BrepSlice.cpp:gates", "reject_no_model_parts", "{}");
        // #endregion
        return false;
    }

    const bool filter_by_solid = parts.size() > 1;
    double unit_scale = 1.0;
    for (const ModelVolume *p : parts)
        unit_scale = std::max(unit_scale, volume_unit_scale_to_mm(*p));

    std::shared_ptr<brepslicer::IShape> shape;
    try {
        shape = brepslicer::readStep(path);
    } catch (const std::exception &e) {
        BOOST_LOG_TRIVIAL(warning) << "BrepSlicer failed to read STEP, falling back to mesh: " << e.what();
        // #region agent log
        agent_dbg("C", "BrepSlice.cpp:readStep", "read_exception",
                  std::string("{\"what\":\"") + e.what() + "\"}");
        // #endregion
        return false;
    }
    if (!shape || shape->empty()) {
        // #region agent log
        agent_dbg("C", "BrepSlice.cpp:readStep", "empty_shape", "{}");
        // #endregion
        return false;
    }

    const brepslicer::BBox bb = shape->bbox();
    const double step_dx = (bb.xmax - bb.xmin) * unit_scale;
    const double step_dy = (bb.ymax - bb.ymin) * unit_scale;
    const double step_dz = (bb.zmax - bb.zmin) * unit_scale;
    const double step_diag = std::sqrt(step_dx * step_dx + step_dy * step_dy + step_dz * step_dz);
    const double step_h   = std::abs(step_dz);
    const double print_h  = zs.size() < 2 ? 0.0 : std::abs(double(zs.back()) - double(zs.front()));

    // Architectural / meter-scale STEP often arrives already scaled down on the mesh
    // for the bed. B-rep still cuts the native CAD once per layer and can hang for
    // minutes/hours (e.g. ~2.5m pier → hundreds of full solid sections). Prefer mesh.
    double min_vol_scale = 1.0;
    for (const ModelVolume *p : parts) {
        const Vec3d s = p->get_scaling_factor();
        min_vol_scale = std::min(min_vol_scale,
                                 std::min({std::abs(s.x()), std::abs(s.y()), std::abs(s.z())}));
    }

    // #region agent log
    {
        std::ostringstream d;
        d << "{\"parts\":" << parts.size()
          << ",\"zs\":" << zs.size() << ",\"z0\":" << (zs.empty() ? 0.0 : zs.front())
          << ",\"zLast\":" << (zs.empty() ? 0.0 : zs.back())
          << ",\"bbox_zmin\":" << bb.zmin << ",\"bbox_zmax\":" << bb.zmax
          << ",\"unit_scale\":" << unit_scale
          << ",\"step_diag_mm\":" << step_diag << ",\"step_h_mm\":" << step_h
          << ",\"print_h\":" << print_h << ",\"min_vol_scale\":" << min_vol_scale
          << ",\"filter_by_solid\":" << (filter_by_solid ? "true" : "false") << "}";
        agent_dbg("H1", "BrepSlice.cpp:start", "shape_and_zs", d.str());
    }
    // #endregion

    const bool oversized_step = step_diag > 400.0;                 // > typical print bed
    const bool heavily_scaled = min_vol_scale < 0.5;              // mesh shrunk vs CAD
    const bool print_vs_cad   = print_h > 1e-6 && step_h > 5.0 * print_h;
    if (oversized_step && (heavily_scaled || print_vs_cad)) {
        BOOST_LOG_TRIVIAL(warning)
            << "BrepSlicer: STEP too large for interactive B-rep (diag=" << step_diag
            << "mm, step_h=" << step_h << "mm, print_h=" << print_h
            << "mm, scale=" << min_vol_scale << "), falling back to scaled mesh";
        // #region agent log
        {
            std::ostringstream d;
            d << "{\"step_diag_mm\":" << step_diag << ",\"step_h_mm\":" << step_h
              << ",\"print_h\":" << print_h << ",\"min_vol_scale\":" << min_vol_scale
              << ",\"unit_scale\":" << unit_scale
              << ",\"oversized\":true,\"heavily_scaled\":" << (heavily_scaled ? "true" : "false")
              << ",\"print_vs_cad\":" << (print_vs_cad ? "true" : "false") << "}";
            agent_dbg("H1", "BrepSlice.cpp:gates", "fallback_oversized_step", d.str());
        }
        // #endregion
        return false;
    }

    out_by_volume.reserve(parts.size());
    bool any = false;

    for (size_t part_idx = 0; part_idx < parts.size(); ++part_idx) {
        const ModelVolume *part = parts[part_idx];
        const double       part_us = volume_unit_scale_to_mm(*part);
        const double       chord_step = chord_error / std::max(part_us, 1e-12);
        const Transform3d  T    = brep_volume_trafo(trafo_centered, *part);
        const Eigen::Matrix3d R = T.linear();
        const Vec3d           t = T.translation();
        Vec3d                 n_step = R.transpose() * Vec3d(0, 0, 1);
        const double          nlen   = n_step.norm();
        if (nlen < 1e-12)
            return false;
        n_step /= nlen;

        brepslicer::SliceOptions opt;
        opt.normal = {n_step.x(), n_step.y(), n_step.z()};
        opt.explicit_heights.reserve(zs.size());
        for (float z : zs)
            opt.explicit_heights.push_back((double(z) - t.z()) / nlen);

        // #region agent log
        {
            std::ostringstream d;
            d << "{\"volume\":\"" << part->name << "\",\"unit_scale\":" << part_us
              << ",\"from_inches\":" << (part->source.is_converted_from_inches ? "true" : "false")
              << ",\"from_meters\":" << (part->source.is_converted_from_meters ? "true" : "false")
              << ",\"mesh_offset\":["
              << part->source.mesh_offset.x() << "," << part->source.mesh_offset.y() << ","
              << part->source.mesh_offset.z() << "],\"t\":[" << t.x() << "," << t.y() << "," << t.z()
              << "],\"nlen\":" << nlen
              << ",\"h0\":" << (opt.explicit_heights.empty() ? 0.0 : opt.explicit_heights.front())
              << ",\"hLast\":" << (opt.explicit_heights.empty() ? 0.0 : opt.explicit_heights.back())
              << ",\"n\":[" << n_step.x() << "," << n_step.y() << "," << n_step.z() << "]}";
            agent_dbg("H2", "BrepSlice.cpp:heights", "mapped_heights", d.str());
        }
        // #endregion

        if (throw_on_cancel)
            throw_on_cancel();

        // #region agent log
        agent_dbg("H1", "BrepSlice.cpp:sliceShape", "sliceShape_begin",
                  std::string("{\"layers\":") + std::to_string(zs.size()) + "}");
        // #endregion

        brepslicer::SliceResult result;
        try {
            BOOST_LOG_TRIVIAL(info) << "BrepSlicer: slicing STEP " << path << " volume=" << part->name
                                    << " layers=" << zs.size()
                                    << " unit_scale=" << part_us
                                    << " mesh_offset=(" << part->source.mesh_offset.x() << ","
                                    << part->source.mesh_offset.y() << "," << part->source.mesh_offset.z() << ")";
            result = brepslicer::sliceShape(shape, opt);
        } catch (const std::exception &e) {
            BOOST_LOG_TRIVIAL(warning) << "BrepSlicer failed, falling back to mesh: " << e.what();
            // #region agent log
            agent_dbg("C", "BrepSlice.cpp:sliceShape", "slice_exception",
                      std::string("{\"what\":\"") + e.what() + "\"}");
            // #endregion
            return false;
        }

        // #region agent log
        agent_dbg("H1", "BrepSlice.cpp:sliceShape", "sliceShape_end",
                  std::string("{\"result_layers\":") + std::to_string(result.layers.size()) + "}");
        // #endregion

        if (result.layers.size() != zs.size()) {
            BOOST_LOG_TRIVIAL(warning) << "BrepSlicer layer count mismatch (" << result.layers.size()
                                       << " vs " << zs.size() << "), falling back to mesh";
            // #region agent log
            {
                std::ostringstream d;
                d << "{\"got\":" << result.layers.size() << ",\"want\":" << zs.size() << "}";
                agent_dbg("C", "BrepSlice.cpp:sliceShape", "layer_count_mismatch", d.str());
            }
            // #endregion
            return false;
        }

        const int solid_id = filter_by_solid ? int(part_idx) : -1;

        // #region agent log
        {
            size_t raw0 = result.layers.empty() ? 0 : result.layers[0].contours.size();
            size_t raw1 = result.layers.size() > 1 ? result.layers[1].contours.size() : 0;
            size_t closed0 = 0, open0 = 0, closed1 = 0, open1 = 0;
            auto tally = [](const brepslicer::Layer &L, size_t &cl, size_t &op) {
                for (const brepslicer::Contour &c : L.contours) {
                    if (c.closed) ++cl;
                    else ++op;
                }
            };
            if (!result.layers.empty())
                tally(result.layers[0], closed0, open0);
            if (result.layers.size() > 1)
                tally(result.layers[1], closed1, open1);
            size_t nonempty_raw = 0;
            for (size_t i = 0; i < std::min<size_t>(result.layers.size(), 8); ++i)
                if (!result.layers[i].contours.empty())
                    ++nonempty_raw;
            std::ostringstream d;
            d << "{\"solid_id_filter\":" << solid_id << ",\"raw_contours_l0\":" << raw0
              << ",\"raw_contours_l1\":" << raw1 << ",\"closed0\":" << closed0 << ",\"open0\":" << open0
              << ",\"closed1\":" << closed1 << ",\"open1\":" << open1
              << ",\"nonempty_raw_first8\":" << nonempty_raw << "}";
            agent_dbg("D", "BrepSlice.cpp:raw", "raw_layer_contours", d.str());
        }
        // #endregion

        BrepVolumeSlices vs;
        vs.volume_id = part->id();
        vs.layers.resize(zs.size());
        for (size_t i = 0; i < zs.size(); ++i) {
            if (throw_on_cancel)
                throw_on_cancel();
            vs.layers[i] = layer_to_expolygons(result.layers[i], T, chord_step, solid_id);
            if (!vs.layers[i].empty())
                any = true;
        }

        // #region agent log
        {
            size_t poly0 = vs.layers.empty() ? 0 : vs.layers[0].size();
            size_t poly1 = vs.layers.size() > 1 ? vs.layers[1].size() : 0;
            size_t poly2 = vs.layers.size() > 2 ? vs.layers[2].size() : 0;
            size_t nonempty = 0;
            size_t first_nonempty = zs.size();
            double area0 = 0;
            for (size_t i = 0; i < vs.layers.size(); ++i) {
                if (!vs.layers[i].empty()) {
                    ++nonempty;
                    if (first_nonempty == zs.size())
                        first_nonempty = i;
                }
            }
            if (!vs.layers.empty()) {
                for (const ExPolygon &ex : vs.layers[0])
                    area0 += unscaled(unscaled(ex.area()));
            }
            std::ostringstream d;
            d << "{\"expoly_l0\":" << poly0 << ",\"expoly_l1\":" << poly1 << ",\"expoly_l2\":" << poly2
              << ",\"nonempty_layers\":" << nonempty << ",\"first_nonempty\":" << first_nonempty
              << ",\"l0_area_mm2\":" << area0
              << ",\"bottom3_empty\":" << (bottom_layers_empty(vs.layers, 3) ? "true" : "false") << "}";
            agent_dbg("B", "BrepSlice.cpp:convert", "converted_expolygons", d.str());
        }
        // #endregion

        if (!layer_has_geometry(vs.layers)) {
            BOOST_LOG_TRIVIAL(warning) << "BrepSlicer produced no contours for volume " << part->name
                                       << ", falling back to mesh";
            // #region agent log
            agent_dbg("B", "BrepSlice.cpp:convert", "fallback_no_geometry", "{}");
            // #endregion
            return false;
        }

        // Empty bottom layers cause "empty initial layer" G-code errors.
        // Require layer 0 specifically (not merely any of the first 3).
        if ((vs.layers.empty() || vs.layers[0].empty()) && layer_has_geometry(vs.layers)) {
            BOOST_LOG_TRIVIAL(warning) << "BrepSlicer bottom layers empty for volume " << part->name
                                       << ", falling back to mesh";
            // #region agent log
            agent_dbg("A", "BrepSlice.cpp:convert", "fallback_bottom_empty",
                      std::string("{\"first_nonempty_fix\":\"l0_empty\"}"));
            // #endregion
            return false;
        }

        for (const std::string &log : result.logs)
            BOOST_LOG_TRIVIAL(info) << "BrepSlicer: " << log;

        out_by_volume.emplace_back(std::move(vs));
    }

    if (!any) {
        BOOST_LOG_TRIVIAL(warning) << "BrepSlicer produced no contours, falling back to mesh";
        out_by_volume.clear();
        // #region agent log
        agent_dbg("B", "BrepSlice.cpp:exit", "fallback_any_false", "{}");
        // #endregion
        return false;
    }

    // #region agent log
    agent_dbg("E", "BrepSlice.cpp:exit", "brep_success",
              std::string("{\"volumes\":") + std::to_string(out_by_volume.size()) + "}");
    // #endregion
    return true;
}

} // namespace Slic3r
