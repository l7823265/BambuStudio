#pragma once

#include <cstdint>
#include <functional>
#include <optional>
#include <string>
#include <vector>

namespace brepslicer {

struct Vec3 {
    double x = 0;
    double y = 0;
    double z = 0;
};

enum class SegmentType {
    Line,
    Arc,
    Ellipse,
    BSpline
};

struct Segment {
    SegmentType type = SegmentType::Line;
    Vec3 start;
    Vec3 end;

    // Arc / ellipse (angles in radians, slice-plane frame: X = Ax2.XDirection of (origin, n)).
    Vec3 center;
    Vec3 normal;
    Vec3 major_axis;
    double radius = 0;   // circle radius, or ellipse semi-major a
    double radius_b = 0; // ellipse semi-minor b
    double start_angle = 0;
    double sweep = 0;

    // L2 placeholder
    int degree = 3;
    std::vector<double> knots;
    std::vector<Vec3> ctrl_pts;
    std::vector<double> weights;
    double fit_error = 0;
};

struct Contour {
    std::string orientation;          // "outer" | "inner"
    std::optional<int> parent;        // index in the same layer, or null
    bool closed = true;
    bool coplanar = false;
    int solid_id = -1;
    int shell_id = -1;
    std::vector<Segment> segments;
};

struct Layer {
    double z = 0;                     // plane offset d = n · p
    std::vector<Contour> contours;
};

struct TopologyEvent {
    double z = 0;
    std::string event;
    std::vector<int> from;
    std::vector<int> to;
    int contour_id = -1;
    Vec3 point;
    bool has_point = false;
};

struct SliceResult {
    std::string unit = "mm";
    double tolerance = 1e-4;
    Vec3 normal{0, 0, 1};
    std::vector<Layer> layers;
    std::vector<TopologyEvent> topology_events;
    std::vector<std::string> logs;
};

struct SliceOptions {
    Vec3 normal{0, 0, 1};
    double layer_height = 0.2;
    double start_height = 0;
    bool start_height_set = false;
    int layer_count = -1;
    std::vector<double> explicit_heights;
    double tolerance = 1e-4;       // endpoint stitch / JSON
    double geom_tolerance = 1e-9;  // analytic equality (L1 circle)
    double angular_tolerance = 1e-10;
    std::string svg_dir;           // directory: one SVG per layer
    std::string dxf_path;          // .dxf file, or directory of per-layer DXF
    // Optional: called once per layer so hosts can cancel a long B-rep slice.
    std::function<void()> throw_on_cancel;
};

struct VerifyReport {
    bool ok = true;
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
    std::vector<std::string> notes;
};

}  // namespace brepslicer
