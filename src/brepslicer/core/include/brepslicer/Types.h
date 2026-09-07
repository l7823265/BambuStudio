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

enum class SeedPointKind { Raw, Kept, Dropped };

struct SeedPoint {
    Vec3 p;
    SeedPointKind kind = SeedPointKind::Raw;
    int face_id = -1;
};

struct FaceSeedStats {
    int face_id = -1;
    int raw = 0;          // iso hits after seeding
    int iso_hits = 0;     // iso curves with >= 1 hit
    int chains = 0;       // chains after chainHits
    int chain_pts = 0;    // points in chains before filter
    int kept = 0;         // points in chains kept after trim filter
    int kept_chains = 0;  // chains kept after trim filter
    int dropped = 0;      // points dropped by trim filter
    int segments = 0;     // fitted segments output
    int neighbor_splits = 0;  // chain clips by neighbor-face signed distance
    bool uvmarch_fallback = false;
};

struct SeedLayer {
    double z = 0;
    std::vector<SeedPoint> points;
    std::vector<std::vector<Vec3>> seed_chains;  // after seeding, before trim filter
    std::vector<std::vector<Vec3>> kept_chains;  // after trim filter
    std::vector<int> kept_chain_face_ids;        // parallel to kept_chains
    std::vector<FaceSeedStats> face_stats;
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

struct FaceRecord;

struct SliceResult {
    std::string unit = "mm";
    double tolerance = 1e-4;
    double stitch_tolerance = 0.01;
    Vec3 normal{0, 0, 1};
    std::vector<Layer> layers;
    std::vector<SeedLayer> seed_layers;
    std::vector<TopologyEvent> topology_events;
    std::vector<std::string> logs;
};

enum class NurbsMethod {
    Auto,     // UVMatch first; fall back to UVMarch if empty
    UvMatch,  // iso-sample / reverse-fit only (no march fallback)
    UvMarch   // UV marching only (skip iso-sample)
};

struct SliceOptions {
    Vec3 normal{0, 0, 1};
    double layer_height = 0.2;
    double start_height = 0;
    bool start_height_set = false;
    int layer_count = -1;
    std::vector<double> explicit_heights;
    double tolerance = 1e-4;       // endpoint stitch / JSON
    double stitch_tolerance = 0.01; // merge segment endpoints when gap < this (mm)
    double geom_tolerance = 1e-9;  // analytic equality (L1 circle)
    double angular_tolerance = 1e-10;
    NurbsMethod nurbs_method = NurbsMethod::Auto;
    std::string svg_dir;           // directory: one SVG per layer
    std::string dxf_path;          // .dxf file, or directory of per-layer DXF
    std::string open_dxf_dir;      // debug: DXF for layers with open contour(s) only
    std::string seed_dxf_path;     // debug: UVMatch seed points / chains
    // Per-layer UVMatch seed recorder (set by slicer during intersect).
    mutable SeedLayer* seed_out = nullptr;
    // Per-layer per-face seed counts (always set during slice).
    mutable std::vector<FaceSeedStats>* face_seed_stats = nullptr;
    // All faces hit by the current slice plane (for neighbor signed-distance clip).
    mutable const std::vector<const FaceRecord*>* plane_faces = nullptr;
    // Per-face chain/constraint audit lines (debug; set when --seed-dxf is used).
    mutable std::vector<std::string>* constraint_audit = nullptr;
    // Optional: called once per layer so hosts (BambuStudio) can cancel a long B-rep slice.
    std::function<void()> throw_on_cancel;
};

struct VerifyReport {
    bool ok = true;
    std::vector<std::string> errors;
    std::vector<std::string> warnings;
    std::vector<std::string> notes;
};

}  // namespace brepslicer
