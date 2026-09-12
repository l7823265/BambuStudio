#include <brepslicer/Slicer.h>
#include <geom/GeomUtil.h>
#include <intersect/BSplineFit.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#include <direct.h>
#else
#include <sys/stat.h>
#endif

namespace brepslicer {
namespace {

void ensureDir(const std::string& dir) {
#ifdef _WIN32
    _mkdir(dir.c_str());
#else
    mkdir(dir.c_str(), 0755);
#endif
}

bool endsWithDxf(const std::string& path) {
    if (path.size() < 4) return false;
    std::string ext = path.substr(path.size() - 4);
    for (char& c : ext) c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
    return ext == ".dxf";
}

std::string layerName(size_t i, double z) {
    std::ostringstream os;
    os << "SLICE_" << std::setw(4) << std::setfill('0') << i;
    (void)z;
    return os.str();
}

void code(std::ostream& os, int c, const std::string& v) { os << c << "\n" << v << "\n"; }
void code(std::ostream& os, int c, int v) { os << c << "\n" << v << "\n"; }
void code(std::ostream& os, int c, double v) {
    os << c << "\n" << std::setprecision(17) << v << "\n";
}

double rad2deg(double r) { return r * 180.0 / kPi; }

double wrapDeg(double d) {
    d = std::fmod(d, 360.0);
    if (d < 0) d += 360.0;
    return d;
}

Vec3 segmentStart3d(const Segment& s) {
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) return s.ctrl_pts.front();
    return s.start;
}

Vec3 segmentEnd3d(const Segment& s) {
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) return s.ctrl_pts.back();
    return s.end;
}

int contourEntityColor(const Contour& c) {
    if (!c.closed) return 6;  // magenta — open contour
    return (c.orientation == "inner") ? 1 : 5;
}

std::string contourEntityLayer(const std::string& base, const Contour& c) {
    if (!c.closed) return base + "_OPEN";
    return base;
}

void writeOpenGapMarker(std::ostream& os, const Contour& c, const std::string& lay, double z,
                        const SliceFrame& frame) {
    if (c.closed || c.segments.empty()) return;
    const Vec3 gap_a = segmentEnd3d(c.segments.back());
    const Vec3 gap_b = segmentStart3d(c.segments.front());
    const Vec2 a = frame.toXY(gap_a);
    const Vec2 b = frame.toXY(gap_b);
    code(os, 0, "LINE");
    code(os, 8, lay + "_OPEN_GAP");
    code(os, 62, 2);  // yellow — shows where the contour fails to close
    code(os, 10, a.x);
    code(os, 20, a.y);
    code(os, 30, z);
    code(os, 11, b.x);
    code(os, 21, b.y);
    code(os, 31, z);
    code(os, 0, "POINT");
    code(os, 8, lay + "_OPEN_GAP");
    code(os, 62, 2);
    code(os, 10, a.x);
    code(os, 20, a.y);
    code(os, 30, z);
    code(os, 0, "POINT");
    code(os, 8, lay + "_OPEN_GAP");
    code(os, 62, 2);
    code(os, 10, b.x);
    code(os, 20, b.y);
    code(os, 30, z);
}

void writeEntities(std::ostream& os, const SliceResult& result) {
    const SliceFrame frame = makeSliceFrame(result.normal);
    auto xy = [&](const Vec3& p) { return frame.toXY(p); };

    for (size_t li = 0; li < result.layers.size(); ++li) {
        const Layer& layer = result.layers[li];
        const std::string lay = layerName(li, layer.z);
        const double z = layer.z;

        for (const Contour& c : layer.contours) {
            const int color = contourEntityColor(c);
            const std::string clayer = contourEntityLayer(lay, c);
            for (const Segment& s : c.segments) {
                const Vec2 a = xy(s.start);
                const Vec2 b = xy(s.end);

                if (s.type == SegmentType::Arc) {
                    const Vec2 ctr = xy(s.center);
                    if (std::abs(std::abs(s.sweep) - kTwoPi) <= 1e-8) {
                        code(os, 0, "CIRCLE");
                        code(os, 8, clayer);
                        code(os, 62, color);
                        code(os, 10, ctr.x);
                        code(os, 20, ctr.y);
                        code(os, 30, z);
                        code(os, 40, s.radius);
                        continue;
                    }
                    // DXF ARC is always CCW from start to end.
                    double a0 = s.start_angle;
                    double a1 = s.start_angle + s.sweep;
                    if (s.sweep < 0) std::swap(a0, a1);
                    code(os, 0, "ARC");
                    code(os, 8, clayer);
                    code(os, 62, color);
                    code(os, 10, ctr.x);
                    code(os, 20, ctr.y);
                    code(os, 30, z);
                    code(os, 40, s.radius);
                    code(os, 50, wrapDeg(rad2deg(a0)));
                    code(os, 51, wrapDeg(rad2deg(a1)));
                    continue;
                }

                if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) {
                    code(os, 0, "SPLINE");
                    code(os, 8, clayer);
                    code(os, 62, color);
                    code(os, 70, c.closed ? 9 : 8);  // planar + optional closed
                    code(os, 71, s.degree);
                    code(os, 72, static_cast<int>(s.knots.size()));
                    code(os, 73, static_cast<int>(s.ctrl_pts.size()));
                    code(os, 74, 0);
                    for (double k : s.knots) code(os, 40, k);
                    for (const Vec3& p : s.ctrl_pts) {
                        const Vec2 q = xy(p);
                        code(os, 10, q.x);
                        code(os, 20, q.y);
                        code(os, 30, z);
                    }
                    if (!s.weights.empty()) {
                        for (double w : s.weights) code(os, 41, w);
                    }
                    continue;
                }

                if (s.type == SegmentType::Ellipse && s.radius > 0) {
                    // Tessellate elliptical arcs as 3D LINE chains (explicit Z), not
                    // ELLIPSE / LWPOLYLINE:
                    // - FreeCAD C++ importer ignores ELLIPSE start/end → full oval
                    // - FreeCAD often ignores LWPOLYLINE elevation (38) → arcs on z=0
                    // - HyperMesh drops ELLIPSE; LINE XYZ works in FreeCAD / Midas / HM
                    Vec3 maj = s.major_axis;
                    const double majLen = length(maj);
                    if (majLen > 1e-30) maj = maj * (1.0 / majLen);
                    Vec3 minv = cross(s.normal, maj);
                    const double minLen = length(minv);
                    if (minLen > 1e-30) minv = minv * (1.0 / minLen);
                    const bool full = std::abs(std::abs(s.sweep) - kTwoPi) <= 1e-8;
                    const int nSamp =
                        full ? 64
                             : std::max(8, static_cast<int>(std::ceil(std::abs(s.sweep) * 48.0 /
                                                                      kPi)));
                    auto ellPt = [&](double t) -> Vec2 {
                        const Vec3 p = s.center + maj * (s.radius * std::cos(t)) +
                                       minv * (s.radius_b * std::sin(t));
                        return xy(p);
                    };
                    for (int i = 0; i < nSamp; ++i) {
                        const double t0 =
                            s.start_angle +
                            s.sweep * (static_cast<double>(i) / static_cast<double>(nSamp));
                        const double t1 =
                            s.start_angle +
                            s.sweep * (static_cast<double>(i + 1) / static_cast<double>(nSamp));
                        const Vec2 p0 = ellPt(t0);
                        const Vec2 p1 = ellPt(t1);
                        code(os, 0, "LINE");
                        code(os, 8, clayer);
                        code(os, 62, color);
                        code(os, 10, p0.x);
                        code(os, 20, p0.y);
                        code(os, 30, z);
                        code(os, 11, p1.x);
                        code(os, 21, p1.y);
                        code(os, 31, z);
                    }
                    continue;
                }

                code(os, 0, "LINE");
                code(os, 8, clayer);
                code(os, 62, color);
                code(os, 10, a.x);
                code(os, 20, a.y);
                code(os, 30, z);
                code(os, 11, b.x);
                code(os, 21, b.y);
                code(os, 31, z);
            }
            writeOpenGapMarker(os, c, lay, z, frame);
        }
    }
}

void writeDxfStream(std::ostream& os, const SliceResult& result) {
    os << std::setprecision(17);
    code(os, 0, "SECTION");
    code(os, 2, "HEADER");
    code(os, 9, "$ACADVER");
    code(os, 1, "AC1014");
    code(os, 0, "ENDSEC");

    code(os, 0, "SECTION");
    code(os, 2, "TABLES");
    code(os, 0, "TABLE");
    code(os, 2, "LAYER");
    code(os, 70, static_cast<int>(result.layers.size() * 3 + 1));
    code(os, 0, "LAYER");
    code(os, 2, "0");
    code(os, 70, 0);
    code(os, 62, 7);
    code(os, 6, "CONTINUOUS");
    for (size_t i = 0; i < result.layers.size(); ++i) {
        const std::string base = layerName(i, result.layers[i].z);
        code(os, 0, "LAYER");
        code(os, 2, base);
        code(os, 70, 0);
        code(os, 62, 7);
        code(os, 6, "CONTINUOUS");
        code(os, 0, "LAYER");
        code(os, 2, base + "_OPEN");
        code(os, 70, 0);
        code(os, 62, 6);  // magenta
        code(os, 6, "CONTINUOUS");
        code(os, 0, "LAYER");
        code(os, 2, base + "_OPEN_GAP");
        code(os, 70, 0);
        code(os, 62, 2);  // yellow
        code(os, 6, "DASHED");
    }
    code(os, 0, "ENDTAB");
    code(os, 0, "ENDSEC");

    code(os, 0, "SECTION");
    code(os, 2, "ENTITIES");
    writeEntities(os, result);
    code(os, 0, "ENDSEC");
    code(os, 0, "EOF");
}

}  // namespace

void writeDxfFile(const SliceResult& result, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    writeDxfStream(f, result);
}

void writeAllLayerDxf(const SliceResult& result, const std::string& path) {
    if (endsWithDxf(path)) {
        writeDxfFile(result, path);
        return;
    }
    ensureDir(path);
    for (size_t i = 0; i < result.layers.size(); ++i) {
        SliceResult one;
        one.unit = result.unit;
        one.tolerance = result.tolerance;
        one.normal = result.normal;
        one.layers.push_back(result.layers[i]);
        std::ostringstream name;
        name << path << "/layer_" << std::setw(4) << std::setfill('0') << i << ".dxf";
        writeDxfFile(one, name.str());
    }
}

void writeOpenLayerDxf(const SliceResult& result,
                       const std::vector<std::pair<double, int>>& unclosed_layers,
                       const std::string& dir) {
    if (unclosed_layers.empty()) return;
    ensureDir(dir);
    for (size_t li = 0; li < result.layers.size(); ++li) {
        const double z = result.layers[li].z;
        int open_count = 0;
        for (const auto& entry : unclosed_layers) {
            if (std::abs(entry.first - z) <= result.tolerance) {
                open_count = entry.second;
                break;
            }
        }
        if (open_count <= 0) continue;

        SliceResult one;
        one.unit = result.unit;
        one.tolerance = result.tolerance;
        one.normal = result.normal;
        one.layers.push_back(result.layers[li]);

        std::ostringstream name;
        name << dir << "/open_z" << std::fixed << std::setprecision(4) << z << "_n" << open_count
             << ".dxf";
        writeDxfFile(one, name.str());
    }
}

namespace {

std::string seedSliceLayerName(size_t i, double z) {
    std::ostringstream os;
    os << "Z_" << std::setw(4) << std::setfill('0') << i << "_" << std::fixed << std::setprecision(4)
       << z;
    return os.str();
}

int seedPointColor(SeedPointKind kind) {
    switch (kind) {
    case SeedPointKind::Raw:
        return 3;  // green
    case SeedPointKind::Kept:
        return 5;  // blue
    case SeedPointKind::Dropped:
        return 1;  // red
    }
    return 7;
}

const char* seedPointLayerSuffix(SeedPointKind kind) {
    switch (kind) {
    case SeedPointKind::Raw:
        return "_RAW";
    case SeedPointKind::Kept:
        return "_KEPT";
    case SeedPointKind::Dropped:
        return "_DROP";
    }
    return "";
}

void writeSeedPoint(std::ostream& os, const std::string& lay, int color, double z, const Vec3& p,
                    const SliceFrame& frame) {
    const Vec2 q = frame.toXY(p);
    code(os, 0, "POINT");
    code(os, 8, lay);
    code(os, 62, color);
    code(os, 10, q.x);
    code(os, 20, q.y);
    code(os, 30, z);
}

void writeSeedChain(std::ostream& os, const std::string& lay, int color, double z,
                    const std::vector<Vec3>& pts, const SliceFrame& frame) {
    if (pts.size() < 2) return;
    auto xy = [&](const Vec3& p) { return frame.toXY(p); };
    for (size_t i = 1; i < pts.size(); ++i) {
        const Vec2 a = xy(pts[i - 1]);
        const Vec2 b = xy(pts[i]);
        code(os, 0, "LINE");
        code(os, 8, lay);
        code(os, 62, color);
        code(os, 10, a.x);
        code(os, 20, a.y);
        code(os, 30, z);
        code(os, 11, b.x);
        code(os, 21, b.y);
        code(os, 31, z);
    }
}

// Layer-name prefix: empty for single-slice files (short CAD layers); Z_####_z for combined.
void writeSeedEntitiesOne(std::ostream& os, const SeedLayer& sl, size_t li,
                          const std::string& base, const SliceFrame& frame) {
    (void)li;
    const double z = sl.z;
    const std::string sep = base.empty() ? "" : "_";

    for (SeedPointKind kind : {SeedPointKind::Raw, SeedPointKind::Dropped}) {
        const char* suf = seedPointLayerSuffix(kind);  // "_RAW" / "_DROP"
        const std::string lay = base.empty() ? std::string(suf + 1) : (base + suf);
        const int color = seedPointColor(kind);
        for (const SeedPoint& sp : sl.points) {
            if (sp.kind != kind) continue;
            writeSeedPoint(os, lay, color, z, sp.p, frame);
        }
    }

    const std::string chain_raw = base + sep + "CHAIN_RAW";
    for (const auto& ch : sl.seed_chains) writeSeedChain(os, chain_raw, 3, z, ch, frame);

    const std::string chain_kept = base + sep + "CHAIN_KEPT";
    const std::string kept_all = base + sep + "KEPT";
    std::map<int, int> face_chain_idx;
    for (size_t ci = 0; ci < sl.kept_chains.size(); ++ci) {
        const auto& ch = sl.kept_chains[ci];
        writeSeedChain(os, chain_kept, 5, z, ch, frame);

        const int face_id =
            ci < sl.kept_chain_face_ids.size() ? sl.kept_chain_face_ids[ci] : -1;
        if (face_id >= 0) {
            const int local = face_chain_idx[face_id]++;
            std::ostringstream kept_lay;
            kept_lay << base << sep << "KEPT_F" << face_id << "_C" << local;
            for (const Vec3& p : ch) writeSeedPoint(os, kept_lay.str(), 5, z, p, frame);
        }
        for (const Vec3& p : ch) writeSeedPoint(os, kept_all, 5, z, p, frame);
    }
}

void collectSeedLayerDefs(const SeedLayer& sl, const std::string& base,
                          std::vector<std::pair<std::string, int>>& layer_defs) {
    const std::string sep = base.empty() ? "" : "_";
    layer_defs.push_back({base + (base.empty() ? "RAW" : "_RAW"), 3});
    layer_defs.push_back({base + (base.empty() ? "KEPT" : "_KEPT"), 5});
    layer_defs.push_back({base + (base.empty() ? "DROP" : "_DROP"), 1});
    layer_defs.push_back({base + (base.empty() ? "CHAIN_RAW" : "_CHAIN_RAW"), 3});
    layer_defs.push_back({base + (base.empty() ? "CHAIN_KEPT" : "_CHAIN_KEPT"), 5});
    std::map<int, int> face_chain_idx;
    for (size_t ci = 0; ci < sl.kept_chains.size(); ++ci) {
        const int face_id =
            ci < sl.kept_chain_face_ids.size() ? sl.kept_chain_face_ids[ci] : -1;
        if (face_id < 0) continue;
        const int local = face_chain_idx[face_id]++;
        std::ostringstream kept_lay;
        kept_lay << base << sep << "KEPT_F" << face_id << "_C" << local;
        layer_defs.push_back({kept_lay.str(), 5});
    }
}

void writeSeedDxfStreamOne(std::ostream& os, const SliceResult& meta, const SeedLayer& sl,
                           size_t li, bool short_layers) {
    os << std::setprecision(17);
    auto addLayer = [&](const std::string& name, int color) {
        code(os, 0, "LAYER");
        code(os, 2, name);
        code(os, 70, 0);
        code(os, 62, color);
        code(os, 6, "CONTINUOUS");
    };

    const std::string base = short_layers ? std::string() : seedSliceLayerName(li, sl.z);

    code(os, 0, "SECTION");
    code(os, 2, "HEADER");
    code(os, 9, "$ACADVER");
    code(os, 1, "AC1014");
    code(os, 0, "ENDSEC");

    code(os, 0, "SECTION");
    code(os, 2, "TABLES");
    std::vector<std::pair<std::string, int>> layer_defs;
    layer_defs.push_back({"0", 7});
    collectSeedLayerDefs(sl, base, layer_defs);

    code(os, 0, "TABLE");
    code(os, 2, "LAYER");
    code(os, 70, static_cast<int>(layer_defs.size()));
    for (const auto& [name, color] : layer_defs) addLayer(name, color);

    code(os, 0, "ENDTAB");
    code(os, 0, "ENDSEC");

    code(os, 0, "SECTION");
    code(os, 2, "ENTITIES");
    writeSeedEntitiesOne(os, sl, li, base, makeSliceFrame(meta.normal));
    code(os, 0, "ENDSEC");
    code(os, 0, "EOF");
}

void writeSeedDxfStream(std::ostream& os, const SliceResult& result) {
    // Combined multi-slice file (legacy): prefixed CAD layers per Z.
    os << std::setprecision(17);
    auto addLayer = [&](const std::string& name, int color) {
        code(os, 0, "LAYER");
        code(os, 2, name);
        code(os, 70, 0);
        code(os, 62, color);
        code(os, 6, "CONTINUOUS");
    };

    code(os, 0, "SECTION");
    code(os, 2, "HEADER");
    code(os, 9, "$ACADVER");
    code(os, 1, "AC1014");
    code(os, 0, "ENDSEC");

    code(os, 0, "SECTION");
    code(os, 2, "TABLES");
    std::vector<std::pair<std::string, int>> layer_defs;
    layer_defs.push_back({"0", 7});
    for (size_t li = 0; li < result.seed_layers.size(); ++li) {
        const std::string base = seedSliceLayerName(li, result.seed_layers[li].z);
        collectSeedLayerDefs(result.seed_layers[li], base, layer_defs);
    }

    code(os, 0, "TABLE");
    code(os, 2, "LAYER");
    code(os, 70, static_cast<int>(layer_defs.size()));
    for (const auto& [name, color] : layer_defs) addLayer(name, color);

    code(os, 0, "ENDTAB");
    code(os, 0, "ENDSEC");

    code(os, 0, "SECTION");
    code(os, 2, "ENTITIES");
    const SliceFrame frame = makeSliceFrame(result.normal);
    for (size_t li = 0; li < result.seed_layers.size(); ++li) {
        const std::string base = seedSliceLayerName(li, result.seed_layers[li].z);
        writeSeedEntitiesOne(os, result.seed_layers[li], li, base, frame);
    }
    code(os, 0, "ENDSEC");
    code(os, 0, "EOF");
}

}  // namespace

void writeSeedDxfFile(const SliceResult& result, const std::string& path) {
    // Directory (or non-.dxf path) → one seed_XXXX.dxf per slice (fast to open).
    // Path ending in .dxf → single combined file (legacy).
    if (!endsWithDxf(path)) {
        ensureDir(path);
        for (size_t i = 0; i < result.seed_layers.size(); ++i) {
            std::ostringstream name;
            name << path << "/seed_" << std::setw(4) << std::setfill('0') << i << ".dxf";
            std::ofstream f(name.str());
            if (!f) throw std::runtime_error("cannot write " + name.str());
            writeSeedDxfStreamOne(f, result, result.seed_layers[i], i, /*short_layers=*/true);
        }
        return;
    }
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    writeSeedDxfStream(f, result);
}

void writeOpenLayerSeedDxf(const SliceResult& result,
                           const std::vector<std::pair<double, int>>& unclosed_layers,
                           const std::string& dir) {
    ensureDir(dir);
    if (unclosed_layers.empty() || result.seed_layers.empty()) return;
    const double ztol = std::max(result.tolerance, 1e-6);
    for (size_t li = 0; li < result.seed_layers.size(); ++li) {
        const SeedLayer& sl = result.seed_layers[li];
        int open_count = 0;
        for (const auto& entry : unclosed_layers) {
            if (std::abs(entry.first - sl.z) <= ztol) {
                open_count = entry.second;
                break;
            }
        }
        if (open_count <= 0) continue;

        std::ostringstream name;
        name << dir << "/seed_" << std::setw(4) << std::setfill('0') << li << "_z" << std::fixed
             << std::setprecision(4) << sl.z << "_n" << open_count << ".dxf";
        std::ofstream f(name.str());
        if (!f) throw std::runtime_error("cannot write " + name.str());
        writeSeedDxfStreamOne(f, result, sl, li, /*short_layers=*/true);
    }
}

void writeSeedSummaryFile(const SliceResult& result, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << std::setprecision(9);
    for (size_t li = 0; li < result.seed_layers.size(); ++li) {
        const SeedLayer& sl = result.seed_layers[li];
        f << "z=" << sl.z << "\n";
        f << "face_stats: id raw kept kept_chains segments neighbor_splits\n";
        for (const FaceSeedStats& s : sl.face_stats) {
            f << s.face_id << " " << s.raw << " " << s.kept << " " << s.kept_chains << " "
              << s.segments << " " << s.neighbor_splits;
            if (s.uvmarch_fallback) f << " uvmarch";
            f << "\n";
        }
        f << "\npoints: face_id kind x y z\n";
        for (const SeedPoint& sp : sl.points) {
            const char* kind = sp.kind == SeedPointKind::Raw      ? "raw"
                               : sp.kind == SeedPointKind::Kept ? "kept"
                                                                : "drop";
            f << sp.face_id << " " << kind << " " << sp.p.x << " " << sp.p.y << " " << sp.p.z
              << "\n";
        }
        f << "\nraw_chains (" << sl.seed_chains.size() << ")\n";
        for (size_t ci = 0; ci < sl.seed_chains.size(); ++ci) {
            f << "  chain[" << ci << "] n=" << sl.seed_chains[ci].size() << "\n";
            for (const Vec3& p : sl.seed_chains[ci]) {
                f << "    " << p.x << " " << p.y << " " << p.z << "\n";
            }
        }
        f << "\nkept_chains (" << sl.kept_chains.size() << ")\n";
        std::map<int, int> face_chain_idx;
        for (size_t ci = 0; ci < sl.kept_chains.size(); ++ci) {
            const int face_id =
                ci < sl.kept_chain_face_ids.size() ? sl.kept_chain_face_ids[ci] : -1;
            const int local = face_id >= 0 ? face_chain_idx[face_id]++ : static_cast<int>(ci);
            f << "  face " << face_id << " chain[" << local << "] n=" << sl.kept_chains[ci].size()
              << "\n";
            for (const Vec3& p : sl.kept_chains[ci]) {
                f << "    " << p.x << " " << p.y << " " << p.z << "\n";
            }
        }
        f << "\n";
    }
}

}  // namespace brepslicer
