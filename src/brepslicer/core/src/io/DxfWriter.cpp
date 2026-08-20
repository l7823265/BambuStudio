#include <brepslicer/Slicer.h>
#include <geom/GeomUtil.h>
#include <intersect/BSplineFit.h>

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
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

void writeEntities(std::ostream& os, const SliceResult& result) {
    const SliceFrame frame = makeSliceFrame(result.normal);
    auto xy = [&](const Vec3& p) { return frame.toXY(p); };

    for (size_t li = 0; li < result.layers.size(); ++li) {
        const Layer& layer = result.layers[li];
        const std::string lay = layerName(li, layer.z);
        const double z = layer.z;

        for (const Contour& c : layer.contours) {
            const int color = (c.orientation == "inner") ? 1 : 5;  // red / blue
            for (const Segment& s : c.segments) {
                const Vec2 a = xy(s.start);
                const Vec2 b = xy(s.end);

                if (s.type == SegmentType::Arc) {
                    const Vec2 ctr = xy(s.center);
                    if (std::abs(std::abs(s.sweep) - kTwoPi) <= 1e-8) {
                        code(os, 0, "CIRCLE");
                        code(os, 8, lay);
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
                    code(os, 8, lay);
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
                    code(os, 8, lay);
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
                    const Vec2 ctr = xy(s.center);
                    const Vec2 maj = xy(s.center + s.major_axis * s.radius);
                    code(os, 0, "ELLIPSE");
                    code(os, 8, lay);
                    code(os, 62, color);
                    code(os, 10, ctr.x);
                    code(os, 20, ctr.y);
                    code(os, 30, z);
                    code(os, 11, maj.x - ctr.x);
                    code(os, 21, maj.y - ctr.y);
                    code(os, 31, 0.0);
                    code(os, 40, (s.radius > 0) ? (s.radius_b / s.radius) : 1.0);
                    double p0 = s.start_angle;
                    double p1 = s.start_angle + s.sweep;
                    if (s.sweep < 0) std::swap(p0, p1);
                    code(os, 41, p0);
                    code(os, 42, p1);
                    continue;
                }

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
    code(os, 70, static_cast<int>(result.layers.size() + 1));
    code(os, 0, "LAYER");
    code(os, 2, "0");
    code(os, 70, 0);
    code(os, 62, 7);
    code(os, 6, "CONTINUOUS");
    for (size_t i = 0; i < result.layers.size(); ++i) {
        code(os, 0, "LAYER");
        code(os, 2, layerName(i, result.layers[i].z));
        code(os, 70, 0);
        code(os, 62, 7);
        code(os, 6, "CONTINUOUS");
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

}  // namespace brepslicer
