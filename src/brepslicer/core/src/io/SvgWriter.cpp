#include <brepslicer/Slicer.h>
#include <geom/GeomUtil.h>
#include <intersect/BSplineFit.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

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

struct Box2 {
    double xmin = 1e300, ymin = 1e300, xmax = -1e300, ymax = -1e300;
    void add(double x, double y) {
        xmin = std::min(xmin, x);
        ymin = std::min(ymin, y);
        xmax = std::max(xmax, x);
        ymax = std::max(ymax, y);
    }
};

}  // namespace

void writeLayerSvg(const Layer& layer, const Vec3& normal, const std::string& path) {
    const SliceFrame frame = makeSliceFrame(normal);
    Box2 b;
    auto xy = [&](const Vec3& p) { return frame.toXY(p); };

    for (const Contour& c : layer.contours) {
        for (const Segment& s : c.segments) {
            const Vec2 a = xy(s.start);
            const Vec2 e = xy(s.end);
            b.add(a.x, a.y);
            b.add(e.x, e.y);
            if (s.type == SegmentType::Arc) {
                const Vec2 cxy = xy(s.center);
                b.add(cxy.x - s.radius, cxy.y - s.radius);
                b.add(cxy.x + s.radius, cxy.y + s.radius);
            } else if (s.type == SegmentType::BSpline) {
                std::vector<Vec3> smp;
                bsplineSample(s, 32, smp);
                for (const Vec3& p : smp) {
                    const Vec2 q = xy(p);
                    b.add(q.x, q.y);
                }
            }
        }
    }
    if (b.xmax < b.xmin) {
        b.xmin = 0;
        b.ymin = 0;
        b.xmax = 1;
        b.ymax = 1;
    }
    const double pad = std::max(1.0, 0.05 * std::max(b.xmax - b.xmin, b.ymax - b.ymin));
    b.xmin -= pad;
    b.ymin -= pad;
    b.xmax += pad;
    b.ymax += pad;
    const double w = b.xmax - b.xmin;
    const double h = b.ymax - b.ymin;

    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << std::setprecision(12);
    f << "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"" << b.xmin << " " << -b.ymax << " "
      << w << " " << h << "\">\n";
    f << "  <g transform=\"scale(1,-1)\">\n";
    f << "    <rect x=\"" << b.xmin << "\" y=\"" << b.ymin << "\" width=\"" << w << "\" height=\""
      << h << "\" fill=\"white\"/>\n";

    for (const Contour& c : layer.contours) {
        if (c.segments.empty()) continue;
        const char* color = (c.orientation == "inner") ? "#c0392b" : "#1f6feb";
        f << "    <path fill=\"none\" stroke=\"" << color << "\" stroke-width=\""
          << (0.01 * std::max(w, h)) << "\" d=\"";
        const Vec2 p0 = xy(c.segments.front().start);
        f << "M " << p0.x << " " << p0.y;
        for (const Segment& s : c.segments) {
            const Vec2 pe = xy(s.end);
            if (s.type == SegmentType::Arc) {
                const double large = (std::abs(s.sweep) > kPi) ? 1 : 0;
                // scale(1,-1) below mirrors Y and reverses visual sweep direction, so
                // invert the SVG sweep-flag or minor arcs render as open major "C"s.
                const double sweepFlag = (s.sweep < 0) ? 1 : 0;
                f << " A " << s.radius << " " << s.radius << " 0 " << large << " " << sweepFlag
                  << " " << pe.x << " " << pe.y;
            } else if (s.type == SegmentType::Ellipse) {
                // Approximate ellipse as polyline in the slice frame (rare in SVG preview).
                const int n = 32;
                const Vec3 maj = (length(s.major_axis) > 1e-15) ? normalized(s.major_axis) : frame.x;
                const Vec3 minv = cross(s.normal, maj);
                for (int i = 1; i <= n; ++i) {
                    const double t = s.start_angle + s.sweep * (static_cast<double>(i) / n);
                    const Vec3 p = s.center + maj * (s.radius * std::cos(t)) +
                                   minv * (s.radius_b * std::sin(t));
                    const Vec2 q = xy(p);
                    f << " L " << q.x << " " << q.y;
                }
            } else if (s.type == SegmentType::BSpline) {
                std::vector<Vec3> smp;
                bsplineSample(s, 48, smp);
                for (size_t k = 1; k < smp.size(); ++k) {
                    const Vec2 q = xy(smp[k]);
                    f << " L " << q.x << " " << q.y;
                }
            } else {
                f << " L " << pe.x << " " << pe.y;
            }
        }
        if (c.closed) f << " Z";
        f << "\"/>\n";
    }
    f << "  </g>\n";
    f << "  <text x=\"" << b.xmin + pad * 0.2 << "\" y=\"" << -b.ymax + pad * 0.6
      << "\" font-size=\"" << pad * 0.5 << "\" fill=\"#333\">z=" << layer.z << "</text>\n";
    f << "</svg>\n";
}

void writeAllLayerSvg(const SliceResult& result, const std::string& dir) {
    ensureDir(dir);
    for (size_t i = 0; i < result.layers.size(); ++i) {
        std::ostringstream name;
        name << dir << "/layer_" << std::setw(4) << std::setfill('0') << i << ".svg";
        writeLayerSvg(result.layers[i], result.normal, name.str());
    }
}

}  // namespace brepslicer
