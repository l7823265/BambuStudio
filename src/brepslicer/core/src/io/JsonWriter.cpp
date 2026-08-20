#include <brepslicer/Slicer.h>

#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace brepslicer {
namespace {

std::string fmt(double v) {
    std::ostringstream os;
    os << std::setprecision(17) << v;
    return os.str();
}

void writeVec(std::ostream& os, const Vec3& v) {
    os << "[" << fmt(v.x) << ", " << fmt(v.y) << ", " << fmt(v.z) << "]";
}

void writeSegment(std::ostream& os, const Segment& s, int indent) {
    const std::string pad(indent, ' ');
    os << pad << "{ ";
    if (s.type == SegmentType::Line) {
        os << "\"type\": \"line\", \"start\": ";
        writeVec(os, s.start);
        os << ", \"end\": ";
        writeVec(os, s.end);
    } else if (s.type == SegmentType::Arc) {
        os << "\"type\": \"arc\", \"center\": ";
        writeVec(os, s.center);
        os << ", \"normal\": ";
        writeVec(os, s.normal);
        os << ", \"radius\": " << fmt(s.radius);
        os << ", \"start_angle\": " << fmt(s.start_angle);
        os << ", \"sweep\": " << fmt(s.sweep);
    } else if (s.type == SegmentType::Ellipse) {
        os << "\"type\": \"ellipse\", \"center\": ";
        writeVec(os, s.center);
        os << ", \"major_axis\": ";
        writeVec(os, s.major_axis);
        os << ", \"a\": " << fmt(s.radius);
        os << ", \"b\": " << fmt(s.radius_b);
        os << ", \"start_param\": " << fmt(s.start_angle);
        os << ", \"sweep\": " << fmt(s.sweep);
    } else {
        os << "\"type\": \"bspline\", \"degree\": " << s.degree << ", \"knots\": [";
        for (size_t i = 0; i < s.knots.size(); ++i) {
            if (i) os << ", ";
            os << fmt(s.knots[i]);
        }
        os << "], \"ctrl_pts\": [";
        for (size_t i = 0; i < s.ctrl_pts.size(); ++i) {
            if (i) os << ", ";
            writeVec(os, s.ctrl_pts[i]);
        }
        os << "], \"weights\": [";
        for (size_t i = 0; i < s.weights.size(); ++i) {
            if (i) os << ", ";
            os << fmt(s.weights[i]);
        }
        os << "], \"fit_error\": " << fmt(s.fit_error);
    }
    os << " }";
}

}  // namespace

std::string writeJson(const SliceResult& result) {
    std::ostringstream os;
    os << "{\n";
    os << "  \"unit\": \"" << result.unit << "\",\n";
    os << "  \"tolerance\": " << fmt(result.tolerance) << ",\n";
    os << "  \"normal\": ";
    writeVec(os, result.normal);
    os << ",\n";
    os << "  \"layers\": [\n";
    for (size_t li = 0; li < result.layers.size(); ++li) {
        const Layer& layer = result.layers[li];
        os << "    {\n";
        os << "      \"z\": " << fmt(layer.z) << ",\n";
        os << "      \"contours\": [\n";
        for (size_t ci = 0; ci < layer.contours.size(); ++ci) {
            const Contour& c = layer.contours[ci];
            os << "        {\n";
            os << "          \"orientation\": \"" << c.orientation << "\",\n";
            os << "          \"parent\": ";
            if (c.parent) os << *c.parent;
            else os << "null";
            os << ",\n";
            os << "          \"closed\": " << (c.closed ? "true" : "false") << ",\n";
            if (c.coplanar) os << "          \"coplanar\": true,\n";
            os << "          \"segments\": [\n";
            for (size_t si = 0; si < c.segments.size(); ++si) {
                writeSegment(os, c.segments[si], 12);
                if (si + 1 < c.segments.size()) os << ",";
                os << "\n";
            }
            os << "          ]\n";
            os << "        }";
            if (ci + 1 < layer.contours.size()) os << ",";
            os << "\n";
        }
        os << "      ]\n";
        os << "    }";
        if (li + 1 < result.layers.size()) os << ",";
        os << "\n";
    }
    os << "  ],\n";
    os << "  \"topology_events\": [\n";
    for (size_t i = 0; i < result.topology_events.size(); ++i) {
        const TopologyEvent& e = result.topology_events[i];
        os << "    { \"z\": " << fmt(e.z) << ", \"event\": \"" << e.event << "\"";
        if (e.has_point) {
            os << ", \"point\": ";
            writeVec(os, e.point);
        }
        if (e.contour_id >= 0) os << ", \"contour_id\": " << e.contour_id;
        if (!e.from.empty()) {
            os << ", \"from\": [";
            for (size_t k = 0; k < e.from.size(); ++k) {
                if (k) os << ", ";
                os << e.from[k];
            }
            os << "]";
        }
        if (!e.to.empty()) {
            os << ", \"to\": [";
            for (size_t k = 0; k < e.to.size(); ++k) {
                if (k) os << ", ";
                os << e.to[k];
            }
            os << "]";
        }
        os << " }";
        if (i + 1 < result.topology_events.size()) os << ",";
        os << "\n";
    }
    os << "  ]";
    if (!result.logs.empty()) {
        os << ",\n  \"logs\": [\n";
        for (size_t i = 0; i < result.logs.size(); ++i) {
            os << "    \"" << result.logs[i] << "\"";
            if (i + 1 < result.logs.size()) os << ",";
            os << "\n";
        }
        os << "  ]";
    }
    os << "\n}\n";
    return os.str();
}

void writeJsonFile(const SliceResult& result, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("cannot write " + path);
    f << writeJson(result);
}

}  // namespace brepslicer
