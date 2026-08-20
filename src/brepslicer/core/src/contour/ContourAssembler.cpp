#include <contour/ContourAssembler.h>
#include <intersect/BSplineFit.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <utility>

namespace brepslicer {
namespace {

void reverseSegment(Segment& s) {
    std::swap(s.start, s.end);
    if (s.type == SegmentType::Arc || s.type == SegmentType::Ellipse) {
        s.start_angle = s.start_angle + s.sweep;
        s.sweep = -s.sweep;
    }
    if (s.type == SegmentType::BSpline) reverseBSpline(s);
    else std::reverse(s.ctrl_pts.begin(), s.ctrl_pts.end());
}

void reverseContour(Contour& c) {
    std::reverse(c.segments.begin(), c.segments.end());
    for (Segment& s : c.segments) reverseSegment(s);
}

Vec2 to2(const Vec3& p, const SliceFrame& frame) { return frame.toXY(p); }

std::vector<Vec2> discretize(const Contour& c, const SliceFrame& frame) {
    std::vector<Vec2> pts;
    auto pushArc = [&](const Segment& s, int n) {
        for (int i = 0; i <= n; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(n);
            const double a = s.start_angle + t * s.sweep;
            const Vec3 p = s.center + frame.x * (s.radius * std::cos(a)) +
                           frame.y * (s.radius * std::sin(a));
            pts.push_back(to2(p, frame));
        }
    };
    for (const Segment& s : c.segments) {
        if (s.type == SegmentType::Arc) {
            const int n = std::max(24, static_cast<int>(std::abs(s.sweep) * 16.0));
            pushArc(s, n);
        } else if (s.type == SegmentType::Ellipse) {
            const int n = 64;
            const Vec3 maj = normalized(s.major_axis);
            const Vec3 minv = cross(s.normal, maj);
            for (int i = 0; i <= n; ++i) {
                const double t = s.start_angle + (static_cast<double>(i) / n) * s.sweep;
                const Vec3 p = s.center + maj * (s.radius * std::cos(t)) +
                               minv * (s.radius_b * std::sin(t));
                pts.push_back(to2(p, frame));
            }
        } else if (s.type == SegmentType::BSpline) {
            std::vector<Vec3> smp;
            bsplineSample(s, 48, smp);
            for (const Vec3& p : smp) pts.push_back(to2(p, frame));
        } else {
            pts.push_back(to2(s.start, frame));
            pts.push_back(to2(s.end, frame));
        }
    }
    return pts;
}

int windingNumber(const Vec2& p, const std::vector<Vec2>& poly) {
    if (poly.size() < 3) return 0;
    int wn = 0;
    for (size_t i = 0, n = poly.size(); i < n; ++i) {
        const Vec2 a = poly[i];
        const Vec2 b = poly[(i + 1) % n];
        if (a.y <= p.y) {
            if (b.y > p.y) {
                const double cross = (b.x - a.x) * (p.y - a.y) - (p.x - a.x) * (b.y - a.y);
                if (cross > 0) ++wn;
            }
        } else if (b.y <= p.y) {
            const double cross = (b.x - a.x) * (p.y - a.y) - (p.x - a.x) * (b.y - a.y);
            if (cross < 0) --wn;
        }
    }
    return wn;
}

bool contourContainsPoint(const Contour& outer, const Vec3& p, const SliceFrame& frame) {
    const std::vector<Vec2> poly = discretize(outer, frame);
    const Vec2 q = to2(p, frame);
    return windingNumber(q, poly) != 0;
}

Vec3 segmentMid(const Segment& s, const SliceFrame& frame) {
    if (s.type == SegmentType::Arc) {
        const double a = s.start_angle + 0.5 * s.sweep;
        return s.center + frame.x * (s.radius * std::cos(a)) + frame.y * (s.radius * std::sin(a));
    }
    if (s.type == SegmentType::BSpline) return bsplineEval(s, 0.5);
    return {(s.start.x + s.end.x) * 0.5, (s.start.y + s.end.y) * 0.5, (s.start.z + s.end.z) * 0.5};
}

Contour makeContour(std::vector<Segment> segs, int solid, int shell, bool coplanar, bool closed) {
    Contour c;
    c.segments = std::move(segs);
    c.solid_id = solid;
    c.shell_id = shell;
    c.coplanar = coplanar;
    c.closed = closed;
    return c;
}

}  // namespace

Vec3 contourStart(const Contour& c) {
    if (c.segments.empty()) return {};
    return c.segments.front().start;
}

Vec3 contourEnd(const Contour& c) {
    if (c.segments.empty()) return {};
    return c.segments.back().end;
}

double contourSignedArea(const Contour& c, const SliceFrame& frame) {
    double area = 0;
    for (const Segment& s : c.segments) {
        if (s.type == SegmentType::BSpline) {
            std::vector<Vec3> smp;
            bsplineSample(s, 48, smp);
            for (size_t i = 0; i + 1 < smp.size(); ++i) {
                const Vec2 a = to2(smp[i], frame);
                const Vec2 b = to2(smp[i + 1], frame);
                area += 0.5 * (a.x * b.y - b.x * a.y);
            }
            continue;
        }
        const Vec2 a = to2(s.start, frame);
        const Vec2 b = to2(s.end, frame);
        area += 0.5 * (a.x * b.y - b.x * a.y);
        if (s.type == SegmentType::Arc) {
            area += 0.5 * s.radius * s.radius * (s.sweep - std::sin(s.sweep));
        } else if (s.type == SegmentType::Ellipse) {
            area += 0.5 * s.radius * s.radius_b * (s.sweep - std::sin(s.sweep));
        }
    }
    return area;
}

Layer assembleLayer(double z, std::vector<RawSegment> segs, const SliceFrame& frame,
                    const SliceOptions& opt, AssembleStats& stats) {
    Layer layer;
    layer.z = z;
    stats = {};

    std::vector<Contour> contours;

    std::vector<char> used(segs.size(), 0);
    for (size_t i = 0; i < segs.size(); ++i) {
        if (segs[i].degenerate) used[i] = 1;
    }

    const double stitch = opt.tolerance;
    auto sameGeom = [&](const Segment& a, const Segment& b) {
        return (dist(a.start, b.start) <= stitch && dist(a.end, b.end) <= stitch) ||
               (dist(a.start, b.end) <= stitch && dist(a.end, b.start) <= stitch);
    };
    for (size_t i = 0; i < segs.size(); ++i) {
        if (used[i]) continue;
        for (size_t j = i + 1; j < segs.size(); ++j) {
            if (used[j] || segs[j].closed_loop) continue;
            if (sameGeom(segs[i].geom, segs[j].geom)) {
                segs[i].coplanar = segs[i].coplanar || segs[j].coplanar;
                used[j] = 1;  // shared-edge / non-manifold
            }
        }
    }

    for (size_t i = 0; i < segs.size(); ++i) {
        if (used[i]) continue;
        if (segs[i].closed_loop) {
            used[i] = 1;
            contours.push_back(makeContour({segs[i].geom}, segs[i].solid_id, segs[i].shell_id,
                                           segs[i].coplanar, true));
        }
    }

    auto sameBody = [&](size_t a, size_t b) {
        return segs[a].solid_id == segs[b].solid_id;
    };
    auto match = [&](const Vec3& a, const Vec3& b) { return dist(a, b) <= stitch; };

    for (size_t seed = 0; seed < segs.size(); ++seed) {
        if (used[seed]) continue;
        used[seed] = 1;
        std::vector<Segment> chain{segs[seed].geom};
        int solid = segs[seed].solid_id;
        int shell = segs[seed].shell_id;
        bool coplanar = segs[seed].coplanar;
        bool grew = true;
        while (grew) {
            grew = false;
            const Vec3 head = chain.front().start;
            const Vec3 tail = chain.back().end;
            for (size_t j = 0; j < segs.size(); ++j) {
                if (used[j] || !sameBody(seed, j)) continue;
                Segment cand = segs[j].geom;
                if (match(tail, cand.start)) {
                    chain.push_back(cand);
                    used[j] = 1;
                    coplanar = coplanar || segs[j].coplanar;
                    grew = true;
                    break;
                }
                if (match(tail, cand.end)) {
                    reverseSegment(cand);
                    chain.push_back(cand);
                    used[j] = 1;
                    coplanar = coplanar || segs[j].coplanar;
                    grew = true;
                    break;
                }
                if (match(head, cand.end)) {
                    chain.insert(chain.begin(), cand);
                    used[j] = 1;
                    coplanar = coplanar || segs[j].coplanar;
                    grew = true;
                    break;
                }
                if (match(head, cand.start)) {
                    reverseSegment(cand);
                    chain.insert(chain.begin(), cand);
                    used[j] = 1;
                    coplanar = coplanar || segs[j].coplanar;
                    grew = true;
                    break;
                }
            }
        }

        const bool closed = !chain.empty() && dist(chain.front().start, chain.back().end) <= stitch;
        if (!closed) {
            const double gap = chain.empty() ? 0.0 : dist(chain.front().start, chain.back().end);
            if (gap > 0.0 && gap <= stitch * 10.0 && gap < 1e-3) {
                Segment line;
                line.type = SegmentType::Line;
                line.start = chain.back().end;
                line.end = chain.front().start;
                chain.push_back(line);
                ++stats.bridged;
                stats.max_bridge = std::max(stats.max_bridge, gap);
                contours.push_back(makeContour(std::move(chain), solid, shell, coplanar, true));
            } else {
                ++stats.open_leftover;
                contours.push_back(makeContour(std::move(chain), solid, shell, coplanar, false));
            }
        } else {
            contours.push_back(makeContour(std::move(chain), solid, shell, coplanar, true));
        }
    }

    // Nesting + orientation (viewed along +n: outer CCW, inner CW).
    const int n = static_cast<int>(contours.size());
    std::vector<int> parent(n, -1);
    std::vector<double> absArea(n, 0);
    for (int i = 0; i < n; ++i) {
        absArea[i] = std::abs(contourSignedArea(contours[i], frame));
    }
    for (int i = 0; i < n; ++i) {
        if (!contours[i].closed || contours[i].segments.empty()) continue;
        const Vec3 probe = segmentMid(contours[i].segments.front(), frame);
        int best = -1;
        double bestA = std::numeric_limits<double>::max();
        for (int j = 0; j < n; ++j) {
            if (i == j || !contours[j].closed) continue;
            if (absArea[j] <= absArea[i] + 1e-18) continue;
            if (contourContainsPoint(contours[j], probe, frame) && absArea[j] < bestA) {
                bestA = absArea[j];
                best = j;
            }
        }
        parent[i] = best;
    }

    std::vector<int> depth(n, 0);
    for (int i = 0; i < n; ++i) {
        int d = 0;
        int p = parent[i];
        int guard = 0;
        while (p >= 0 && guard++ < n + 2) {
            ++d;
            p = parent[p];
        }
        depth[i] = d;
        const double a = contourSignedArea(contours[i], frame);
        const bool wantCcw = (d % 2 == 0);
        if (wantCcw && a < 0) reverseContour(contours[i]);
        if (!wantCcw && a > 0) reverseContour(contours[i]);
        contours[i].orientation = (d % 2 == 0) ? "outer" : "inner";
        if (parent[i] >= 0) contours[i].parent = parent[i];
    }

    layer.contours = std::move(contours);
    return layer;
}

}  // namespace brepslicer
