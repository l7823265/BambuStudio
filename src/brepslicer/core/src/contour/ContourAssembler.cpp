#include <contour/ContourAssembler.h>
#include <intersect/BSplineFit.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
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

Vec3 segStart(const Segment& s) {
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) return s.ctrl_pts.front();
    return s.start;
}
Vec3 segEnd(const Segment& s) {
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) return s.ctrl_pts.back();
    return s.end;
}

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

bool sameCircleSeg(const Segment& a, const Segment& b, double tol) {
    if (a.type != SegmentType::Arc || b.type != SegmentType::Arc) return false;
    return dist(a.center, b.center) <= tol && std::abs(a.radius - b.radius) <= tol;
}

bool isAngleWrap(double a) {
    a = wrapTwoPi(a);
    return a < 1e-3 || a > kTwoPi - 1e-3;
}

// Merge arcs split at 0/360 (e.g. after analytic intersection snap) into one arc.
void coalesceWrapSplitArcs(Contour& c, double tol) {
    if (c.segments.size() < 2) return;
    std::vector<Segment> out;
    out.reserve(c.segments.size());
    for (const Segment& s : c.segments) {
        if (!out.empty() && s.type == SegmentType::Arc && out.back().type == SegmentType::Arc &&
            sameCircleSeg(out.back(), s, tol) && dist(segEnd(out.back()), segStart(s)) <= tol) {
            const double tail = out.back().start_angle + out.back().sweep;
            if (isAngleWrap(tail) || isAngleWrap(s.start_angle)) {
                const bool sameDir =
                    (out.back().sweep >= 0.0 && s.sweep >= 0.0) ||
                    (out.back().sweep < 0.0 && s.sweep < 0.0);
                if (sameDir) {
                    out.back().sweep += s.sweep;
                    out.back().end = s.end;
                    continue;
                }
            }
        }
        out.push_back(s);
    }
    c.segments = std::move(out);
}

void nestAndOrientContours(std::vector<Contour>& contours, const SliceFrame& frame) {
    const int n = static_cast<int>(contours.size());
    if (n == 0) return;
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

    for (int i = 0; i < n; ++i) {
        int d = 0;
        int p = parent[i];
        int guard = 0;
        while (p >= 0 && guard++ < n + 2) {
            ++d;
            p = parent[p];
        }
        const double a = contourSignedArea(contours[i], frame);
        const bool wantCcw = (d % 2 == 0);
        if (wantCcw && a < 0) reverseContour(contours[i]);
        if (!wantCcw && a > 0) reverseContour(contours[i]);
        contours[i].orientation = (d % 2 == 0) ? "outer" : "inner";
        if (parent[i] >= 0) contours[i].parent = parent[i];
    }
}

void setSegStart(Segment& s, const Vec3& p) {
    s.start = p;
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) s.ctrl_pts.front() = p;
}

void setSegEnd(Segment& s, const Vec3& p) {
    s.end = p;
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) s.ctrl_pts.back() = p;
}

void weldChainJoints(std::vector<Segment>& chain, double tol) {
    if (chain.size() < 2) return;
    for (size_t i = 0; i + 1 < chain.size(); ++i) {
        Segment& a = chain[i];
        Segment& b = chain[i + 1];
        const Vec3 tail = segEnd(a);
        const Vec3 head = segStart(b);
        if (dist(tail, head) > tol) continue;
        const Vec3 m{(tail.x + head.x) * 0.5, (tail.y + head.y) * 0.5, (tail.z + head.z) * 0.5};
        setSegEnd(a, m);
        setSegStart(b, m);
    }
    const Vec3 tail = segEnd(chain.back());
    const Vec3 head = segStart(chain.front());
    if (dist(tail, head) <= tol) {
        const Vec3 m{(tail.x + head.x) * 0.5, (tail.y + head.y) * 0.5, (tail.z + head.z) * 0.5};
        setSegEnd(chain.back(), m);
        setSegStart(chain.front(), m);
    }
}

bool tryMergeOpen(Contour& a, Contour& b, double stitch) {
    if (a.closed || b.closed || a.solid_id != b.solid_id) return false;
    auto match = [&](const Vec3& x, const Vec3& y) { return dist(x, y) <= stitch; };

    const Vec3 aHead = segStart(a.segments.front());
    const Vec3 aTail = segEnd(a.segments.back());
    const Vec3 bHead = segStart(b.segments.front());
    const Vec3 bTail = segEnd(b.segments.back());

    if (match(aTail, bHead)) {
        a.segments.insert(a.segments.end(), b.segments.begin(), b.segments.end());
        a.coplanar = a.coplanar || b.coplanar;
        b.segments.clear();
        return true;
    }
    if (match(aTail, bTail)) {
        reverseContour(b);
        a.segments.insert(a.segments.end(), b.segments.begin(), b.segments.end());
        a.coplanar = a.coplanar || b.coplanar;
        b.segments.clear();
        return true;
    }
    if (match(aHead, bTail)) {
        a.segments.insert(a.segments.begin(), b.segments.begin(), b.segments.end());
        a.coplanar = a.coplanar || b.coplanar;
        b.segments.clear();
        return true;
    }
    if (match(aHead, bHead)) {
        reverseContour(b);
        a.segments.insert(a.segments.begin(), b.segments.begin(), b.segments.end());
        a.coplanar = a.coplanar || b.coplanar;
        b.segments.clear();
        return true;
    }
    return false;
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

    const double stitch = std::max(opt.tolerance, opt.stitch_tolerance);
    auto sameArcFrame = [&](const Segment& a, const Segment& b) {
        if (a.type != SegmentType::Arc || b.type != SegmentType::Arc) return true;
        return dist(a.center, b.center) <= stitch && std::abs(a.radius - b.radius) <= stitch;
    };
    auto onBothCircles = [&](const Vec3& p, const Segment& a, const Segment& b) {
        if (a.type != SegmentType::Arc || b.type != SegmentType::Arc) return true;
        if (sameArcFrame(a, b)) return true;
        const double tol = stitch * 10.0;
        return std::abs(dist(p, a.center) - a.radius) <= tol &&
               std::abs(dist(p, b.center) - b.radius) <= tol;
    };
    auto sameGeom = [&](const Segment& a, const Segment& b) {
        const bool fwd = dist(segStart(a), segStart(b)) <= stitch && dist(segEnd(a), segEnd(b)) <= stitch;
        const bool rev = dist(segStart(a), segEnd(b)) <= stitch && dist(segEnd(a), segStart(b)) <= stitch;
        if (!fwd && !rev) return false;
        return sameArcFrame(a, b);
    };
    auto canStitch = [&](const Segment& from, const Vec3& joint, const Segment& to) {
        return onBothCircles(joint, from, to);
    };
    for (size_t i = 0; i < segs.size(); ++i) {
        if (used[i]) continue;
        for (size_t j = i + 1; j < segs.size(); ++j) {
            if (used[j] || segs[j].closed_loop) continue;
            if (segs[i].solid_id != segs[j].solid_id) continue;
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
            const Vec3 head = segStart(chain.front());
            const Vec3 tail = segEnd(chain.back());
            for (size_t j = 0; j < segs.size(); ++j) {
                if (used[j] || !sameBody(seed, j)) continue;
                Segment cand = segs[j].geom;
                if (match(tail, segStart(cand)) && canStitch(chain.back(), tail, cand)) {
                    chain.push_back(cand);
                    used[j] = 1;
                    coplanar = coplanar || segs[j].coplanar;
                    grew = true;
                    break;
                }
                if (match(tail, segEnd(cand)) && canStitch(chain.back(), tail, cand)) {
                    reverseSegment(cand);
                    chain.push_back(cand);
                    used[j] = 1;
                    coplanar = coplanar || segs[j].coplanar;
                    grew = true;
                    break;
                }
                if (match(head, segEnd(cand)) && canStitch(cand, head, chain.front())) {
                    chain.insert(chain.begin(), cand);
                    used[j] = 1;
                    coplanar = coplanar || segs[j].coplanar;
                    grew = true;
                    break;
                }
                if (match(head, segStart(cand)) && canStitch(cand, head, chain.front())) {
                    reverseSegment(cand);
                    chain.insert(chain.begin(), cand);
                    used[j] = 1;
                    coplanar = coplanar || segs[j].coplanar;
                    grew = true;
                    break;
                }
            }
        }

        weldChainJoints(chain, stitch);

        const bool closed =
            !chain.empty() && dist(segStart(chain.front()), segEnd(chain.back())) <= stitch;
        if (!closed) {
            const double gap =
                chain.empty() ? 0.0 : dist(segStart(chain.front()), segEnd(chain.back()));
            if (gap > 0.0 && gap <= std::max(stitch * 10.0, 0.5)) {
                Segment line;
                line.type = SegmentType::Line;
                line.start = segEnd(chain.back());
                line.end = segStart(chain.front());
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

    const double mergeTol = std::max(stitch * 10.0, 1e-3);
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < contours.size(); ++i) {
            if (contours[i].segments.empty()) continue;
            for (size_t j = i + 1; j < contours.size(); ++j) {
                if (contours[j].segments.empty()) continue;
                if (tryMergeOpen(contours[i], contours[j], mergeTol)) {
                    const bool closed =
                        dist(segStart(contours[i].segments.front()),
                             segEnd(contours[i].segments.back())) <= stitch;
                    contours[i].closed = closed;
                    merged = true;
                    break;
                }
            }
            if (merged) break;
        }
    }
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [](const Contour& c) { return c.segments.empty(); }),
                   contours.end());

    for (Contour& c : contours) weldChainJoints(c.segments, stitch);

    for (Contour& c : contours) coalesceWrapSplitArcs(c, stitch);

    // Nesting + orientation per solid (multi-solid STEP: each body has its own outer/inner tree).
    std::map<int, std::vector<Contour>> bySolid;
    for (Contour& c : contours) bySolid[c.solid_id].push_back(std::move(c));
    contours.clear();
    for (auto& kv : bySolid) {
        nestAndOrientContours(kv.second, frame);
        const int base = static_cast<int>(contours.size());
        for (Contour& c : kv.second) {
            if (c.parent) *c.parent += base;
            contours.push_back(std::move(c));
        }
    }

    layer.contours = std::move(contours);
    return layer;
}

}  // namespace brepslicer
