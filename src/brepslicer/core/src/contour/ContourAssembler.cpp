#include <contour/ContourAssembler.h>
#include <intersect/BSplineFit.h>
#include <topo/SolidAdjacency.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <utility>

namespace brepslicer {
namespace {

constexpr double kThinWallMin = 1e-3;  // DESIGN: keep distinct contours below this offset
// Hard invent-bridge / adjacency / stub join length cap (mm). Gaps longer than this
// are never invented as straight chords — topology still required above designSnap.
constexpr double kMaxInventBridge = 0.15;
// Free/NURBS tip ↔ analytic tip after fitting-edge drift (ROBOT_4 z≈278 ~0.27 mm).
// Tip replace onto the analytic endpoint only — not an invented chord.
constexpr double kFittingEdgeSnap = 0.5;

bool thinWallConflict(const Contour& a, const Contour& b, const SliceFrame& frame,
                      bool aAtHead, bool bAtHead, double joinGap);

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

void setSegStart(Segment& s, const Vec3& p);
void setSegEnd(Segment& s, const Vec3& p);

Vec2 to2(const Vec3& p, const SliceFrame& frame) { return frame.toXY(p); }

Vec3 segStart(const Segment& s) {
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) return s.ctrl_pts.front();
    return s.start;
}
Vec3 segEnd(const Segment& s) {
    if (s.type == SegmentType::BSpline && !s.ctrl_pts.empty()) return s.ctrl_pts.back();
    return s.end;
}

double pointChordDist(const Vec3& p, const Vec3& a, const Vec3& b) {
    const Vec3 ab{b.x - a.x, b.y - a.y, b.z - a.z};
    const double L2 = length2(ab);
    if (L2 <= 1e-18) return dist(p, a);
    double t = dot({p.x - a.x, p.y - a.y, p.z - a.z}, ab) / L2;
    t = std::max(0.0, std::min(1.0, t));
    return dist(p, {a.x + ab.x * t, a.y + ab.y * t, a.z + ab.z * t});
}

double distPointToSegmentGeom(const Vec3& p, const Segment& s) {
    double best = std::min(dist(p, segStart(s)), dist(p, segEnd(s)));
    best = std::min(best, pointChordDist(p, segStart(s), segEnd(s)));
    if (s.type == SegmentType::BSpline && s.ctrl_pts.size() >= 2) {
        for (size_t i = 0; i + 1 < s.ctrl_pts.size(); ++i)
            best = std::min(best, pointChordDist(p, s.ctrl_pts[i], s.ctrl_pts[i + 1]));
    }
    return best;
}

double approxContourLen(const Contour& c) {
    double L = 0;
    for (const Segment& s : c.segments) {
        if (s.type == SegmentType::Arc) {
            L += std::abs(s.radius * s.sweep);
        } else if (s.type == SegmentType::Ellipse) {
            // Chord-length sample — closed ellipse arcs have start≈end so tip chord is useless.
            const Vec3 maj = normalized(s.major_axis);
            const Vec3 minv = cross(s.normal, maj);
            constexpr int n = 16;
            Vec3 prev = s.center + maj * (s.radius * std::cos(s.start_angle)) +
                        minv * (s.radius_b * std::sin(s.start_angle));
            for (int i = 1; i <= n; ++i) {
                const double t =
                    s.start_angle + s.sweep * (static_cast<double>(i) / static_cast<double>(n));
                const Vec3 p = s.center + maj * (s.radius * std::cos(t)) +
                               minv * (s.radius_b * std::sin(t));
                L += dist(prev, p);
                prev = p;
            }
        } else if (s.type == SegmentType::BSpline && s.ctrl_pts.size() >= 2) {
            for (size_t i = 0; i + 1 < s.ctrl_pts.size(); ++i)
                L += dist(s.ctrl_pts[i], s.ctrl_pts[i + 1]);
        } else {
            L += dist(segStart(s), segEnd(s));
        }
    }
    return L;
}

// UVMatch / fillet crumbs that close on themselves (~0.3–1 mm). They poison
// dropOpenHuggingClosed: real body opens whose tips sit near the crumb get deleted
// (ROBOT_7 z≈521 layer_0027 — only four micro-bsplines left).
void dropMicroClosedLoops(std::vector<Contour>& contours, double minLen) {
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [&](const Contour& c) {
                                      if (!c.closed || c.segments.empty()) return false;
                                      return approxContourLen(c) < minLen;
                                  }),
                   contours.end());
}

// Drop tip crumbs and short opens that hug an already-closed contour (cross-face /
// cross-solid fillet duplicates just outside stitch).
void dropOpenHuggingClosed(std::vector<Contour>& contours, double stitch) {
    // Ignore micro closed crumbs as hug targets (see dropMicroClosedLoops).
    const double minClosedLen = 2.0;
    auto onClosedTip = [&](const Vec3& p) {
        for (const Contour& c : contours) {
            if (!c.closed || approxContourLen(c) < minClosedLen) continue;
            for (const Segment& s : c.segments) {
                if (dist(p, segStart(s)) <= stitch || dist(p, segEnd(s)) <= stitch) return true;
            }
        }
        return false;
    };
    auto distToClosedGeom = [&](const Vec3& p) {
        double best = std::numeric_limits<double>::infinity();
        for (const Contour& c : contours) {
            if (!c.closed || approxContourLen(c) < minClosedLen) continue;
            for (const Segment& s : c.segments) best = std::min(best, distPointToSegmentGeom(p, s));
        }
        return best;
    };
    const double hug = std::max(0.35, stitch * 40.0);
    // Near-tangent torus nick: short open whose tips sit just off a closed loop (~0.5–1 mm).
    const double graze_hug = std::max(1.0, hug);
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [&](const Contour& c) {
                                      if (c.closed || c.segments.empty()) return false;
                                      const double L = approxContourLen(c);
                                      if (L <= 0.5) return true;
                                      const Vec3 a = segStart(c.segments.front());
                                      const Vec3 b = segEnd(c.segments.back());
                                      if (L <= 1.0 && (onClosedTip(a) || onClosedTip(b))) return true;
                                      // Blind-hole U: tips sit on both ends of a short mouth edge —
                                      // keep for splice/repair (ROBOT_7 z≈603).
                                      auto isBlindHoleU = [&]() {
                                          if (L < 8.0 || c.segments.size() < 2) return false;
                                          for (const Contour& host : contours) {
                                              if (&host == &c || !host.closed) continue;
                                              if (host.solid_id != c.solid_id) continue;
                                              for (const Segment& s : host.segments) {
                                                  const Vec3 ha = segStart(s), hb = segEnd(s);
                                                  const double el = dist(ha, hb);
                                                  if (el < 1.0 || el > 12.0) continue;
                                                  const bool fwd =
                                                      dist(a, ha) <= 0.5 && dist(b, hb) <= 0.5;
                                                  const bool rev =
                                                      dist(a, hb) <= 0.5 && dist(b, ha) <= 0.5;
                                                  if ((fwd || rev) && L > el * 1.2) return true;
                                              }
                                          }
                                          return false;
                                      };
                                      if (isBlindHoleU()) return false;
                                      // Open whose tips both lie on/near a closed loop — mate/fillet
                                      // duplicate (any length; R20 bridged mid-gap left a 15 mm offset).
                                      if (distToClosedGeom(a) <= hug && distToClosedGeom(b) <= hug)
                                          return true;
                                      // Almost-tangent torus: short chord open next to closed body.
                                      if (L <= 4.0 && c.segments.size() <= 2 &&
                                          distToClosedGeom(a) <= graze_hug &&
                                          distToClosedGeom(b) <= graze_hug)
                                          return true;
                                      // Leftover ear after dual stub-join (ROBOT_7 z≈521: 2-seg
                                      // open ~3 mm off the closed body once the main loop closed).
                                      if (L <= 5.0 && c.segments.size() <= 3 &&
                                          distToClosedGeom(a) <= 3.5 && distToClosedGeom(b) <= 3.5)
                                          return true;
                                      // Do NOT drop blind-hole U opens (tips on a short mouth
                                      // chord) — splice/repair still need them (z≈603).
                                      return false;
                                  }),
                   contours.end());
}

// Short open fragments that never meet any other contour tip — typically //cylinder
// generators emitted past the trimmed wire (ROBOT_4 z≈202). Dropping them leaves the
// already-closed body and clears false "unclosed layer" reports.
void dropFloatingOpenOrphans(std::vector<Contour>& contours, double stitch) {
    const double join = std::max(2.0, 80.0 * stitch);
    auto tipNearAnyOther = [&](const Contour& self, const Vec3& tip) {
        for (const Contour& o : contours) {
            if (&o == &self || o.segments.empty()) continue;
            if (o.solid_id != self.solid_id) continue;
            if (dist(tip, segStart(o.segments.front())) <= join) return true;
            if (dist(tip, segEnd(o.segments.back())) <= join) return true;
            if (o.closed) {
                for (const Segment& s : o.segments) {
                    if (dist(tip, segStart(s)) <= join || dist(tip, segEnd(s)) <= join)
                        return true;
                }
            }
        }
        return false;
    };
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [&](const Contour& c) {
                                      if (c.closed || c.segments.empty()) return false;
                                      const double L = approxContourLen(c);
                                      if (L > 12.0 || c.segments.size() > 2) return false;
                                      const Vec3 a = segStart(c.segments.front());
                                      const Vec3 b = segEnd(c.segments.back());
                                      return !tipNearAnyOther(c, a) && !tipNearAnyOther(c, b);
                                  }),
                   contours.end());
}

// Blind-hole walls left as two opens whose tips sit on a short outer edge and whose
// free tips face each other across the pocket: replace the edge with wall→bridge→wall
// (ROBOT_7 z≈603.35 — F1 walls otherwise vanish into a flat chord / digon).
void repairBlindHoleNotches(std::vector<Contour>& contours, double stitch, AssembleStats* stats) {
    const double join = std::max(stitch, 0.05);
    const double maxEdge = 8.0;
    const double maxCap = kMaxInventBridge;  // hard invent-bridge cap
    bool progressed = true;
    while (progressed) {
        progressed = false;
        for (size_t hi = 0; hi < contours.size(); ++hi) {
            Contour& host = contours[hi];
            // Closed outer, or long open that already contains the mouth chord.
            if (host.segments.size() < 3) continue;
            if (!host.closed && approxContourLen(host) < 30.0) continue;
            for (size_t e = 0; e < host.segments.size(); ++e) {
                const Vec3 a = segStart(host.segments[e]);
                const Vec3 b = segEnd(host.segments[e]);
                const double elen = dist(a, b);
                if (elen < join || elen > maxEdge) continue;

                int iWall = -1, jWall = -1;
                bool iAtA = true, jAtA = true;
                for (size_t oi = 0; oi < contours.size(); ++oi) {
                    if (oi == hi) continue;
                    Contour& o = contours[oi];
                    if (o.closed || o.segments.empty() || o.solid_id != host.solid_id) continue;
                    if (o.segments.size() > 4) continue;
                    if (approxContourLen(o) > 40.0) continue;
                    const Vec3 oh = segStart(o.segments.front());
                    const Vec3 ot = segEnd(o.segments.back());
                    const bool hA = dist(oh, a) <= join, hB = dist(oh, b) <= join;
                    const bool tA = dist(ot, a) <= join, tB = dist(ot, b) <= join;
                    if (!(hA || hB || tA || tB)) continue;
                    if ((hA || tA) && (hB || tB)) continue;
                    if (iWall < 0) {
                        iWall = static_cast<int>(oi);
                        iAtA = hA || tA;
                    } else if (jWall < 0 && static_cast<int>(oi) != iWall) {
                        jWall = static_cast<int>(oi);
                        jAtA = hA || tA;
                    }
                }
                if (iWall < 0 || jWall < 0 || iAtA == jAtA) continue;

                Contour& w0 = contours[static_cast<size_t>(iWall)];
                Contour& w1 = contours[static_cast<size_t>(jWall)];
                auto freeTip = [&](Contour& w, bool onA) {
                    const Vec3 oh = segStart(w.segments.front());
                    const Vec3 ot = segEnd(w.segments.back());
                    const Vec3 edgePt = onA ? a : b;
                    if (dist(oh, edgePt) <= join) return ot;
                    return oh;
                };
                const Vec3 f0 = freeTip(w0, iAtA);
                const Vec3 f1 = freeTip(w1, jAtA);
                const double cap = dist(f0, f1);
                if (cap < join || cap > maxCap) continue;
                if (approxContourLen(w0) + approxContourLen(w1) + cap < elen * 1.5) continue;

                Contour* wa = iAtA ? &w0 : &w1;
                Contour* wb = iAtA ? &w1 : &w0;
                std::vector<Segment> left = wa->segments;
                std::vector<Segment> right = wb->segments;
                if (dist(segStart(left.front()), a) > join) {
                    for (Segment& s : left) reverseSegment(s);
                    std::reverse(left.begin(), left.end());
                }
                if (dist(segStart(left.front()), a) > join) continue;
                if (dist(segEnd(right.back()), b) > join) {
                    for (Segment& s : right) reverseSegment(s);
                    std::reverse(right.begin(), right.end());
                }
                if (dist(segEnd(right.back()), b) > join) continue;

                Segment bridge;
                bridge.type = SegmentType::Line;
                bridge.start = segEnd(left.back());
                bridge.end = segStart(right.front());
                bridge.face_id = left.back().face_id;
                if (dist(bridge.start, bridge.end) > maxCap) continue;

                std::vector<Segment> insert;
                insert.insert(insert.end(), left.begin(), left.end());
                if (dist(bridge.start, bridge.end) > join) insert.push_back(bridge);
                insert.insert(insert.end(), right.begin(), right.end());
                setSegStart(insert.front(), a);
                setSegEnd(insert.back(), b);

                std::vector<Segment> rebuilt;
                rebuilt.reserve(host.segments.size() + insert.size());
                for (size_t k = 0; k < e; ++k) rebuilt.push_back(host.segments[k]);
                rebuilt.insert(rebuilt.end(), insert.begin(), insert.end());
                for (size_t k = e + 1; k < host.segments.size(); ++k)
                    rebuilt.push_back(host.segments[k]);
                host.segments = std::move(rebuilt);
                if (host.closed ||
                    dist(segStart(host.segments.front()), segEnd(host.segments.back())) <= join)
                    host.closed = true;
                wa->segments.clear();
                wb->segments.clear();
                if (stats) {
                    ++stats->bridged;
                    stats->max_bridge = std::max(stats->max_bridge, cap);
                }
                progressed = true;
                break;
            }
            if (progressed) break;
        }
        contours.erase(std::remove_if(contours.begin(), contours.end(),
                                      [](const Contour& c) { return c.segments.empty(); }),
                       contours.end());
    }
}

// Open pocket U whose tips sit on both ends of a short edge of a closed host:
// replace that edge with the U (ROBOT_7 z≈603.35: outer kept the blind-hole
// opening chord while F1 walls+cap formed a separate open U / digon cap).
void spliceOpenReplacingShortHostEdge(std::vector<Contour>& contours, double stitch) {
    const double join = std::max(stitch, 0.5);  // wall tips vs mouth chord (z≈603)
    const double maxEdge = 12.0;
    bool progressed = true;
    while (progressed) {
        progressed = false;
        for (size_t oi = 0; oi < contours.size(); ++oi) {
            Contour& open = contours[oi];
            if (open.closed || open.segments.size() < 2) continue;
            const Vec3 oh = segStart(open.segments.front());
            const Vec3 ot = segEnd(open.segments.back());
            if (dist(oh, ot) <= join) continue;

            for (size_t hi = 0; hi < contours.size(); ++hi) {
                if (hi == oi) continue;
                Contour& host = contours[hi];
                if (!host.closed || host.segments.size() < 3) continue;
                if (host.solid_id != open.solid_id) continue;

                for (size_t e = 0; e < host.segments.size(); ++e) {
                    const Vec3 a = segStart(host.segments[e]);
                    const Vec3 b = segEnd(host.segments[e]);
                    const double elen = dist(a, b);
                    if (elen < join || elen > maxEdge) continue;

                    bool fwd = dist(oh, a) <= join && dist(ot, b) <= join;
                    bool rev = dist(oh, b) <= join && dist(ot, a) <= join;
                    if (!fwd && !rev) continue;
                    // Prefer a real detour: open path should be longer than the chord.
                    if (approxContourLen(open) < elen * 1.25) continue;

                    std::vector<Segment> insert = open.segments;
                    if (rev) {
                        for (Segment& s : insert) reverseSegment(s);
                        std::reverse(insert.begin(), insert.end());
                    }
                    // Snap tips onto host edge ends.
                    setSegStart(insert.front(), a);
                    setSegEnd(insert.back(), b);

                    std::vector<Segment> rebuilt;
                    rebuilt.reserve(host.segments.size() + insert.size());
                    for (size_t k = 0; k < e; ++k) rebuilt.push_back(host.segments[k]);
                    rebuilt.insert(rebuilt.end(), insert.begin(), insert.end());
                    for (size_t k = e + 1; k < host.segments.size(); ++k)
                        rebuilt.push_back(host.segments[k]);
                    host.segments = std::move(rebuilt);
                    host.closed = true;
                    open.segments.clear();
                    progressed = true;
                    break;
                }
                if (progressed) break;
            }
            if (progressed) break;
        }
        contours.erase(std::remove_if(contours.begin(), contours.end(),
                                      [](const Contour& c) { return c.segments.empty(); }),
                       contours.end());
    }
}

// Remove invented pocket diagonals (~11 mm on ROBOT_7 z≈603) when walls+cap already
// form the U — leftover from one-hop / wrong tip pairing.
// Do NOT touch real outer chamfers (ROBOT layer_0001: ~9.9 mm 45° face diagonals);
// those are near-isosceles in Δx/Δy, while invented pocket chords span a rectangular U.
void dropPocketDiagonalChords(std::vector<Contour>& contours, double stitch) {
    const double join = std::max(stitch, 0.5);
    for (Contour& c : contours) {
        if (!c.closed || c.segments.size() < 6) continue;
        for (size_t i = 0; i < c.segments.size();) {
            Segment& s = c.segments[i];
            if (s.type != SegmentType::Line) {
                ++i;
                continue;
            }
            const Vec3 a = segStart(s), b = segEnd(s);
            const double L = dist(a, b);
            if (L < 9.0 || L > 14.0) {
                ++i;
                continue;
            }
            const double dx = std::abs(a.x - b.x), dy = std::abs(a.y - b.y);
            // Pocket sides are axis-aligned; a diagonal has both Δx and Δy large.
            if (dx > 2.0 && dy > 2.0) {
                const double aspect = std::min(dx, dy) / std::max(dx, dy);
                if (aspect >= 0.65) {
                    ++i;
                    continue;
                }
                const size_t n = c.segments.size();
                const size_t ip = (i + n - 1) % n;
                const size_t in = (i + 1) % n;
                const Vec3 p = segEnd(c.segments[ip]);
                const Vec3 q = segStart(c.segments[in]);
                // Replace diagonal with the missing axis-aligned wall through (p.x,q.y)
                // or (q.x,p.y) — the corner that completes the blind-hole U.
                const Vec3 c1{p.x, q.y, p.z};
                const Vec3 c2{q.x, p.y, p.z};
                const Vec3 corner =
                    (dist(c1, p) + dist(c1, q) <= dist(c2, p) + dist(c2, q)) ? c1 : c2;
                c.segments.erase(c.segments.begin() + static_cast<std::ptrdiff_t>(i));
                // If next seg is the mouth chord (q → corner), drop it — wall ends at corner.
                if (i < c.segments.size()) {
                    const Vec3 ns = segStart(c.segments[i]), ne = segEnd(c.segments[i]);
                    if ((dist(ns, q) <= join && dist(ne, corner) <= join) ||
                        (dist(ns, corner) <= join && dist(ne, q) <= join)) {
                        c.segments.erase(c.segments.begin() + static_cast<std::ptrdiff_t>(i));
                    }
                }
                Segment wall;
                wall.type = SegmentType::Line;
                wall.start = p;
                wall.end = corner;
                wall.face_id = c.segments[ip].face_id;
                c.segments.insert(c.segments.begin() + static_cast<std::ptrdiff_t>(i), wall);
                if (dist(corner, segStart(c.segments[(i + 1) % c.segments.size()])) > join) {
                    Segment toQ;
                    toQ.type = SegmentType::Line;
                    toQ.start = corner;
                    toQ.end = segStart(c.segments[(i + 1) % c.segments.size()]);
                    toQ.face_id = wall.face_id;
                    if (dist(toQ.start, toQ.end) > join)
                        c.segments.insert(c.segments.begin() + static_cast<std::ptrdiff_t>(i + 1),
                                          toQ);
                }
                continue;
            }
            ++i;
        }
    }
}

// Closed digon (2 segs, same endpoints both ways): false blind-hole "cap" hole from an
// untrimmed F1 edge chord self-closed with its reverse (ROBOT_7 z≈603.35 layer_0041).
// Also thin 4-seg lids (ROBOT_7 L0047: ~2.7×0.5 rectangle sitting on the outer rim).
void dropDigonClosedContours(std::vector<Contour>& contours, double stitch) {
    const double join = std::max(stitch, 0.05);
    auto isReversePair = [&](const Segment& a, const Segment& b) {
        const Vec3 a0 = segStart(a), a1 = segEnd(a);
        const Vec3 b0 = segStart(b), b1 = segEnd(b);
        const bool fwd = dist(a0, b1) <= join && dist(a1, b0) <= join;
        const bool same = dist(a0, b0) <= join && dist(a1, b1) <= join;
        if (!fwd && !same) return false;
        const double la = dist(a0, a1), lb = dist(b0, b1);
        if (la < 1e-6 || lb < 1e-6) return true;
        const double ratio = la < lb ? la / lb : lb / la;
        return ratio > 0.85;
    };
    auto isLineCurveDigon = [](const Segment& a, const Segment& b) {
        // Plane generator + cylinder/cone ellipse (or arc) share tips — real slice digon
        // (ROBOT z≈525.7255 F883∥F896). Chord lengths match so isReversePair alone is wrong.
        const bool aLine = a.type == SegmentType::Line;
        const bool bLine = b.type == SegmentType::Line;
        const bool aCurve = a.type == SegmentType::Ellipse || a.type == SegmentType::Arc;
        const bool bCurve = b.type == SegmentType::Ellipse || b.type == SegmentType::Arc;
        return (aLine && bCurve) || (bLine && aCurve);
    };
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [&](const Contour& c) {
                                      if (!c.closed) return false;
                                      if (c.segments.size() == 2) {
                                          if (isLineCurveDigon(c.segments[0], c.segments[1]))
                                              return false;
                                          return isReversePair(c.segments[0], c.segments[1]);
                                      }
                                      if (c.segments.size() != 4) return false;
                                      // Opposite sides reverse each other → thin digon lid.
                                      if (isReversePair(c.segments[0], c.segments[2]) &&
                                          isReversePair(c.segments[1], c.segments[3]))
                                          return true;
                                      // Or two consecutive reverse pairs (folded digon).
                                      if (isReversePair(c.segments[0], c.segments[1]) &&
                                          isReversePair(c.segments[2], c.segments[3]))
                                          return true;
                                      // Thin rectangular lid (~2.7×0.5 on ROBOT_7 L0047).
                                      // Do not catch real shelf ears (~7×1) or mid digons
                                      // (~8.8×1) on huapingdun z≈±11.62 — those match OCC.
                                      Vec3 mn = segStart(c.segments[0]), mx = mn;
                                      for (const Segment& s : c.segments) {
                                          for (const Vec3& p : {segStart(s), segEnd(s)}) {
                                              mn.x = std::min(mn.x, p.x);
                                              mn.y = std::min(mn.y, p.y);
                                              mn.z = std::min(mn.z, p.z);
                                              mx.x = std::max(mx.x, p.x);
                                              mx.y = std::max(mx.y, p.y);
                                              mx.z = std::max(mx.z, p.z);
                                          }
                                      }
                                      const double dx = mx.x - mn.x, dy = mx.y - mn.y;
                                      return std::min(dx, dy) <= 0.7 && std::max(dx, dy) <= 4.0 &&
                                             approxContourLen(c) <= 12.0;
                                  }),
                   contours.end());
}

// False pocket "caps" that landed as digon ears on a closed outer: consecutive
// segment + reverse (bspline↔line). Drop both; joints reconnect (ROBOT_7 L0037/L0041).
void stripEmbeddedDigonEars(std::vector<Contour>& contours, double stitch) {
    const double join = std::max(stitch, 0.05);
    for (Contour& c : contours) {
        if (!c.closed || c.segments.size() < 4) continue;
        bool progressed = true;
        while (progressed) {
            progressed = false;
            const size_t n = c.segments.size();
            if (n < 4) break;
            for (size_t i = 0; i < n; ++i) {
                const size_t j = (i + 1) % n;
                const Vec3 a0 = segStart(c.segments[i]), a1 = segEnd(c.segments[i]);
                const Vec3 b0 = segStart(c.segments[j]), b1 = segEnd(c.segments[j]);
                if (!(dist(a0, b1) <= join && dist(a1, b0) <= join)) continue;
                const double la = dist(a0, a1), lb = dist(b0, b1);
                if (la < 1.0 || lb < 1.0) continue;
                const double ratio = la < lb ? la / lb : lb / la;
                if (ratio < 0.85) continue;
                // Prefer dropping the shorter / line-looking ear of a real boundary arc.
                if (j > i) {
                    c.segments.erase(c.segments.begin() + static_cast<std::ptrdiff_t>(j));
                    c.segments.erase(c.segments.begin() + static_cast<std::ptrdiff_t>(i));
                } else {
                    c.segments.erase(c.segments.begin() + static_cast<std::ptrdiff_t>(i));
                    c.segments.erase(c.segments.begin() + static_cast<std::ptrdiff_t>(j));
                }
                progressed = true;
                break;
            }
        }
    }
}

// F1→F1 U-turn pockets that closed as separate loops sit on the outer boundary as
// false holes (ROBOT_7 z≈544.866: two ~16 mm 4-seg loops). Drop short closed loops
// whose every vertex lies near a much longer closed contour.
void dropShortClosedOnLongBoundary(std::vector<Contour>& contours, double stitch,
                                   double maxLen) {
    // Mid vertices sit ~3.5 mm off the outer (only tips coincide).
    const double hug = std::max(4.0, stitch * 50.0);
    auto distTo = [&](const Contour& host, const Vec3& p) {
        double best = std::numeric_limits<double>::infinity();
        for (const Segment& s : host.segments)
            best = std::min(best, distPointToSegmentGeom(p, s));
        return best;
    };
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [&](const Contour& c) {
                                      if (!c.closed || c.segments.empty()) return false;
                                      const double L = approxContourLen(c);
                                      if (L > maxLen || c.segments.size() > 8) return false;
                                      for (const Contour& host : contours) {
                                          if (&host == &c || !host.closed) continue;
                                          if (host.solid_id != c.solid_id) continue;
                                          if (approxContourLen(host) < L * 3.0) continue;
                                          bool allOn = true;
                                          for (const Segment& s : c.segments) {
                                              if (distTo(host, segStart(s)) > hug ||
                                                  distTo(host, segEnd(s)) > hug) {
                                                  allOn = false;
                                                  break;
                                              }
                                          }
                                          if (allOn) return true;
                                      }
                                      return false;
                                  }),
                   contours.end());
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

Segment stampFace(const RawSegment& rs) {
    Segment s = rs.geom;
    s.face_id = rs.face_id;
    s.chain_idx = rs.chain_idx >= 0 ? rs.chain_idx : rs.geom.chain_idx;
    return s;
}

// Faces are linkable if they are the same face, share an edge, or share a
// one-hop neighbor (missing slick face). Same-face must stay topo-preferred in
// the seed chain so continuous F1 walks are not stolen by a nearby F2 tip
// (ROBOT_7 z≈521). F1→F1 U-turn mistakes are repaired after chaining
// (spliceOpenOrphansAtMidJoints / stripShortClosedSubloops).
bool topoLinkableFaces(int fa, int fb, const SolidAdjacency& adj) {
    if (fa < 0 || fb < 0) return false;
    if (fa == fb) return true;
    if (adj.facesShareEdge(fa, fb)) return true;
    for (int n : adj.neighbors(fa)) {
        if (n == fb) return true;
        if (adj.facesShareEdge(n, fb)) return true;
    }
    return false;
}

// Same-face U-turn pockets embedded in an open chain (F1→F1 while F2 is a separate
// open). Strip short closed sub-loops, then tip-link the leftover opens
// (ROBOT_7 z≈544.866).
double approxSegmentLen(const Segment& s) {
    Contour tmp;
    tmp.segments = {s};
    return approxContourLen(tmp);
}

void stripShortClosedSubloops(std::vector<Contour>& contours, double stitch, double maxLoopLen) {
    const double join = std::max(stitch, 0.05);
    for (Contour& c : contours) {
        // Open chains only — closed touching loops are real pockets to extract, not delete
        // (ROBOT_7 z≈626.747: deleting the ear removed the pocket).
        if (c.closed || c.segments.size() < 3) continue;
        bool changed = true;
        while (changed && c.segments.size() >= 3) {
            changed = false;
            const size_t n = c.segments.size();
            std::vector<Vec3> node(n + 1);
            for (size_t i = 0; i < n; ++i) node[i] = segStart(c.segments[i]);
            node[n] = segEnd(c.segments.back());
            for (size_t i = 0; i + 2 <= n; ++i) {
                double plen = 0;
                for (size_t j = i + 1; j <= n; ++j) {
                    plen += approxSegmentLen(c.segments[j - 1]);
                    if (j < i + 2) continue;
                    if (plen > maxLoopLen) break;
                    if (dist(node[i], node[j]) > join) continue;
                    // Only strip true F1→F1 U-turns (all segs same face).
                    const int face0 = c.segments[i].face_id;
                    bool sameFace = face0 >= 0;
                    for (size_t k = i; k < j && sameFace; ++k) {
                        if (c.segments[k].face_id != face0) sameFace = false;
                    }
                    if (!sameFace) continue;
                    if (plen * 2.0 > approxContourLen(c)) continue;
                    c.segments.erase(c.segments.begin() + static_cast<std::ptrdiff_t>(i),
                                     c.segments.begin() + static_cast<std::ptrdiff_t>(j));
                    changed = true;
                    break;
                }
                if (changed) break;
            }
        }
        if (!c.segments.empty()) {
            const double gap = dist(segStart(c.segments.front()), segEnd(c.segments.back()));
            if (gap <= std::max(stitch, 0.06)) c.closed = true;
        }
    }
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [](const Contour& c) { return c.segments.empty(); }),
                   contours.end());
}

// Open chain whose tip gap equals the first/last segment (digon cord): that segment is
// a false pocket "cap". Drop the cord; the remaining self-touch loop becomes a notch
// on the closed outer (ROBOT_7 z≈626.747: ~5 mm bridge duplicated seg0).
void dropOpenDigonCapChords(std::vector<Contour>& contours, double stitch) {
    const double join = std::max(stitch, 0.05);
    for (Contour& c : contours) {
        if (c.closed || c.segments.size() < 4) continue;
        const Vec3 h = segStart(c.segments.front());
        const Vec3 t = segEnd(c.segments.back());
        const double tipGap = dist(h, t);
        if (tipGap <= join || tipGap > 25.0) continue;

        auto tryDropFront = [&]() -> bool {
            Segment& s0 = c.segments.front();
            if (dist(segStart(s0), h) > join || dist(segEnd(s0), t) > join) return false;
            if (std::abs(approxSegmentLen(s0) - tipGap) > std::max(0.5, 0.25 * tipGap))
                return false;
            // Require a self-touch near h after the cord (pocket notch remains).
            bool touch = false;
            for (size_t i = 1; i + 1 < c.segments.size(); ++i) {
                if (dist(segEnd(c.segments[i]), h) <= join) {
                    touch = true;
                    break;
                }
            }
            if (!touch) return false;
            c.segments.erase(c.segments.begin());
            c.closed = dist(segStart(c.segments.front()), segEnd(c.segments.back())) <= join;
            if (c.closed) {
                const Vec3 m = segStart(c.segments.front());
                setSegStart(c.segments.front(), m);
                setSegEnd(c.segments.back(), m);
            }
            return true;
        };
        auto tryDropBack = [&]() -> bool {
            Segment& sN = c.segments.back();
            if (dist(segEnd(sN), t) > join || dist(segStart(sN), h) > join) return false;
            if (std::abs(approxSegmentLen(sN) - tipGap) > std::max(0.5, 0.25 * tipGap))
                return false;
            bool touch = false;
            for (size_t i = 0; i + 1 < c.segments.size(); ++i) {
                if (dist(segStart(c.segments[i]), t) <= join) {
                    touch = true;
                    break;
                }
            }
            if (!touch) return false;
            c.segments.pop_back();
            c.closed = dist(segStart(c.segments.front()), segEnd(c.segments.back())) <= join;
            if (c.closed) {
                const Vec3 m = segStart(c.segments.front());
                setSegStart(c.segments.front(), m);
                setSegEnd(c.segments.back(), m);
            }
            return true;
        };
        if (!tryDropFront()) tryDropBack();
    }
}

// Orphan open B whose head/tail sit on mid-joints of open A (A already took F1→F1
// pockets). Rebuild A as body ∪ B using the pocket edge spurs as bridges
// (ROBOT_7 z≈544.866). Requires short closed sub-loops — arbitrary mid-span replace
// invents wrong bridges when one face has multiple chain locations (z≈574).
bool spliceOpenOrphanAtMidJoints(Contour& a, Contour& b, double stitch, AssembleStats* stats) {
    if (a.closed || b.closed || a.segments.size() < 6 || b.segments.empty()) return false;
    if (a.solid_id != b.solid_id) return false;
    if (b.segments.size() >= a.segments.size()) return false;

    auto endIndexNear = [&](const Contour& c, const Vec3& p) -> int {
        if (dist(p, segStart(c.segments.front())) <= stitch) return -1;
        for (size_t i = 0; i < c.segments.size(); ++i) {
            if (dist(p, segEnd(c.segments[i])) <= stitch) return static_cast<int>(i);
        }
        return -999;
    };

    const Vec3 aHead = segStart(a.segments.front());
    const Vec3 aTail = segEnd(a.segments.back());

    for (int orient = 0; orient < 2; ++orient) {
        Contour bUse = b;
        if (orient) reverseContour(bUse);
        const Vec3 bh = segStart(bUse.segments.front());
        const Vec3 bt = segEnd(bUse.segments.back());
        if (dist(bh, aHead) <= stitch || dist(bh, aTail) <= stitch) continue;
        if (dist(bt, aHead) <= stitch || dist(bt, aTail) <= stitch) continue;

        const int iHead = endIndexNear(a, bh);
        const int iTail = endIndexNear(a, bt);
        if (iHead < 0 || iTail < 0 || iHead == iTail) continue;

        const size_t n = a.segments.size();
        std::vector<Vec3> node(n + 1);
        for (size_t i = 0; i < n; ++i) node[i] = segStart(a.segments[i]);
        node[n] = segEnd(a.segments.back());

        auto findLoopCovering = [&](int idx, size_t& L0, size_t& L1) -> bool {
            for (size_t i = 0; i + 2 <= n; ++i) {
                double plen = 0;
                for (size_t j = i + 1; j <= n; ++j) {
                    plen += approxSegmentLen(a.segments[j - 1]);
                    if (j < i + 2) continue;
                    if (plen > 25.0) break;
                    if (dist(node[i], node[j]) > stitch) continue;
                    if (static_cast<int>(i) <= idx && idx < static_cast<int>(j)) {
                        L0 = i;
                        L1 = j;
                        return true;
                    }
                }
            }
            return false;
        };

        size_t L0 = 0, L1 = 0, H0 = 0, H1 = 0;
        if (!findLoopCovering(iTail, L0, L1)) continue;
        if (!findLoopCovering(iHead, H0, H1)) continue;
        if (L1 > H0) continue;

        std::vector<Segment> body(a.segments.begin() + static_cast<std::ptrdiff_t>(L1),
                                  a.segments.begin() + static_cast<std::ptrdiff_t>(H0));
        if (body.size() < 2) continue;

        const Vec3 bodyH = segStart(body.front());
        const Vec3 bodyT = segEnd(body.back());

        auto makeBridge = [&](const Vec3& from, const Vec3& to, size_t p0, size_t p1,
                              int faceHint) -> Segment {
            for (size_t i = p0; i < p1; ++i) {
                Segment s = a.segments[i];
                if (dist(segStart(s), from) <= stitch && dist(segEnd(s), to) <= stitch) {
                    setSegStart(s, from);
                    setSegEnd(s, to);
                    return s;
                }
                reverseSegment(s);
                if (dist(segStart(s), from) <= stitch && dist(segEnd(s), to) <= stitch) {
                    setSegStart(s, from);
                    setSegEnd(s, to);
                    return s;
                }
            }
            Segment line;
            line.type = SegmentType::Line;
            line.start = from;
            line.end = to;
            line.face_id = faceHint;
            return line;
        };

        Segment bridgeTop = makeBridge(bodyT, bh, H0, H1, body.back().face_id);
        Segment bridgeBot = makeBridge(bt, bodyH, L0, L1, body.front().face_id);
        // Reject invented long fills unless they are short tip snaps.
        // Reused pocket-edge geometry may be longer (real curve).
        auto isInvented = [](const Segment& s) {
            return s.type == SegmentType::Line && s.ctrl_pts.empty();
        };
        if (isInvented(bridgeTop) &&
            dist(segStart(bridgeTop), segEnd(bridgeTop)) > kMaxInventBridge)
            continue;
        if (isInvented(bridgeBot) &&
            dist(segStart(bridgeBot), segEnd(bridgeBot)) > kMaxInventBridge)
            continue;

        const double inventLen = std::max(dist(segStart(bridgeTop), segEnd(bridgeTop)),
                                          dist(segStart(bridgeBot), segEnd(bridgeBot)));

        std::vector<Segment> out = std::move(body);
        out.push_back(std::move(bridgeTop));
        out.insert(out.end(), bUse.segments.begin(), bUse.segments.end());
        out.push_back(std::move(bridgeBot));
        a.segments = std::move(out);
        a.closed = dist(segStart(a.segments.front()), segEnd(a.segments.back())) <=
                   kMaxInventBridge;
        if (!a.closed) {
            const double gap = dist(segStart(a.segments.front()), segEnd(a.segments.back()));
            if (gap <= kMaxInventBridge) {
                const Vec3 h = segStart(a.segments.front());
                const Vec3 t = segEnd(a.segments.back());
                const Vec3 m{(h.x + t.x) * 0.5, (h.y + t.y) * 0.5, (h.z + t.z) * 0.5};
                setSegStart(a.segments.front(), m);
                setSegEnd(a.segments.back(), m);
                a.closed = true;
            }
        }
        b.segments.clear();
        if (stats) {
            ++stats->bridged;
            stats->max_bridge = std::max(stats->max_bridge, inventLen);
        }
        return true;
    }
    return false;
}

void spliceOpenOrphansAtMidJoints(std::vector<Contour>& contours, double stitch,
                                  AssembleStats* stats) {
    bool progressed = true;
    while (progressed) {
        progressed = false;
        for (size_t i = 0; i < contours.size(); ++i) {
            if (contours[i].closed || contours[i].segments.empty()) continue;
            for (size_t j = 0; j < contours.size(); ++j) {
                if (i == j || contours[j].closed || contours[j].segments.empty()) continue;
                Contour& a = contours[i];
                Contour& b = contours[j];
                Contour& longC = (a.segments.size() >= b.segments.size()) ? a : b;
                Contour& shortC = (a.segments.size() >= b.segments.size()) ? b : a;
                if (spliceOpenOrphanAtMidJoints(longC, shortC, stitch, stats)) {
                    progressed = true;
                    break;
                }
            }
            if (progressed) break;
        }
        contours.erase(std::remove_if(contours.begin(), contours.end(),
                                      [](const Contour& c) { return c.segments.empty(); }),
                       contours.end());
    }
}

// Connect open contours whose tips meet at a shared-edge / one-hop face junction.
// Hard bridge length cap kMaxInventBridge (debug / WithSphere wrong long joins). Topology
// still required above designSnap; nothing longer than the hard cap is invented.
// designSnap: DESIGN geometry-only tip snap = min(1e-3, 10×tol); longer joins need topo.
void linkOpenByFaceAdjacency(std::vector<Contour>& contours, const SolidAdjacency* adj,
                             double stitch, double designSnap, const SliceFrame& frame,
                             AssembleStats* stats) {
    if (!adj || adj->empty()) return;
    const double snapGap = std::max(designSnap, 1e-9);  // geometry-only tip snap (D1)
    const double adjGap = kMaxInventBridge;  // hard cap — no long invented bridges
    const double closeTol = snapGap;

    auto tipFace = [](const Contour& c, bool at_head) {
        if (c.segments.empty()) return -1;
        return at_head ? c.segments.front().face_id : c.segments.back().face_id;
    };
    auto tipChain = [](const Contour& c, bool at_head) {
        if (c.segments.empty()) return -1;
        return at_head ? c.segments.front().chain_idx : c.segments.back().chain_idx;
    };

    bool progressed = true;
    while (progressed) {
        progressed = false;
        for (size_t i = 0; i < contours.size(); ++i) {
            if (contours[i].closed || contours[i].segments.empty()) continue;
            for (size_t j = i + 1; j < contours.size(); ++j) {
                if (contours[j].closed || contours[j].segments.empty()) continue;
                if (contours[i].solid_id != contours[j].solid_id) continue;

                Contour& a = contours[i];
                Contour& b = contours[j];
                struct Cand {
                    bool a_tail = true;
                    bool b_head = true;
                    double gap = 0;
                    int rank = 2;  // 0 shared-edge, 1 one-hop/same-face, 2 snap-only
                    int fa = -1;
                    int fb = -1;
                    Vec3 pa{};
                    Vec3 pb{};
                };
                std::vector<Cand> cands;
                auto consider = [&](bool a_tail, bool b_head, const Vec3& pa, const Vec3& pb, int fa,
                                    int fb, int ca, int cb) {
                    const double gap = dist(pa, pb);
                    if (gap > adjGap) return;
                    const bool shared = adj->facesShareEdge(fa, fb);
                    const bool topo = topoLinkableFaces(fa, fb, *adj);
                    if (gap > snapGap && !topo) return;
                    // Same face + same chain_idx: one curve instance — only tip snap.
                    if (fa == fb && ca >= 0 && ca == cb && gap > snapGap) return;
                    // Same face + different chain_idx: pocket wall↔cap corner only.
                    // Reject long diagonal joins (ROBOT_7 z≈603: walls bridged at 11 mm
                    // diagonal instead of 4.8 mm blind-end cap).
                    if (fa == fb && ca >= 0 && cb >= 0 && ca != cb && gap > snapGap) {
                        if (gap > 8.0) return;
                    }
                    // Long one-hop fills invent diagonal pocket chords; require a real
                    // shared edge past ~8 mm (z≈603 max bridge 11.08).
                    if (gap > 8.0 && !shared && fa != fb) return;
                    const int rank = shared ? 0 : (topo ? 1 : 2);
                    cands.push_back({a_tail, b_head, gap, rank, fa, fb, pa, pb});
                };
                consider(true, true, segEnd(a.segments.back()), segStart(b.segments.front()),
                         tipFace(a, false), tipFace(b, true), tipChain(a, false), tipChain(b, true));
                consider(true, false, segEnd(a.segments.back()), segEnd(b.segments.back()),
                         tipFace(a, false), tipFace(b, false), tipChain(a, false),
                         tipChain(b, false));
                consider(false, true, segStart(a.segments.front()), segStart(b.segments.front()),
                         tipFace(a, true), tipFace(b, true), tipChain(a, true), tipChain(b, true));
                consider(false, false, segStart(a.segments.front()), segEnd(b.segments.back()),
                         tipFace(a, true), tipFace(b, false), tipChain(a, true), tipChain(b, false));
                if (cands.empty()) continue;
                std::sort(cands.begin(), cands.end(), [](const Cand& x, const Cand& y) {
                    if (x.rank != y.rank) return x.rank < y.rank;
                    return x.gap < y.gap;
                });
                // D3: skip antiparallel / parallel thin-wall tip pairs.
                Cand c{};
                bool picked = false;
                for (const Cand& cand : cands) {
                    const bool aAtHead = !cand.a_tail;
                    const bool bAtHead = cand.b_head;
                    if (thinWallConflict(a, b, frame, aAtHead, bAtHead, cand.gap)) continue;
                    c = cand;
                    picked = true;
                    break;
                }
                if (!picked) continue;

                std::vector<Segment> left = a.segments;
                if (!c.a_tail) {
                    for (Segment& s : left) reverseSegment(s);
                    std::reverse(left.begin(), left.end());
                }
                std::vector<Segment> right = b.segments;
                if (!c.b_head) {
                    for (Segment& s : right) reverseSegment(s);
                    std::reverse(right.begin(), right.end());
                }

                if (c.gap > stitch) {
                    Segment line;
                    line.type = SegmentType::Line;
                    line.start = segEnd(left.back());
                    line.end = segStart(right.front());
                    line.face_id = (c.fa >= 0) ? c.fa : c.fb;
                    left.push_back(line);
                    if (stats) {
                        ++stats->bridged;
                        stats->max_bridge = std::max(stats->max_bridge, c.gap);
                    }
                }
                left.insert(left.end(), right.begin(), right.end());
                a.segments = std::move(left);
                a.coplanar = a.coplanar || b.coplanar;
                a.closed =
                    dist(segStart(a.segments.front()), segEnd(a.segments.back())) <= closeTol;
                b.segments.clear();
                progressed = true;
                break;
            }
            if (progressed) break;
        }
        contours.erase(std::remove_if(contours.begin(), contours.end(),
                                      [](const Contour& c) { return c.segments.empty(); }),
                       contours.end());
    }
}

// Single open loop whose head/tail are a fillet tip and a topo-neighbor chain tip
// (ROBOT_4 z≈167.121: ~0.8 mm short of shared-edge link after UVMatch constraint snap).
//
// Do NOT self-close short analytic stubs (F2↔F3) in an F1→F2→F3→F1 cycle — that makes a
// tiny false pocket and dropOpenHuggingClosed then deletes the real F1 opens (z≈293.559).
// designSnap: DESIGN geometry-only self-close; larger gaps need tip-face adjacency.
void closeOpenSelfByFaceAdjacency(std::vector<Contour>& contours, const SolidAdjacency* adj,
                                  double stitch, double designSnap, AssembleStats* stats) {
    const double snapGap = std::max(designSnap, 1e-9);
    const double adjGap = kMaxInventBridge;  // hard cap — match linkOpenByFaceAdjacency
    const double trivialClose = snapGap;

    for (Contour& c : contours) {
        if (c.closed || c.segments.size() < 2) continue;
        const Vec3 h = segStart(c.segments.front());
        const Vec3 t = segEnd(c.segments.back());
        const double gap = dist(h, t);
        if (gap > trivialClose) continue;
        if (gap > stitch) {
            const Vec3 m{(h.x + t.x) * 0.5, (h.y + t.y) * 0.5, (h.z + t.z) * 0.5};
            setSegStart(c.segments.front(), m);
            setSegEnd(c.segments.back(), m);
        }
        c.closed = true;
    }

    if (!adj || adj->empty()) return;

    // Only defer self-close when another open tip sits at this junction (not within adjGap).
    auto tipJoinsOtherOpen = [&](const Contour& self, const Vec3& tip, double gap) {
        const double reach = std::max(gap + stitch, 0.5);
        for (const Contour& o : contours) {
            if (&o == &self || o.closed || o.segments.empty()) continue;
            if (o.solid_id != self.solid_id) continue;
            if (dist(tip, segStart(o.segments.front())) <= reach) return true;
            if (dist(tip, segEnd(o.segments.back())) <= reach) return true;
        }
        return false;
    };

    for (Contour& c : contours) {
        if (c.closed || c.segments.size() < 2) continue;
        const Vec3 h = segStart(c.segments.front());
        const Vec3 t = segEnd(c.segments.back());
        const double gap = dist(h, t);
        if (gap <= trivialClose) continue;

        Segment& headSeg = c.segments.front();
        Segment& tailSeg = c.segments.back();
        const bool fillet_analytic =
            (tailSeg.type == SegmentType::BSpline && headSeg.type != SegmentType::BSpline) ||
            (headSeg.type == SegmentType::BSpline && tailSeg.type != SegmentType::BSpline);
        const bool fillet_fillet =
            headSeg.type == SegmentType::BSpline && tailSeg.type == SegmentType::BSpline;
        // BSpline→analytic tip replace may exceed invent-bridge; invented chords stay capped.
        const double allowGap = fillet_analytic ? std::max(adjGap, kFittingEdgeSnap) : adjGap;
        if (gap > allowGap) continue;

        const int fa = headSeg.face_id;
        const int fb = tailSeg.face_id;
        const int ca = headSeg.chain_idx;
        const int cb = tailSeg.chain_idx;
        // Same chain instance: only tip snap (already handled).
        if (fa == fb && ca >= 0 && ca == cb) continue;
        // Same-face pocket corner only — not the long diagonal (z≈603 ~11 mm).
        if (fa == fb && ca >= 0 && cb >= 0 && ca != cb && gap > 8.0) continue;
        if (!topoLinkableFaces(fa, fb, *adj)) continue;

        // Short pure line/line self-close → z293 false pocket. Allow only when the
        // tips are shared-edge neighbors (or different location on same face) and
        // the chain is substantial (narrow face∩plane tip gap on a real body).
        if (!fillet_analytic && !fillet_fillet) {
            const bool shared = adj->facesShareEdge(fa, fb) ||
                                (fa == fb && (ca < 0 || cb < 0 || ca != cb));
            if (!shared) continue;
            if (c.segments.size() < 4 || approxContourLen(c) < 10.0) continue;
        }

        // Another open tip is waiting on this junction — link there instead of self-close.
        if (tipJoinsOtherOpen(c, h, gap) || tipJoinsOtherOpen(c, t, gap)) continue;

        // Prefer snapping a fillet tip onto the neighbor tip.
        if (tailSeg.type == SegmentType::BSpline && headSeg.type != SegmentType::BSpline) {
            setSegEnd(tailSeg, h);
        } else if (headSeg.type == SegmentType::BSpline &&
                   tailSeg.type != SegmentType::BSpline) {
            setSegStart(headSeg, t);
        } else {
            if (gap > adjGap) continue;
            Segment line;
            line.type = SegmentType::Line;
            line.start = t;
            line.end = h;
            line.face_id = (fb >= 0) ? fb : fa;
            c.segments.push_back(line);
        }
        c.closed = true;
        if (stats) {
            ++stats->bridged;
            stats->max_bridge = std::max(stats->max_bridge, gap);
        }
    }
}

bool sameCircleSeg(const Segment& a, const Segment& b, double tol) {
    if (a.type != SegmentType::Arc || b.type != SegmentType::Arc) return false;
    return dist(a.center, b.center) <= tol && std::abs(a.radius - b.radius) <= tol &&
           almostParallel(a.normal, b.normal, std::max(tol, 1e-8));
}

bool coversAngle(const Segment& s, double a_test) {
    constexpr double eps = 1e-8;
    double d = a_test - s.start_angle;
    while (d <= -kPi) d += kTwoPi;
    while (d > kPi) d -= kTwoPi;
    if (s.sweep >= 0.0) return d >= -eps && d <= s.sweep + eps;
    return d <= eps && d >= s.sweep - eps;
}

bool arcCoversArc(const Segment& outer, const Segment& inner) {
    if (std::abs(inner.sweep) >= kTwoPi - 1e-6) return std::abs(outer.sweep) >= kTwoPi - 1e-6;
    const double a0 = inner.start_angle;
    const double a1 = inner.start_angle + 0.5 * inner.sweep;
    const double a2 = inner.start_angle + inner.sweep;
    return coversAngle(outer, a0) && coversAngle(outer, a1) && coversAngle(outer, a2);
}

bool isAngleWrap(double a) {
    a = wrapTwoPi(a);
    return a < 1e-3 || a > kTwoPi - 1e-3;
}

Segment arcFromAngles(const Segment& proto, double start_angle, double sweep, const SliceFrame& frame) {
    Segment s = proto;
    s.type = SegmentType::Arc;
    s.start_angle = start_angle;
    s.sweep = sweep;
    s.start = s.center + frame.x * (s.radius * std::cos(start_angle)) +
              frame.y * (s.radius * std::sin(start_angle));
    const double a1 = start_angle + sweep;
    s.end = s.center + frame.x * (s.radius * std::cos(a1)) + frame.y * (s.radius * std::sin(a1));
    return s;
}

void normalizeArcCcw(Segment& s) {
    if (s.type != SegmentType::Arc) return;
    if (s.sweep < 0.0) {
        s.start_angle = s.start_angle + s.sweep;
        s.sweep = -s.sweep;
        std::swap(s.start, s.end);
    }
    s.start_angle = wrapTwoPi(s.start_angle);
    if (s.sweep > kTwoPi) s.sweep = kTwoPi;
}

// Prefer coplanar dump, else lower face_id (A1 shared-edge owner).
bool preferOwner(const RawSegment& keep, const RawSegment& drop) {
    if (keep.coplanar != drop.coplanar) return keep.coplanar;
    if (keep.face_id != drop.face_id) return keep.face_id <= drop.face_id;
    return true;
}

void applyOwnerMeta(RawSegment& dst, const RawSegment& src) {
    if (src.coplanar) dst.coplanar = true;
    if (src.closed_loop) dst.closed_loop = true;
    if (!preferOwner(dst, src)) {
        dst.face_id = src.face_id;
        dst.chain_idx = src.chain_idx;
        dst.shell_id = src.shell_id;
        dst.geom.face_id = src.face_id;
        dst.geom.chain_idx = src.chain_idx;
    }
}

std::vector<std::pair<double, double>> arcToLinearParts(double start, double sweep) {
    start = wrapTwoPi(start);
    if (sweep >= kTwoPi - 1e-9) return {{0.0, kTwoPi}};
    const double end = start + sweep;
    if (end <= kTwoPi + 1e-12) return {{start, std::min(end, kTwoPi)}};
    return {{start, kTwoPi}, {0.0, end - kTwoPi}};
}

std::vector<std::pair<double, double>> mergeLinearIntervals(
    std::vector<std::pair<double, double>> parts, double gapTol) {
    if (parts.empty()) return parts;
    std::sort(parts.begin(), parts.end());
    std::vector<std::pair<double, double>> out;
    out.push_back(parts.front());
    for (size_t i = 1; i < parts.size(); ++i) {
        if (parts[i].first <= out.back().second + gapTol) {
            out.back().second = std::max(out.back().second, parts[i].second);
        } else {
            out.push_back(parts[i]);
        }
    }
    // Merge wrap across 0 if first and last touch.
    if (out.size() >= 2 && out.front().first <= gapTol &&
        out.back().second >= kTwoPi - gapTol) {
        out.front().first = out.back().first - kTwoPi;
        out.pop_back();
        // Represent as single wrapping interval via start+sweep later.
    }
    return out;
}

// Phase 1 B2: same-circle partial arc interval union (solid_id scoped).
void unionOverlayArcs(std::vector<RawSegment>& segs, const SliceFrame& frame, double tol) {
    const double match = std::max(tol, 1e-6);
    std::vector<char> used(segs.size(), 0);

    for (size_t i = 0; i < segs.size(); ++i) {
        if (used[i] || segs[i].degenerate || segs[i].geom.type != SegmentType::Arc) continue;
        std::vector<size_t> group{i};
        for (size_t j = i + 1; j < segs.size(); ++j) {
            if (used[j] || segs[j].degenerate || segs[j].geom.type != SegmentType::Arc) continue;
            if (segs[j].solid_id != segs[i].solid_id) continue;
            if (!sameCircleSeg(segs[i].geom, segs[j].geom, match)) continue;
            group.push_back(j);
        }
        if (group.size() < 2) continue;

        bool anyFull = false;
        RawSegment owner = segs[group.front()];
        std::vector<std::pair<double, double>> parts;
        const double r = std::max(segs[group.front()].geom.radius, 1e-9);
        const double gap = std::max(match / r, 1e-8);

        // Drop exact angle-subsets via coversAngle before interval union.
        std::vector<char> drop(group.size(), 0);
        for (size_t gi = 0; gi < group.size(); ++gi) {
            if (drop[gi]) continue;
            Segment ai = segs[group[gi]].geom;
            normalizeArcCcw(ai);
            for (size_t gj = 0; gj < group.size(); ++gj) {
                if (gi == gj || drop[gj]) continue;
                Segment aj = segs[group[gj]].geom;
                normalizeArcCcw(aj);
                if (arcCoversArc(ai, aj) && !arcCoversArc(aj, ai)) drop[gj] = 1;
            }
        }

        for (size_t gi = 0; gi < group.size(); ++gi) {
            if (drop[gi]) {
                used[group[gi]] = 1;
                continue;
            }
            const size_t idx = group[gi];
            Segment a = segs[idx].geom;
            normalizeArcCcw(a);
            if (segs[idx].closed_loop || std::abs(a.sweep) >= kTwoPi - 1e-6) anyFull = true;
            applyOwnerMeta(owner, segs[idx]);
            auto lp = arcToLinearParts(a.start_angle, a.sweep);
            parts.insert(parts.end(), lp.begin(), lp.end());
            used[idx] = 1;
        }
        if (parts.empty()) continue;

        if (anyFull) {
            Segment full = owner.geom;
            full.sweep = kTwoPi;
            full.start_angle = 0.0;
            full.start = full.center + frame.x * full.radius;
            full.end = full.start;
            owner.geom = full;
            owner.closed_loop = true;
            owner.degenerate = false;
            segs[group.front()] = owner;
            used[group.front()] = 0;
            continue;
        }

        parts = mergeLinearIntervals(std::move(parts), gap);
        std::vector<RawSegment> merged;
        for (const auto& pr : parts) {
            const double a0 = pr.first;
            const double a1 = pr.second;
            const double sweep = a1 - a0;
            if (sweep <= 1e-12) continue;
            // Negative a0 marks a wrap-merged interval spanning 0.
            RawSegment rs = owner;
            rs.geom = arcFromAngles(owner.geom, wrapTwoPi(a0), sweep, frame);
            rs.closed_loop = std::abs(sweep - kTwoPi) <= 1e-6;
            rs.degenerate = false;
            merged.push_back(rs);
        }
        if (merged.empty()) continue;
        segs[group.front()] = merged.front();
        used[group.front()] = 0;
        for (size_t k = 1; k < merged.size(); ++k) {
            segs.push_back(merged[k]);
            used.push_back(0);  // keep in sync with segs
        }
    }

    std::vector<RawSegment> kept;
    kept.reserve(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        if (!used[i] && !segs[i].degenerate) kept.push_back(std::move(segs[i]));
    }
    segs = std::move(kept);
}

bool collinearLines(const Segment& a, const Segment& b, double tol) {
    if (a.type != SegmentType::Line || b.type != SegmentType::Line) return false;
    const Vec3 da = a.end - a.start;
    const Vec3 db = b.end - b.start;
    const double la = length(da);
    const double lb = length(db);
    if (la <= tol || lb <= tol) return false;
    if (!almostParallel(da, db, 1e-6)) return false;
    const Vec3 n = normalized(da);
    const Vec3 d0 = b.start - a.start;
    const Vec3 perp = d0 - n * dot(d0, n);
    return length(perp) <= tol;
}

// Phase 1 B3: collinear overlapping line merge (solid_id scoped).
void unionOverlayLines(std::vector<RawSegment>& segs, double tol) {
    const double match = std::max(tol, 1e-6);
    std::vector<char> used(segs.size(), 0);

    for (size_t i = 0; i < segs.size(); ++i) {
        if (used[i] || segs[i].degenerate || segs[i].geom.type != SegmentType::Line) continue;
        std::vector<size_t> group{i};
        for (size_t j = i + 1; j < segs.size(); ++j) {
            if (used[j] || segs[j].degenerate || segs[j].geom.type != SegmentType::Line) continue;
            if (segs[j].solid_id != segs[i].solid_id) continue;
            if (!collinearLines(segs[i].geom, segs[j].geom, match)) continue;
            group.push_back(j);
        }
        if (group.size() < 2) continue;

        const Vec3 origin = segs[group.front()].geom.start;
        Vec3 axis = segs[group.front()].geom.end - segs[group.front()].geom.start;
        if (length(axis) <= match) continue;
        axis = normalized(axis);

        RawSegment owner = segs[group.front()];
        std::vector<std::pair<double, double>> ivals;
        for (size_t idx : group) {
            applyOwnerMeta(owner, segs[idx]);
            const double t0 = dot(segs[idx].geom.start - origin, axis);
            const double t1 = dot(segs[idx].geom.end - origin, axis);
            ivals.push_back({std::min(t0, t1), std::max(t0, t1)});
            used[idx] = 1;
        }
        std::sort(ivals.begin(), ivals.end());
        std::vector<std::pair<double, double>> merged;
        for (const auto& iv : ivals) {
            if (merged.empty() || iv.first > merged.back().second + match) {
                merged.push_back(iv);
            } else {
                merged.back().second = std::max(merged.back().second, iv.second);
            }
        }

        segs[group.front()] = owner;
        segs[group.front()].geom.start = origin + axis * merged.front().first;
        segs[group.front()].geom.end = origin + axis * merged.front().second;
        segs[group.front()].degenerate = false;
        used[group.front()] = 0;
        for (size_t k = 1; k < merged.size(); ++k) {
            RawSegment extra = owner;
            extra.geom.start = origin + axis * merged[k].first;
            extra.geom.end = origin + axis * merged[k].second;
            extra.degenerate = false;
            segs.push_back(extra);
            used.push_back(0);  // keep in sync with segs (extras were dropped when used[] OOB)
        }
    }

    std::vector<RawSegment> kept;
    kept.reserve(segs.size());
    for (size_t i = 0; i < segs.size(); ++i) {
        if (!used[i] && !segs[i].degenerate) kept.push_back(std::move(segs[i]));
    }
    segs = std::move(kept);
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

// Open tips lie on a free-face full circle (sphere∩plane kept as closed_loop): splice the
// gap-closing arc into the open and drop the full circle (interior remainder discarded).
void spliceOpenThroughClosedCircles(std::vector<Contour>& contours, const SliceFrame& frame,
                                    double on_tol) {
    const double onCirc = std::max(on_tol, 0.05);
    for (size_t oi = 0; oi < contours.size(); ++oi) {
        Contour& open = contours[oi];
        if (open.closed || open.segments.empty()) continue;
        const Vec3 head = segStart(open.segments.front());
        const Vec3 tail = segEnd(open.segments.back());
        const double tipGap = dist(head, tail);
        if (tipGap <= on_tol) continue;

        for (size_t ci = 0; ci < contours.size(); ++ci) {
            if (ci == oi) continue;
            Contour& circ = contours[ci];
            if (!circ.closed || circ.segments.size() != 1) continue;
            const Segment& full = circ.segments.front();
            if (full.type != SegmentType::Arc || full.radius <= 1e-9) continue;
            if (std::abs(full.sweep) < kTwoPi - 1e-3) continue;
            if (open.solid_id != circ.solid_id) continue;

            const double dh = std::abs(dist(head, full.center) - full.radius);
            const double dt = std::abs(dist(tail, full.center) - full.radius);
            if (dh > onCirc || dt > onCirc) continue;

            const double aH = wrapTwoPi(frame.angle(full.center, head));
            const double aT = wrapTwoPi(frame.angle(full.center, tail));
            double sw = aH - aT;
            while (sw <= -kPi) sw += kTwoPi;
            while (sw > kPi) sw -= kTwoPi;
            // Prefer the shorter arc for tipGap (chord ≈ tip gap).
            if (std::abs(sw) * full.radius > tipGap * 1.25 &&
                (kTwoPi - std::abs(sw)) * full.radius <= tipGap * 1.25) {
                sw = sw > 0.0 ? sw - kTwoPi : sw + kTwoPi;
            }
            if (std::abs(sw) < 1e-9) continue;

            Segment bridge = arcFromAngles(full, aT, sw, frame);
            if (dist(segStart(bridge), tail) > dist(segEnd(bridge), tail)) reverseSegment(bridge);
            setSegStart(bridge, tail);
            setSegEnd(bridge, head);
            open.segments.push_back(bridge);
            open.closed = true;
            circ.segments.clear();
            break;
        }
    }
    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [](const Contour& c) { return c.segments.empty(); }),
                   contours.end());
}

// Open chain already picked the long sphere arc (wrong way around) while both tips lie on
// that circle: replace trailing/leading arc with the short tip-to-tip bridge and close.
void repairOpenWrongCircleArc(std::vector<Contour>& contours, const SliceFrame& frame,
                              double on_tol) {
    const double onCirc = std::max(on_tol, 0.05);
    for (Contour& c : contours) {
        if (c.closed || c.segments.size() < 2) continue;

        auto tryFix = [&](bool trailing) -> bool {
            Segment& arc = trailing ? c.segments.back() : c.segments.front();
            if (arc.type != SegmentType::Arc || arc.radius <= 1e-9) return false;

            const Vec3 otherTip = trailing ? segStart(c.segments.front()) : segEnd(c.segments.back());
            const Vec3 joint = trailing ? segEnd(c.segments[c.segments.size() - 2])
                                       : segStart(c.segments[1]);
            if (std::abs(dist(otherTip, arc.center) - arc.radius) > onCirc) return false;
            if (std::abs(dist(joint, arc.center) - arc.radius) > onCirc) return false;

            const double tipGap = dist(joint, otherTip);
            if (tipGap <= on_tol) {
                // Tips already meet — drop the stray arc.
                if (trailing) c.segments.pop_back();
                else c.segments.erase(c.segments.begin());
                c.closed = true;
                return true;
            }

            const double aJ = wrapTwoPi(frame.angle(arc.center, joint));
            const double aO = wrapTwoPi(frame.angle(arc.center, otherTip));
            double sw = aO - aJ;
            while (sw <= -kPi) sw += kTwoPi;
            while (sw > kPi) sw -= kTwoPi;
            // Prefer shorter arc for tip gap.
            if (std::abs(sw) * arc.radius > tipGap * 1.25 &&
                (kTwoPi - std::abs(sw)) * arc.radius <= tipGap * 1.25) {
                sw = sw > 0.0 ? sw - kTwoPi : sw + kTwoPi;
            }
            if (std::abs(sw) < 1e-9) return false;

            Segment bridge = arcFromAngles(arc, aJ, sw, frame);
            if (dist(segStart(bridge), joint) > dist(segEnd(bridge), joint)) reverseSegment(bridge);
            setSegStart(bridge, joint);
            setSegEnd(bridge, otherTip);
            arc = bridge;
            c.closed = true;
            return true;
        };

        if (!tryFix(true)) tryFix(false);
    }
}

Vec3 segmentTangentForward(const Segment& s, const SliceFrame& frame, bool atStart) {
    if (s.type == SegmentType::Line) {
        return normalized(s.end - s.start);
    }
    if (s.type == SegmentType::Arc) {
        const Vec3 p = atStart ? segStart(s) : segEnd(s);
        Vec3 t = cross(s.normal, p - s.center);
        if (length(t) < 1e-18) t = cross(frame.n, p - s.center);
        t = normalized(t);
        if (s.sweep < 0.0) t = t * -1.0;
        return t;
    }
    if (s.type == SegmentType::BSpline && s.ctrl_pts.size() >= 2) {
        if (atStart) return normalized(s.ctrl_pts[1] - s.ctrl_pts[0]);
        const size_t n = s.ctrl_pts.size();
        return normalized(s.ctrl_pts[n - 1] - s.ctrl_pts[n - 2]);
    }
    return normalized(segEnd(s) - segStart(s));
}

Vec3 contourEndTangent(const Contour& c, const SliceFrame& frame, bool atHead) {
    if (c.segments.empty()) return {};
    if (atHead) return segmentTangentForward(c.segments.front(), frame, true);
    return segmentTangentForward(c.segments.back(), frame, false);
}

Vec3 pointBackAlong(const Contour& c, const SliceFrame& frame, bool fromHead, double travel) {
    if (c.segments.empty()) return {};
    if (fromHead) {
        const Segment& s = c.segments.front();
        const Vec3 t = segmentTangentForward(s, frame, true);
        return segStart(s) + t * travel;
    }
    const Segment& s = c.segments.back();
    const Vec3 t = segmentTangentForward(s, frame, false);
    return segEnd(s) - t * travel;
}

// D3: refuse joining antiparallel / parallel thin-wall opens across a small offset.
bool thinWallConflict(const Contour& a, const Contour& b, const SliceFrame& frame,
                      bool aAtHead, bool bAtHead, double joinGap) {
    if (joinGap <= 1e-9) return false;
    const Vec3 ta = contourEndTangent(a, frame, aAtHead);
    const Vec3 tb = contourEndTangent(b, frame, bAtHead);
    if (length(ta) < 1e-12 || length(tb) < 1e-12) return false;
    const double c = std::abs(dot(ta, tb));
    if (c < 0.85) return false;  // not (anti)parallel

    const double sample = std::max(5.0 * joinGap, 0.05);
    const Vec3 pa = pointBackAlong(a, frame, aAtHead, sample);
    const Vec3 pb = pointBackAlong(b, frame, bAtHead, sample);
    const Vec3 axis = ta;
    const Vec3 delta = pb - pa;
    const Vec3 perp = delta - axis * dot(delta, axis);
    const double offset = length(perp);
    return offset > 1e-6 && offset < kThinWallMin;
}

bool tryMergeOpen(Contour& a, Contour& b, double stitch, const SliceFrame& frame) {
    if (a.closed || b.closed || a.solid_id != b.solid_id) return false;
    auto match = [&](const Vec3& x, const Vec3& y) { return dist(x, y) <= stitch; };

    const Vec3 aHead = segStart(a.segments.front());
    const Vec3 aTail = segEnd(a.segments.back());
    const Vec3 bHead = segStart(b.segments.front());
    const Vec3 bTail = segEnd(b.segments.back());

    auto join = [&](bool aAtHead, bool bAtHead, const Vec3& pa, const Vec3& pb) {
        const double gap = dist(pa, pb);
        if (!match(pa, pb)) return false;
        if (thinWallConflict(a, b, frame, aAtHead, bAtHead, gap)) return false;
        return true;
    };

    if (join(false, true, aTail, bHead)) {
        a.segments.insert(a.segments.end(), b.segments.begin(), b.segments.end());
        a.coplanar = a.coplanar || b.coplanar;
        b.segments.clear();
        return true;
    }
    if (join(false, false, aTail, bTail)) {
        reverseContour(b);
        a.segments.insert(a.segments.end(), b.segments.begin(), b.segments.end());
        a.coplanar = a.coplanar || b.coplanar;
        b.segments.clear();
        return true;
    }
    if (join(true, false, aHead, bTail)) {
        a.segments.insert(a.segments.begin(), b.segments.begin(), b.segments.end());
        a.coplanar = a.coplanar || b.coplanar;
        b.segments.clear();
        return true;
    }
    if (join(true, true, aHead, bHead)) {
        reverseContour(b);
        a.segments.insert(a.segments.begin(), b.segments.begin(), b.segments.end());
        a.coplanar = a.coplanar || b.coplanar;
        b.segments.clear();
        return true;
    }
    return false;
}

// Fill a short missing mid-fillet span: two open tips nearly collinear / opposing within
// maxGap (R20: ~2.5 mm hole between UVMatch halves after neighbor/gap split).
// Geometry/tangent gated only — tip faces need not share a B-rep edge
// (ROBOT_7 z≈521: ~3.1 mm tip pairs after micro-ring drop).
bool tryBridgeFilletGap(Contour& a, Contour& b, double maxGap, double stitch,
                        const SliceFrame& frame, AssembleStats* stats,
                        const std::vector<Contour>* all = nullptr) {
    if (a.closed || b.closed || a.segments.empty() || b.segments.empty()) return false;

    struct Cand {
        bool aAtHead = false;
        bool bAtHead = false;
        Vec3 pa;
        Vec3 pb;
        double gap = 0;
    };
    const Vec3 aHead = segStart(a.segments.front());
    const Vec3 aTail = segEnd(a.segments.back());
    const Vec3 bHead = segStart(b.segments.front());
    const Vec3 bTail = segEnd(b.segments.back());
    const Cand cands[] = {
        {false, true, aTail, bHead, dist(aTail, bHead)},
        {false, false, aTail, bTail, dist(aTail, bTail)},
        {true, false, aHead, bTail, dist(aHead, bTail)},
        {true, true, aHead, bHead, dist(aHead, bHead)},
    };

    // Hole mouth arcs already spanning these tips (often on another closed lobe) must
    // not be replaced by an invented chord (ROBOT z≈545 F97 through top hole, gap≈2.49).
    auto arcOwnsTips = [&](const Vec3& pa, const Vec3& pb) {
        if (!all) return false;
        for (const Contour& c : *all) {
            if (&c == &a || &c == &b) continue;
            for (const Segment& s : c.segments) {
                if (s.type == SegmentType::Line) continue;
                const bool fwd =
                    dist(segStart(s), pa) <= stitch && dist(segEnd(s), pb) <= stitch;
                const bool rev =
                    dist(segStart(s), pb) <= stitch && dist(segEnd(s), pa) <= stitch;
                if (fwd || rev) return true;
            }
        }
        return false;
    };

    const Cand* best = nullptr;
    for (const Cand& c : cands) {
        if (c.gap <= stitch || c.gap > maxGap) continue;
        if (arcOwnsTips(c.pa, c.pb)) continue;
        const Vec3 ta = contourEndTangent(a, frame, c.aAtHead);
        const Vec3 tb = contourEndTangent(b, frame, c.bAtHead);
        if (length(ta) < 1e-12 || length(tb) < 1e-12) continue;
        // Outward (into the gap) directions from each tip.
        const Vec3 aOut = c.aAtHead ? (ta * -1.0) : ta;
        const Vec3 bOut = c.bAtHead ? (tb * -1.0) : tb;
        const Vec3 join = normalized(c.pb - c.pa);
        if (length(join) < 1e-12) continue;
        // Both should aim into the gap toward each other.
        const double along = -dot(aOut, bOut);
        const double aim = dot(aOut, join) + dot(bOut, join * -1.0);
        // Aligned stubs on one axis (ROBOT_7 z≈521): tips share a line; one contour's
        // forward may already point into the gap so along/aim fail.
        const Vec3 delta = c.pb - c.pa;
        const Vec3 taN = normalized(ta);
        const Vec3 tbN = normalized(tb);
        const double latA = length(delta - taN * dot(delta, taN));
        const double latB = length(delta - tbN * dot(delta, tbN));
        const double latTol = std::max(0.25, 0.1 * c.gap);
        const bool alignedStubs = std::min(latA, latB) <= latTol;
        if (!alignedStubs && (along < 0.5 || aim < 1.2)) continue;
        if (thinWallConflict(a, b, frame, c.aAtHead, c.bAtHead, c.gap)) continue;
        if (!best || c.gap < best->gap) best = &c;
    }
    if (!best) return false;

    Segment bridge;
    bridge.type = SegmentType::Line;
    bridge.start = best->pa;
    bridge.end = best->pb;

    Contour bb = b;
    if (best->aAtHead == false && best->bAtHead == true) {
        // aTail -> bHead
        a.segments.push_back(bridge);
        a.segments.insert(a.segments.end(), bb.segments.begin(), bb.segments.end());
    } else if (best->aAtHead == false && best->bAtHead == false) {
        reverseContour(bb);
        bridge.end = segStart(bb.segments.front());
        a.segments.push_back(bridge);
        a.segments.insert(a.segments.end(), bb.segments.begin(), bb.segments.end());
    } else if (best->aAtHead == true && best->bAtHead == false) {
        // bTail -> aHead: prepend b then bridge
        a.segments.insert(a.segments.begin(), bridge);
        a.segments.insert(a.segments.begin(), bb.segments.begin(), bb.segments.end());
    } else {
        reverseContour(bb);
        bridge.start = segEnd(bb.segments.back());
        bridge.end = aHead;
        a.segments.insert(a.segments.begin(), bridge);
        a.segments.insert(a.segments.begin(), bb.segments.begin(), bb.segments.end());
    }
    a.coplanar = a.coplanar || b.coplanar;
    if (a.solid_id < 0) a.solid_id = b.solid_id;
    b.segments.clear();
    if (stats) {
        ++stats->bridged;
        stats->max_bridge = std::max(stats->max_bridge, best->gap);
    }
    const double closeGap = dist(segStart(a.segments.front()), segEnd(a.segments.back()));
    a.closed = closeGap <= stitch;
    return true;
}

void bridgeFilletOpenGaps(std::vector<Contour>& contours, const SliceFrame& frame, double maxGap,
                          double stitch, AssembleStats* stats) {
    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < contours.size(); ++i) {
            if (contours[i].closed || contours[i].segments.empty()) continue;
            for (size_t j = i + 1; j < contours.size(); ++j) {
                if (contours[j].closed || contours[j].segments.empty()) continue;
                if (tryBridgeFilletGap(contours[i], contours[j], maxGap, stitch, frame, stats,
                                       &contours)) {
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
}

double pointSegDist2(const Vec2& p, const Vec2& a, const Vec2& b) {
    const double vx = b.x - a.x;
    const double vy = b.y - a.y;
    const double len2 = vx * vx + vy * vy;
    if (len2 < 1e-30) {
        const double dx = p.x - a.x;
        const double dy = p.y - a.y;
        return dx * dx + dy * dy;
    }
    double t = ((p.x - a.x) * vx + (p.y - a.y) * vy) / len2;
    t = std::max(0.0, std::min(1.0, t));
    const double dx = p.x - (a.x + t * vx);
    const double dy = p.y - (a.y + t * vy);
    return dx * dx + dy * dy;
}

double directedPolyDist(const std::vector<Vec2>& a, const std::vector<Vec2>& b) {
    if (a.empty() || b.size() < 2) return std::numeric_limits<double>::infinity();
    double maxd = 0.0;
    const size_t step = std::max<size_t>(1, a.size() / 32);
    for (size_t i = 0; i < a.size(); i += step) {
        double mind = std::numeric_limits<double>::infinity();
        for (size_t j = 0; j + 1 < b.size(); ++j) {
            mind = std::min(mind, pointSegDist2(a[i], b[j], b[j + 1]));
        }
        mind = std::min(mind, pointSegDist2(a[i], b.back(), b.front()));
        maxd = std::max(maxd, mind);
    }
    return std::sqrt(maxd);
}

bool closedLoopsSimilar(const Contour& a, const Contour& b, const SliceFrame& frame, double tol) {
    if (!a.closed || !b.closed) return false;
    if (a.solid_id != b.solid_id) return false;
    const double aa = std::abs(contourSignedArea(a, frame));
    const double ab = std::abs(contourSignedArea(b, frame));
    if (aa < 1e-18 || ab < 1e-18) return false;
    // Require near-equal area so congruent-but-displaced holes are not dropped.
    const double ratio = aa < ab ? aa / ab : ab / aa;
    if (ratio < 0.95) return false;

    const std::vector<Vec2> pa = discretize(a, frame);
    const std::vector<Vec2> pb = discretize(b, frame);
    if (pa.size() < 3 || pb.size() < 3) return false;

    auto centroid = [](const std::vector<Vec2>& p) {
        Vec2 c{};
        for (const Vec2& q : p) {
            c.x += q.x;
            c.y += q.y;
        }
        const double n = static_cast<double>(p.size());
        return Vec2{c.x / n, c.y / n};
    };
    const Vec2 ca = centroid(pa);
    const Vec2 cb = centroid(pb);
    const double cdx = ca.x - cb.x;
    const double cdy = ca.y - cb.y;
    // Overlapping duplicates share a centroid; displaced copies do not.
    if (cdx * cdx + cdy * cdy > tol * tol * 4.0) return false;

    return directedPolyDist(pa, pb) <= tol && directedPolyDist(pb, pa) <= tol;
}

int dropDuplicateClosedLoops(std::vector<Contour>& contours, const SliceFrame& frame, double tol) {
    const int n = static_cast<int>(contours.size());
    std::vector<char> drop(n, 0);
    std::vector<double> area(n, 0);
    for (int i = 0; i < n; ++i) {
        if (contours[i].closed) area[i] = std::abs(contourSignedArea(contours[i], frame));
    }
    int removed = 0;
    for (int i = 0; i < n; ++i) {
        if (drop[i] || !contours[i].closed || contours[i].segments.empty()) continue;
        for (int j = i + 1; j < n; ++j) {
            if (drop[j] || !contours[j].closed) continue;
            if (!closedLoopsSimilar(contours[i], contours[j], frame, tol)) continue;
            // Keep larger area (more complete outer); break ties by orientation outer.
            int kill = j;
            if (area[j] > area[i] * 1.01) kill = i;
            else if (std::abs(area[j] - area[i]) <= 1e-12 &&
                     contours[j].orientation == "outer" && contours[i].orientation != "outer")
                kill = i;
            drop[kill] = 1;
            ++removed;
            if (kill == i) break;
        }
    }
    if (removed == 0) return 0;
    std::vector<Contour> kept;
    kept.reserve(contours.size());
    for (int i = 0; i < n; ++i) {
        if (!drop[i]) kept.push_back(std::move(contours[i]));
    }
    contours = std::move(kept);
    return removed;
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

void fixContours(std::vector<Contour>& contours, double gap_tol) {
    std::vector<Vec3> openEnds;
    for (const Contour& c : contours) {
        if (c.closed || c.segments.empty()) continue;
        openEnds.push_back(segStart(c.segments.front()));
        openEnds.push_back(segEnd(c.segments.back()));
    }
    if (openEnds.empty()) return;

    auto nearOpenEnd = [&](const Vec3& p) {
        for (const Vec3& e : openEnds) {
            if (dist(p, e) <= gap_tol) return true;
        }
        return false;
    };

    contours.erase(std::remove_if(contours.begin(), contours.end(),
                                  [&](const Contour& c) {
                                      if (!c.closed || c.segments.empty()) return false;
                                      // Spurious pocket: every vertex sits on an open endpoint.
                                      for (const Segment& s : c.segments) {
                                          if (!nearOpenEnd(segStart(s))) return false;
                                          if (!nearOpenEnd(segEnd(s))) return false;
                                      }
                                      return true;
                                  }),
                   contours.end());
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
    // D1: DESIGN invent-bridge cap = min(1e-3, 10×tol). Geometry-only tip close uses
    // stitch (endpoint coincidence), not the old hard 0.05. Adjacency/fillet passes
    // still handle longer topo-backed gaps.
    const double designBridge = std::min(1e-3, 10.0 * opt.tolerance);
    const double closeTol = stitch;
    const double bridgeMax = designBridge;
    // Open tip merge at stitch (endpoint coincidence); D3 thinWallConflict guards walls.
    const double mergeTol = stitch;
    // Geometry-only tip snap for adjacency helpers (was 0.05).
    const double designSnap = stitch;

    // Phase 1: overlay curve union before stitch (arcs + collinear lines).
    {
        std::vector<RawSegment> live;
        live.reserve(segs.size());
        for (size_t i = 0; i < segs.size(); ++i) {
            if (!used[i]) live.push_back(std::move(segs[i]));
        }
        unionOverlayArcs(live, frame, stitch);
        unionOverlayLines(live, stitch);
        segs = std::move(live);
        used.assign(segs.size(), 0);
    }

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
    // Exact twins only: same type + matching tips. A digon (plane line + cylinder
    // ellipse sharing both ends) must NOT collapse — that dropped F896's ellipse and
    // left F883's generator as an orphan open (ROBOT z≈525.7255).
    auto sameGeom = [&](const Segment& a, const Segment& b) {
        if (a.type != b.type) return false;
        const bool fwd = dist(segStart(a), segStart(b)) <= stitch && dist(segEnd(a), segEnd(b)) <= stitch;
        const bool rev = dist(segStart(a), segEnd(b)) <= stitch && dist(segEnd(a), segStart(b)) <= stitch;
        if (!fwd && !rev) return false;
        return sameArcFrame(a, b);
    };
    auto canStitch = [&](const Segment& from, const Vec3& joint, const Segment& to) {
        return onBothCircles(joint, from, to);
    };

    // Exact twin drop with A1 ownership: coplanar wins, else lower face_id.
    for (size_t i = 0; i < segs.size(); ++i) {
        if (used[i]) continue;
        for (size_t j = i + 1; j < segs.size(); ++j) {
            if (used[j] || segs[j].closed_loop) continue;
            if (segs[i].solid_id != segs[j].solid_id) continue;
            if (!sameGeom(segs[i].geom, segs[j].geom)) continue;
            if (preferOwner(segs[i], segs[j])) {
                segs[i].coplanar = segs[i].coplanar || segs[j].coplanar;
                used[j] = 1;
            } else {
                segs[j].coplanar = segs[j].coplanar || segs[i].coplanar;
                used[i] = 1;
                break;
            }
        }
    }

    // Phase 3 A1: shared-edge near-twin ownership — coincident edge∩plane with snap drift
    // (endpoints within stitch) from faces that share an edge → keep one owner.
    if (opt.solid_adjacency && !opt.solid_adjacency->empty()) {
        const double endTol = std::max(stitch, designBridge);
        auto nearTwin = [&](const Segment& a, const Segment& b) {
            const bool fwd =
                dist(segStart(a), segStart(b)) <= endTol && dist(segEnd(a), segEnd(b)) <= endTol;
            const bool rev =
                dist(segStart(a), segEnd(b)) <= endTol && dist(segEnd(a), segStart(b)) <= endTol;
            if (!fwd && !rev) return false;
            if (a.type != b.type) return false;
            if (a.type == SegmentType::Arc) return sameArcFrame(a, b);
            if (a.type == SegmentType::Line) {
                // Offset must stay at designBridge — do not collapse thin walls.
                return collinearLines(a, b, designBridge);
            }
            return true;
        };
        for (size_t i = 0; i < segs.size(); ++i) {
            if (used[i] || segs[i].closed_loop) continue;
            for (size_t j = i + 1; j < segs.size(); ++j) {
                if (used[j] || segs[j].closed_loop) continue;
                if (segs[i].solid_id != segs[j].solid_id) continue;
                if (segs[i].face_id < 0 || segs[j].face_id < 0) continue;
                if (segs[i].face_id == segs[j].face_id) continue;
                if (!opt.solid_adjacency->facesShareEdge(segs[i].face_id, segs[j].face_id))
                    continue;
                if (!nearTwin(segs[i].geom, segs[j].geom)) continue;
                if (preferOwner(segs[i], segs[j])) {
                    segs[i].coplanar = segs[i].coplanar || segs[j].coplanar;
                    used[j] = 1;
                } else {
                    segs[j].coplanar = segs[j].coplanar || segs[i].coplanar;
                    used[i] = 1;
                    break;
                }
            }
        }
    }

    for (size_t i = 0; i < segs.size(); ++i) {
        if (used[i]) continue;
        if (segs[i].closed_loop) {
            used[i] = 1;
            contours.push_back(makeContour({stampFace(segs[i])}, segs[i].solid_id, segs[i].shell_id,
                                           segs[i].coplanar, true));
        }
    }

    auto sameBody = [&](size_t a, size_t b) {
        return segs[a].solid_id == segs[b].solid_id;
    };
    auto match = [&](const Vec3& a, const Vec3& b) { return dist(a, b) <= stitch; };
    const SolidAdjacency* adj = opt.solid_adjacency;

    for (size_t seed = 0; seed < segs.size(); ++seed) {
        if (used[seed]) continue;
        used[seed] = 1;
        std::vector<Segment> chain{stampFace(segs[seed])};
        int solid = segs[seed].solid_id;
        int shell = segs[seed].shell_id;
        bool coplanar = segs[seed].coplanar;
        bool grew = true;
        while (grew) {
            grew = false;
            const Vec3 head = segStart(chain.front());
            const Vec3 tail = segEnd(chain.back());
            // Already a closed loop — do not attach dangling stubs past the join,
            // unless unused same-face fragments still touch the joint (huapingdun
            // z≈±11.62: shelf+one rim piece self-closes into a digon ear while the
            // rest of face 4's MST rim fragments remain).
            if (chain.size() >= 2 && match(head, tail)) {
                bool same_face_cont = false;
                const int tip_face = chain.back().face_id;
                if (tip_face >= 0) {
                    for (size_t j = 0; j < segs.size(); ++j) {
                        if (used[j] || !sameBody(seed, j)) continue;
                        if (segs[j].face_id != tip_face) continue;
                        Segment cand = stampFace(segs[j]);
                        if (match(tail, segStart(cand)) || match(tail, segEnd(cand)) ||
                            match(head, segStart(cand)) || match(head, segEnd(cand))) {
                            same_face_cont = true;
                            break;
                        }
                    }
                }
                if (!same_face_cont) break;
            }

            // Prefer candidates on faces that share an edge with the tip face.
            // Same face + same chain_idx = continue one location instance.
            // Same face + different chain_idx = pocket wall↔cap on one face
            // (F1→F2→F3→F2→F1 / wall→cap→wall). Rank below shared-edge neighbors so a
            // real F2 tip still wins, but above "no topo" so corners join
            // (ROBOT_7 z≈603.35: face 225/329 blind-hole U). Do NOT use rank 3 — that
            // left the cap as a digon hole while the outer kept only the opening chord.
            size_t best_j = segs.size();
            int best_mode = -1;  // 0 tail-start, 1 tail-end, 2 head-end, 3 head-start
            int best_rank = 100;
            int best_same = -1;
            auto rankOf = [&](const Segment& from, const Segment& to) {
                if (!adj || adj->empty() || from.face_id < 0 || to.face_id < 0) return 2;
                if (from.face_id == to.face_id) {
                    if (from.chain_idx >= 0 && to.chain_idx >= 0 &&
                        from.chain_idx != to.chain_idx)
                        return 0;  // pocket wall↔cap corner — equal to shared-edge
                    return 0;
                }
                if (adj->facesShareEdge(from.face_id, to.face_id)) return 0;
                if (topoLinkableFaces(from.face_id, to.face_id, *adj)) return 1;
                return 2;
            };
            // Refuse attaching a pocket wall across faces when its free tip faces a
            // same-face partner wall (blind-hole U). Those walls must stay separate
            // until linkOpen/repair forms the U (ROBOT_7 z≈603).
            auto wouldPocketSpur = [&](const Segment& cand, bool jointAtCandStart) {
                const Vec3 free = jointAtCandStart ? segEnd(cand) : segStart(cand);
                for (size_t k = 0; k < segs.size(); ++k) {
                    if (used[k] || !sameBody(seed, k)) continue;
                    Segment o = stampFace(segs[k]);
                    if (o.face_id != cand.face_id || o.face_id < 0) continue;
                    if (cand.chain_idx >= 0 && o.chain_idx >= 0 &&
                        cand.chain_idx == o.chain_idx)
                        continue;
                    if (dist(free, segStart(o)) <= 8.0 || dist(free, segEnd(o)) <= 8.0)
                        return true;
                }
                return false;
            };
            auto take = [&](size_t j, int mode, int rank, int tip_face, const Segment& cand,
                            bool jointAtCandStart) {
                if (cand.face_id != tip_face && wouldPocketSpur(cand, jointAtCandStart))
                    return;
                const int same = (cand.face_id == tip_face) ? 1 : 0;
                if (rank < best_rank || (rank == best_rank && same > best_same)) {
                    best_rank = rank;
                    best_same = same;
                    best_j = j;
                    best_mode = mode;
                }
            };

            for (size_t j = 0; j < segs.size(); ++j) {
                if (used[j] || !sameBody(seed, j)) continue;
                Segment cand = stampFace(segs[j]);
                if (match(tail, segStart(cand)) && canStitch(chain.back(), tail, cand))
                    take(j, 0, rankOf(chain.back(), cand), chain.back().face_id, cand, true);
                if (match(tail, segEnd(cand))) {
                    Segment rev = cand;
                    reverseSegment(rev);
                    if (canStitch(chain.back(), tail, rev))
                        take(j, 1, rankOf(chain.back(), cand), chain.back().face_id, cand, false);
                }
                if (match(head, segEnd(cand)) && canStitch(cand, head, chain.front()))
                    take(j, 2, rankOf(cand, chain.front()), chain.front().face_id, cand, false);
                if (match(head, segStart(cand))) {
                    Segment rev = cand;
                    reverseSegment(rev);
                    if (canStitch(rev, head, chain.front()))
                        take(j, 3, rankOf(cand, chain.front()), chain.front().face_id, cand, true);
                }
            }
            if (best_j >= segs.size()) break;

            Segment cand = stampFace(segs[best_j]);
            if (best_mode == 0) {
                chain.push_back(cand);
            } else if (best_mode == 1) {
                reverseSegment(cand);
                chain.push_back(cand);
            } else if (best_mode == 2) {
                chain.insert(chain.begin(), cand);
            } else {
                reverseSegment(cand);
                chain.insert(chain.begin(), cand);
            }
            used[best_j] = 1;
            coplanar = coplanar || segs[best_j].coplanar;
            grew = true;
        }

        // Trim trailing stubs attached after the chain already returned to its start.
        while (chain.size() >= 3) {
            const Vec3 head = segStart(chain.front());
            if (dist(head, segEnd(chain.back())) <= stitch) break;
            if (dist(head, segEnd(chain[chain.size() - 2])) <= stitch) {
                chain.pop_back();
                continue;
            }
            break;
        }

        weldChainJoints(chain, stitch);

        const bool closed =
            !chain.empty() && dist(segStart(chain.front()), segEnd(chain.back())) <= closeTol;
        if (!closed) {
            const double gap =
                chain.empty() ? 0.0 : dist(segStart(chain.front()), segEnd(chain.back()));
            if (gap > 0.0 && gap <= bridgeMax) {
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
            // Snap head/tail when closing within closeTol but above exact stitch.
            if (!chain.empty()) {
                const Vec3 h = segStart(chain.front());
                const Vec3 t = segEnd(chain.back());
                if (dist(h, t) > stitch) {
                    const Vec3 m{(h.x + t.x) * 0.5, (h.y + t.y) * 0.5, (h.z + t.z) * 0.5};
                    setSegStart(chain.front(), m);
                    setSegEnd(chain.back(), m);
                }
            }
            contours.push_back(makeContour(std::move(chain), solid, shell, coplanar, true));
        }
    }

    // Blind-hole walls + mouth chord: expand notch before tip-merge absorbs walls as spurs.
    repairBlindHoleNotches(contours, stitch, &stats);

    bool merged = true;
    while (merged) {
        merged = false;
        for (size_t i = 0; i < contours.size(); ++i) {
            if (contours[i].segments.empty()) continue;
            for (size_t j = i + 1; j < contours.size(); ++j) {
                if (contours[j].segments.empty()) continue;
                if (tryMergeOpen(contours[i], contours[j], mergeTol, frame)) {
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

    // Drop UVMatch micro-rings before any gap bridge so tips can meet the real body
    // (ROBOT_7 z≈521: four ~0.5 mm closed bsplines blocked a 3.2 mm close).
    dropMicroClosedLoops(contours, 2.0);

    // Blind-hole walls → outer notch BEFORE mid-joint splice steals them as orphans
    // (ROBOT_7 z≈603.35: splice reported 11 mm and left a flat chord).
    repairBlindHoleNotches(contours, stitch, &stats);

    // F1→F1 U-turn pockets left an orphan F2 span (ROBOT_7 z≈544.866).
    // Splice before stripping — splice uses the pocket loops as landmarks.
    spliceOpenOrphansAtMidJoints(contours, stitch, &stats);
    stripShortClosedSubloops(contours, stitch, 25.0);

    // Missing mid-fillet spans — hard-capped with other invent bridges.
    bridgeFilletOpenGaps(contours, frame, kMaxInventBridge, stitch, &stats);

    // Shared-edge / one-hop face adjacency: close gaps left by slick faces with no seeds.
    linkOpenByFaceAdjacency(contours, opt.solid_adjacency, stitch, designSnap, frame, &stats);
    // Same loop, open head↔tail: fillet tip short of neighbor chain (ROBOT_4 z≈167).
    closeOpenSelfByFaceAdjacency(contours, opt.solid_adjacency, stitch, designSnap, &stats);
    repairBlindHoleNotches(contours, stitch, &stats);
    spliceOpenReplacingShortHostEdge(contours, stitch);

    // Sphere∩plane full circle left while vase open tips sit on it (huapingdun+sphere).
    spliceOpenThroughClosedCircles(contours, frame, stitch);
    repairOpenWrongCircleArc(contours, frame, stitch);

    // Micro UVMatch rings first — otherwise dropOpenHuggingClosed treats them as the
    // body and deletes the real opens (ROBOT_7 z≈521).
    dropMicroClosedLoops(contours, 2.0);
    // Tips that were blocked by micro-rings — hard-capped invent bridge.
    bridgeFilletOpenGaps(contours, frame, kMaxInventBridge, stitch, &stats);
    linkOpenByFaceAdjacency(contours, opt.solid_adjacency, stitch, designSnap, frame, &stats);
    closeOpenSelfByFaceAdjacency(contours, opt.solid_adjacency, stitch, designSnap, &stats);
    repairBlindHoleNotches(contours, stitch, &stats);
    spliceOpenReplacingShortHostEdge(contours, stitch);
    // Dual tip-pair leftovers: one aligned stub join leaves a matching self-gap.
    {
        constexpr double maxStub = kMaxInventBridge;  // hard invent-bridge cap
        const double snapGap = designSnap;
        for (Contour& c : contours) {
            if (c.closed || c.segments.size() < 4) continue;
            const Vec3 h = segStart(c.segments.front());
            const Vec3 t = segEnd(c.segments.back());
            const double gap = dist(h, t);
            if (gap <= snapGap) {
                c.closed = true;
                continue;
            }
            if (gap > maxStub) continue;
            const Vec3 th = contourEndTangent(c, frame, true);
            const Vec3 tt = contourEndTangent(c, frame, false);
            if (length(th) < 1e-12 || length(tt) < 1e-12) continue;
            const Vec3 delta = h - t;
            const Vec3 thN = normalized(th);
            const Vec3 ttN = normalized(tt);
            const double latH = length(delta - thN * dot(delta, thN));
            const double latT = length(delta - ttN * dot(delta, ttN));
            if (std::min(latH, latT) > std::max(0.25, 0.1 * gap)) continue;
            Segment bridge;
            bridge.type = SegmentType::Line;
            bridge.start = t;
            bridge.end = h;
            bridge.face_id = c.segments.back().face_id;
            c.segments.push_back(bridge);
            c.closed = true;
            ++stats.bridged;
            stats.max_bridge = std::max(stats.max_bridge, gap);
        }
    }
    // Substantial body open: tip gap only if head/tail faces are adjacent chains.
    {
        const double snapGap = designSnap;
        const double adjGap = kMaxInventBridge;  // hard invent-bridge cap
        const SolidAdjacency* adj = opt.solid_adjacency;
        for (Contour& c : contours) {
            if (c.closed || c.segments.size() < 8) continue;
            if (approxContourLen(c) < 30.0) continue;
            const Vec3 h = segStart(c.segments.front());
            const Vec3 t = segEnd(c.segments.back());
            const double gap = dist(h, t);
            if (gap <= snapGap) {
                c.closed = true;
                continue;
            }
            const int fa = c.segments.front().face_id;
            const int fb = c.segments.back().face_id;
            const int ca = c.segments.front().chain_idx;
            const int cb = c.segments.back().chain_idx;
            if (fa == fb && ca >= 0 && ca == cb) continue;
            // Same-face pocket corner — block diagonal (~11 mm on z≈603).
            if (fa == fb && ca >= 0 && cb >= 0 && ca != cb && gap > 8.0) continue;
            if (!adj || adj->empty() || !topoLinkableFaces(fa, fb, *adj)) continue;
            if (gap > adjGap) continue;
            Segment bridge;
            bridge.type = SegmentType::Line;
            bridge.start = t;
            bridge.end = h;
            bridge.face_id = fb >= 0 ? fb : fa;
            c.segments.push_back(bridge);
            c.closed = true;
            ++stats.bridged;
            stats.max_bridge = std::max(stats.max_bridge, gap);
        }
    }
    // Drop tiny open stubs / fillet extras hugging an already-closed loop.
    // Blind-hole walls must splice into the outer notch first (z≈603) — hug-drop
    // would otherwise delete them as near-boundary opens.
    repairBlindHoleNotches(contours, stitch, &stats);
    spliceOpenReplacingShortHostEdge(contours, stitch);
    dropOpenHuggingClosed(contours, stitch);
    // //cylinder generator orphans past trim wire (ROBOT_4 z≈202).
    dropFloatingOpenOrphans(contours, stitch);
    // F1→F1 pockets closed as false holes on the outer (ROBOT_7 z≈544.866, L≈16 mm).
    dropShortClosedOnLongBoundary(contours, stitch, 20.0);
    // Blind-hole U (walls+cap) replaces the short opening chord on the outer.
    spliceOpenReplacingShortHostEdge(contours, stitch);
    // Blind-hole cap digons (bspline + reverse line) after failed F1 chain join.
    dropDigonClosedContours(contours, stitch);
    // Digon ears on the closed outer (ROBOT_7 L0037/L0041: segment + reverse).
    stripEmbeddedDigonEars(contours, stitch);
    dropPocketDiagonalChords(contours, stitch);
    // Digon pocket-cap chord (= tip gap): drop chord, keep notch (ROBOT_7 z≈626.747).
    dropOpenDigonCapChords(contours, stitch);
    // Thin lids / short false holes that appear after pocket repairs.
    dropShortClosedOnLongBoundary(contours, stitch, 20.0);
    dropDigonClosedContours(contours, stitch);

    // Blind-hole walls that fitted as bowed cubics: collapse near-axis-aligned
    // BSplines to Lines (ROBOT_7 z≈603: one wall max_dev≈3 mm on a 10 mm span).
    for (Contour& c : contours) {
        for (Segment& s : c.segments) {
            if (s.type != SegmentType::BSpline || s.ctrl_pts.size() < 2) continue;
            const Vec3 a = segStart(s), b = segEnd(s);
            const double dx = std::abs(a.x - b.x), dy = std::abs(a.y - b.y);
            if (!(dx <= 0.15 || dy <= 0.15)) continue;
            const double chord = dist(a, b);
            if (chord < 2.0 || chord > 30.0) continue;
            const Vec3 ab = b - a;
            const double L2 = length2(ab);
            double maxd = 0.0;
            for (const Vec3& p : s.ctrl_pts) {
                double t = (L2 > 0) ? dot(p - a, ab) / L2 : 0.0;
                t = std::max(0.0, std::min(1.0, t));
                maxd = std::max(maxd, dist(p, a + ab * t));
            }
            if (maxd <= 0.5) continue;
            s.type = SegmentType::Line;
            s.ctrl_pts.clear();
            s.knots.clear();
            s.weights.clear();
            s.start = a;
            s.end = b;
        }
    }

    for (Contour& c : contours) weldChainJoints(c.segments, stitch);

    for (Contour& c : contours) coalesceWrapSplitArcs(c, stitch);

    // Phase 2 C1: drop closed pockets that only duplicate open endpoints.
    const size_t beforeFix = contours.size();
    fixContours(contours, std::max(opt.tolerance, 1e-3));
    stats.spurious_closed_removed +=
        static_cast<int>(beforeFix) - static_cast<int>(contours.size());

    // Nesting + orientation per solid (multi-solid STEP: each body has its own outer/inner tree).
    std::map<int, std::vector<Contour>> bySolid;
    for (Contour& c : contours) bySolid[c.solid_id].push_back(std::move(c));
    contours.clear();
    for (auto& kv : bySolid) {
        nestAndOrientContours(kv.second, frame);
        // Phase 2 C2: drop near-duplicate closed loops after nest.
        stats.spurious_closed_removed +=
            dropDuplicateClosedLoops(kv.second, frame, std::max(stitch, 1e-3));
        nestAndOrientContours(kv.second, frame);  // refresh parents after drops
        const int base = static_cast<int>(contours.size());
        for (Contour& c : kv.second) {
            if (c.parent) *c.parent += base;
            contours.push_back(std::move(c));
        }
    }

    // Recount opens after cleanup merges/drops.
    stats.open_leftover = 0;
    for (const Contour& c : contours) {
        if (!c.closed) ++stats.open_leftover;
    }

    layer.contours = std::move(contours);
    return layer;
}

void reassembleCrossSolidContours(std::vector<Contour>& contours, const SliceFrame& frame,
                                  const SliceOptions& opt, AssembleStats& stats) {
    if (contours.empty()) return;
    const double tol = std::max(opt.tolerance, opt.stitch_tolerance);

    // Cross-solid overlay: near-identical closed loops from different solids → keep one.
    stats.spurious_closed_removed += dropDuplicateClosedLoops(contours, frame, std::max(tol, 1e-3));

    // Cross-solid / same-layer fillet mid-gaps left as two open halves.
    bridgeFilletOpenGaps(contours, frame, 3.5, tol, &stats);

    spliceOpenThroughClosedCircles(contours, frame, tol);
    repairOpenWrongCircleArc(contours, frame, tol);

    dropMicroClosedLoops(contours, 2.0);
    bridgeFilletOpenGaps(contours, frame, 3.5, tol, &stats);
    // Fillet / mate extras from another solid hugging a closed body contour.
    dropOpenHuggingClosed(contours, tol);
    dropFloatingOpenOrphans(contours, tol);

    for (Contour& c : contours) c.parent = std::nullopt;

    // Global nest/orient across all solids (containment / hole relationships).
    nestAndOrientContours(contours, frame);

    stats.open_leftover = 0;
    for (const Contour& c : contours) {
        if (!c.closed) ++stats.open_leftover;
    }
}

}  // namespace brepslicer
