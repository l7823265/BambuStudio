#include <intersect/BSplineFit.h>
#include <geom/GeomUtil.h>

#include <algorithm>
#include <cmath>
#include <vector>

namespace brepslicer {
namespace {

int findSpan(int nctrl, int p, double t, const std::vector<double>& knots) {
    if (t >= knots[nctrl]) return nctrl - 1;
    if (t <= knots[p]) return p;
    int low = p;
    int high = nctrl;
    int mid = (low + high) / 2;
    while (t < knots[mid] || t >= knots[mid + 1]) {
        if (t < knots[mid]) high = mid;
        else low = mid;
        mid = (low + high) / 2;
    }
    return mid;
}

bool solveLinear(std::vector<std::vector<double>>& A, std::vector<double>& b);

void basisFuns(int span, double t, int p, const std::vector<double>& U, std::vector<double>& N) {
    N.assign(p + 1, 0.0);
    std::vector<double> left(p + 1, 0.0), right(p + 1, 0.0);
    N[0] = 1.0;
    for (int j = 1; j <= p; ++j) {
        left[j] = t - U[span + 1 - j];
        right[j] = U[span + j] - t;
        double saved = 0.0;
        for (int r = 0; r < j; ++r) {
            const double den = right[r + 1] + left[j - r];
            const double tmp = (den != 0.0) ? N[r] / den : 0.0;
            N[r] = saved + right[r + 1] * tmp;
            saved = left[j - r] * tmp;
        }
        N[j] = saved;
    }
}

// Piegl/Tiller averaging knots for interpolating a cubic at chord parameters t.
std::vector<double> interpolatingKnots(const std::vector<double>& t, int p) {
    const int n = static_cast<int>(t.size()) - 1;
    const int m = n + p + 1;
    std::vector<double> U(m + 1, 0.0);
    for (int i = 0; i <= p; ++i) U[i] = 0;
    for (int i = m - p; i <= m; ++i) U[i] = 1;
    for (int j = 1; j <= n - p; ++j) {
        double s = 0;
        for (int i = j; i <= j + p - 1; ++i) s += t[i];
        U[j + p] = s / p;
    }
    return U;
}

std::vector<Vec3> interpolate(const std::vector<Vec3>& pts, const std::vector<double>& t, int p,
                              const std::vector<double>& knots) {
    const int nctrl = static_cast<int>(pts.size());
    if (nctrl < p + 1) return {};
    std::vector<std::vector<double>> A(nctrl, std::vector<double>(nctrl, 0.0));
    for (int j = 0; j < nctrl; ++j) {
        const int span = findSpan(nctrl, p, t[j], knots);
        std::vector<double> N;
        basisFuns(span, t[j], p, knots, N);
        for (int r = 0; r <= p; ++r) {
            const int i = span - p + r;
            if (i >= 0 && i < nctrl) A[j][i] = N[r];
        }
    }
    std::vector<Vec3> ctrl(nctrl);
    for (int c = 0; c < 3; ++c) {
        auto M = A;
        std::vector<double> b(nctrl);
        for (int i = 0; i < nctrl; ++i) {
            b[i] = (c == 0) ? pts[i].x : (c == 1) ? pts[i].y : pts[i].z;
        }
        if (!solveLinear(M, b)) return {};
        for (int i = 0; i < nctrl; ++i) {
            if (c == 0) ctrl[i].x = b[i];
            else if (c == 1) ctrl[i].y = b[i];
            else ctrl[i].z = b[i];
        }
    }
    return ctrl;
}

std::vector<Vec3> subsample(const std::vector<Vec3>& pts, int nkeep) {
    if (nkeep >= static_cast<int>(pts.size())) return pts;
    nkeep = std::max(nkeep, 4);
    std::vector<Vec3> out;
    out.reserve(nkeep);
    const int last = static_cast<int>(pts.size()) - 1;
    for (int i = 0; i < nkeep; ++i) {
        const int idx = (i == nkeep - 1) ? last : (i * last) / (nkeep - 1);
        out.push_back(pts[idx]);
    }
    return out;
}

Vec3 deBoor(const std::vector<Vec3>& ctrl, int p, double t, const std::vector<double>& knots) {
    const int n = static_cast<int>(ctrl.size());
    t = std::min(std::max(t, knots[p]), knots[n]);
    const int span = findSpan(n, p, t, knots);
    std::vector<Vec3> d(p + 1);
    for (int j = 0; j <= p; ++j) d[j] = ctrl[span - p + j];
    for (int r = 1; r <= p; ++r) {
        for (int j = p; j >= r; --j) {
            const int i = span - p + j;
            const double den = knots[i + p - r + 1] - knots[i];
            const double alpha = (den > 0) ? (t - knots[i]) / den : 0;
            d[j] = d[j - 1] * (1.0 - alpha) + d[j] * alpha;
        }
    }
    return d[p];
}

std::vector<double> chordParams(const std::vector<Vec3>& pts) {
    std::vector<double> t(pts.size(), 0.0);
    for (size_t i = 1; i < pts.size(); ++i) t[i] = t[i - 1] + dist(pts[i - 1], pts[i]);
    const double L = t.back();
    if (L <= 1e-18) return t;
    for (double& v : t) v /= L;
    t.back() = 1.0;
    // Collapse zero-length runs so the collocation matrix stays invertible.
    for (size_t i = 1; i < t.size(); ++i) {
        if (t[i] <= t[i - 1]) t[i] = std::min(1.0, t[i - 1] + 1e-12);
    }
    t.back() = 1.0;
    return t;
}

bool solveLinear(std::vector<std::vector<double>>& A, std::vector<double>& b) {
    const int n = static_cast<int>(b.size());
    for (int k = 0; k < n; ++k) {
        int piv = k;
        for (int i = k + 1; i < n; ++i) {
            if (std::abs(A[i][k]) > std::abs(A[piv][k])) piv = i;
        }
        std::swap(A[k], A[piv]);
        std::swap(b[k], b[piv]);
        if (std::abs(A[k][k]) < 1e-18) return false;
        for (int i = k + 1; i < n; ++i) {
            const double f = A[i][k] / A[k][k];
            b[i] -= f * b[k];
            for (int j = k; j < n; ++j) A[i][j] -= f * A[k][j];
        }
    }
    for (int i = n - 1; i >= 0; --i) {
        double s = b[i];
        for (int j = i + 1; j < n; ++j) s -= A[i][j] * b[j];
        b[i] = s / A[i][i];
    }
    return true;
}

bool isPiecewiseBezier(const Segment& s) {
    if (s.degree != 3) return false;
    const int nctrl = static_cast<int>(s.ctrl_pts.size());
    if (nctrl < 4 || (nctrl - 1) % 3 != 0) return false;
    const int nseg = (nctrl - 1) / 3;
    if (static_cast<int>(s.knots.size()) != nctrl + 4) return false;
    if (nseg == 1) return true;
    // Interior knots of a cubic Bezier spline are triple.
    return std::abs(s.knots[4] - s.knots[5]) < 1e-14 && std::abs(s.knots[5] - s.knots[6]) < 1e-14;
}

Vec3 evalBezierSeg(const std::vector<Vec3>& ctrl, int i, double u) {
    const Vec3 b0 = ctrl[3 * i];
    const Vec3 b1 = ctrl[3 * i + 1];
    const Vec3 b2 = ctrl[3 * i + 2];
    const Vec3 b3 = ctrl[3 * i + 3];
    const Vec3 c0 = b0 * (1 - u) + b1 * u;
    const Vec3 c1 = b1 * (1 - u) + b2 * u;
    const Vec3 c2 = b2 * (1 - u) + b3 * u;
    const Vec3 d0 = c0 * (1 - u) + c1 * u;
    const Vec3 d1 = c1 * (1 - u) + c2 * u;
    return d0 * (1 - u) + d1 * u;
}

int bezierSegCount(const Segment& s) {
    return (static_cast<int>(s.ctrl_pts.size()) - 1) / 3;
}

int evalSampleCount(const Segment& s) {
    const int nctrl = static_cast<int>(s.ctrl_pts.size());
    int nseg = 1;
    if (isPiecewiseBezier(s)) nseg = std::max(1, bezierSegCount(s));
    else nseg = std::max(1, nctrl - 3);
    return std::max(256, nseg * 8);
}

Segment finishCubic(std::vector<Vec3> ctrl, std::vector<double> knots) {
    Segment s;
    s.type = SegmentType::BSpline;
    s.degree = 3;
    s.ctrl_pts = std::move(ctrl);
    s.knots = std::move(knots);
    s.weights.assign(s.ctrl_pts.size(), 1.0);
    if (!s.ctrl_pts.empty()) {
        s.start = s.ctrl_pts.front();
        s.end = s.ctrl_pts.back();
    }
    return s;
}

// Encode the polyline as C0 cubic Bezier (linear cubics). Interpolates every vertex.
Segment chordBezier(const std::vector<Vec3>& pts) {
    Segment empty;
    empty.type = SegmentType::BSpline;
    empty.degree = 3;
    if (pts.size() < 2) return empty;
    const int nseg = static_cast<int>(pts.size()) - 1;
    std::vector<Vec3> ctrl;
    ctrl.reserve(3 * nseg + 1);
    ctrl.push_back(pts[0]);
    for (int i = 0; i < nseg; ++i) {
        const Vec3 d = pts[i + 1] - pts[i];
        ctrl.push_back(pts[i] + d * (1.0 / 3.0));
        ctrl.push_back(pts[i] + d * (2.0 / 3.0));
        ctrl.push_back(pts[i + 1]);
    }
    std::vector<double> knots;
    knots.insert(knots.end(), 4, 0.0);
    for (int i = 1; i < nseg; ++i) knots.insert(knots.end(), 3, static_cast<double>(i));
    knots.insert(knots.end(), 4, static_cast<double>(nseg));
    return finishCubic(std::move(ctrl), std::move(knots));
}

Segment interpolateCubic(const std::vector<Vec3>& pts) {
    Segment empty;
    empty.type = SegmentType::BSpline;
    empty.degree = 3;
    if (pts.size() < 4) return empty;
    const std::vector<double> t = chordParams(pts);
    const std::vector<double> knots = interpolatingKnots(t, 3);
    const std::vector<Vec3> ctrl = interpolate(pts, t, 3, knots);
    if (ctrl.size() < 4) return empty;
    return finishCubic(ctrl, knots);
}

double pointToSpline(const Segment& s, const Vec3& p, int nseg_samples) {
    double best = 1e300;
    if (isPiecewiseBezier(s)) {
        const int nseg = bezierSegCount(s);
        const int ns = std::max(8, nseg_samples);
        for (int i = 0; i < nseg; ++i) {
            for (int k = 0; k <= ns; ++k) {
                const double u = static_cast<double>(k) / ns;
                best = std::min(best, dist(p, evalBezierSeg(s.ctrl_pts, i, u)));
            }
        }
        return best;
    }
    const int n = std::max(64, nseg_samples * 8);
    for (int i = 0; i <= n; ++i) {
        const double t = static_cast<double>(i) / n;
        best = std::min(best, dist(p, bsplineEval(s, t)));
    }
    return best;
}

}  // namespace

Vec3 bsplineEval(const Segment& s, double t01) {
    if (s.ctrl_pts.empty()) return s.start;
    t01 = std::max(0.0, std::min(1.0, t01));
    if (isPiecewiseBezier(s)) {
        const int nseg = bezierSegCount(s);
        double T = t01 * nseg;
        int i = static_cast<int>(T);
        if (i >= nseg) i = nseg - 1;
        return evalBezierSeg(s.ctrl_pts, i, T - i);
    }
    const int p = std::max(1, s.degree);
    const auto& knots = s.knots;
    const int nctrl = static_cast<int>(s.ctrl_pts.size());
    if (static_cast<int>(knots.size()) < nctrl + p + 1) {
        return s.start * (1.0 - t01) + s.end * t01;
    }
    const double t = knots.front() + t01 * (knots.back() - knots.front());
    return deBoor(s.ctrl_pts, p, t, knots);
}

void bsplineSample(const Segment& s, int n, std::vector<Vec3>& out) {
    n = std::max(n, 2);
    out.clear();
    out.reserve(n);
    for (int i = 0; i < n; ++i) {
        const double t = static_cast<double>(i) / (n - 1);
        out.push_back(bsplineEval(s, t));
    }
}

double bsplineLength(const Segment& s) {
    std::vector<Vec3> pts;
    bsplineSample(s, evalSampleCount(s), pts);
    double L = 0;
    for (size_t i = 1; i < pts.size(); ++i) L += dist(pts[i - 1], pts[i]);
    return L;
}

double bsplineHausdorff(const Segment& s, const std::vector<Vec3>& poly) {
    if (poly.size() < 2 || s.ctrl_pts.empty()) return 0;
    const int ns = isPiecewiseBezier(s) ? 12 : 32;
    double h = 0;
    for (const Vec3& p : poly) h = std::max(h, pointToSpline(s, p, ns));
    // Spline → polyline (sampled).
    std::vector<Vec3> spl;
    bsplineSample(s, evalSampleCount(s), spl);
    for (const Vec3& q : spl) {
        double best = 1e300;
        for (const Vec3& p : poly) best = std::min(best, dist(p, q));
        // Also distance to polyline segments.
        for (size_t i = 1; i < poly.size(); ++i) {
            const Vec3 a = poly[i - 1];
            const Vec3 ab = poly[i] - a;
            const double L2 = length2(ab);
            double t = (L2 > 0) ? dot(q - a, ab) / L2 : 0;
            t = std::max(0.0, std::min(1.0, t));
            best = std::min(best, dist(q, a + ab * t));
        }
        h = std::max(h, best);
    }
    return h;
}

Segment fitCubicBSpline(const std::vector<Vec3>& pts_in, bool closed, double max_err) {
    std::vector<Vec3> pts = pts_in;
    if (pts.size() < 2) {
        Segment s;
        s.type = SegmentType::BSpline;
        s.degree = 3;
        return s;
    }
    if (closed && dist(pts.front(), pts.back()) > 1e-12) pts.push_back(pts.front());

    Segment best;
    best.type = SegmentType::BSpline;
    best.degree = 3;
    double bestErr = 1e300;
    int bestCtrl = 1e9;

    auto consider = [&](Segment s) {
        if (s.ctrl_pts.size() < 4) return;
        const double err = bsplineHausdorff(s, pts);
        s.fit_error = err;
        const int nc = static_cast<int>(s.ctrl_pts.size());
        const bool ok = err <= max_err;
        const bool bestOk = bestErr <= max_err;
        if (ok && (!bestOk || nc < bestCtrl || (nc == bestCtrl && err < bestErr))) {
            best = s;
            bestErr = err;
            bestCtrl = nc;
        } else if (!ok && !bestOk && err < bestErr) {
            best = s;
            bestErr = err;
            bestCtrl = nc;
        }
    };

    const int nmax = std::max(4, static_cast<int>(pts.size()));
    const int ninterp = std::min(nmax, 48);
    for (int nkeep = 4; nkeep <= ninterp; nkeep = (nkeep < 32) ? nkeep + 4 : nkeep * 2) {
        const int nk = std::min(nkeep, ninterp);
        const std::vector<Vec3> samp = subsample(pts, nk);
        consider(interpolateCubic(samp));
        if (bestErr <= max_err) break;
        if (nk == ninterp) break;
    }
    // Polyline-as-cubic always interpolates the march vertices (true-curve samples).
    consider(chordBezier(pts));
    if (best.ctrl_pts.size() < 4) best = chordBezier(pts);
    best.fit_error = (bestErr < 1e299) ? bestErr : bsplineHausdorff(best, pts);
    return best;
}

Segment cubicFromPolyline(const std::vector<Vec3>& pts_in, bool closed) {
    std::vector<Vec3> pts = pts_in;
    if (pts.size() < 2) {
        Segment s;
        s.type = SegmentType::BSpline;
        s.degree = 3;
        return s;
    }
    if (closed && dist(pts.front(), pts.back()) > 1e-12) pts.push_back(pts.front());
    Segment s = chordBezier(pts);
    s.fit_error = 0;
    return s;
}

void reverseBSpline(Segment& s) {
    std::reverse(s.ctrl_pts.begin(), s.ctrl_pts.end());
    std::reverse(s.weights.begin(), s.weights.end());
    if (s.knots.size() >= 2) {
        const double a = s.knots.front();
        const double b = s.knots.back();
        std::reverse(s.knots.begin(), s.knots.end());
        for (double& k : s.knots) k = a + b - k;
    }
    std::swap(s.start, s.end);
}

}  // namespace brepslicer
