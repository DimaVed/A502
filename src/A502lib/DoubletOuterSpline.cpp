// DoubletOuterSpline.cpp
// ----------------------
// Doublet analysis outer spline for deep interior panels.
//
// Port of the interior branches of blccal()/blcal():
//   * five-by-five grid stencil and tangent-plane frame (msrotm);
//   * 4x4 block of panel-centre doublet parameters;
//   * weighted quadratic least-squares fit for panel corners (16 points)
//     and for the two adjacent edge midpoints (12 points);
//   * assembly of the 9x25 B_D matrix.

#include "OuterSplines.h"

#include "AeroModel.h"
#include "PanelGeometry.h"

#include <algorithm>
#include <array>
#include <cmath>

namespace a502 {

namespace {

constexpr double kDeltf = 1000.0;

struct TangentFrame {
    Eigen::Matrix3d toLocal{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d toGlobal{Eigen::Matrix3d::Identity()};
    Eigen::Vector3d normal{Eigen::Vector3d::UnitZ()};
    bool valid{false};
};

Eigen::Vector3d networkPoint(const Network& net, int m, int n) {
    // m,n are 1-based.
    const std::size_t g =
        static_cast<std::size_t>(n - 1) * static_cast<std::size_t>(net.nm) +
        static_cast<std::size_t>(m - 1);
    return Eigen::Vector3d(net.coordinates[g * 3],
                           net.coordinates[g * 3 + 1],
                           net.coordinates[g * 3 + 2]);
}

bool buildStencil5(const Network& net, int m, int n,
                   std::array<Eigen::Vector3d, 25>& out) {
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 5; ++i) {
            const int mm = m + i - 2;
            const int nn = n + j - 2;
            if (mm < 1 || mm > net.nm || nn < 1 || nn > net.nn)
                return false;
            out[static_cast<std::size_t>(i + 5 * j)] =
                networkPoint(net, mm, nn);
        }
    }
    return true;
}

bool buildParameterBlock4(const Network& net, int m, int n,
                          std::array<Eigen::Vector3d, 16>& points,
                          std::array<int, 16>& panels) {
    const int nPanelM = net.nm - 1;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            const int pi = m + i - 2; // 1-based panel row
            const int pj = n + j - 2; // 1-based panel column
            if (pi < 1 || pi > nPanelM || pj < 1 || pj > net.nn - 1)
                return false;
            const PanelGeometry g = buildPanelGeometry(
                net.coordinates, net.nm, net.nn, pi - 1, pj - 1);
            if (g.degenerate)
                return false;
            const std::size_t k = static_cast<std::size_t>(i + 4 * j);
            points[k] = g.cp[8];
            panels[k] = (pi - 1) + nPanelM * (pj - 1);
        }
    }
    return true;
}

TangentFrame makeTangentFrame(const Eigen::Vector3d* z, int nm, int nn,
                              int m, int n) {
    TangentFrame out;
    auto at = [&](int i, int j) -> Eigen::Vector3d {
        return z[(i - 1) + nm * (j - 1)];
    };

    const int mm1 = std::max(1, m - 1);
    const int mp1 = std::min(nm, m + 1);
    const int nm1 = std::max(1, n - 1);
    const int np1 = std::min(nn, n + 1);

    Eigen::Vector3d u =
        0.125 * (at(mm1, np1) + 2.0 * at(m, np1) + at(mp1, np1) -
                 at(mm1, nm1) - 2.0 * at(m, nm1) - at(mp1, nm1));
    Eigen::Vector3d v =
        0.125 * (at(mp1, nm1) + 2.0 * at(mp1, n) + at(mp1, np1) -
                 at(mm1, nm1) - 2.0 * at(mm1, n) - at(mm1, np1));
    Eigen::Vector3d w = u.cross(v);
    if (w.norm() <= 1.0e-14)
        return out;
    v = w.cross(u);
    if (u.norm() <= 1.0e-14 || v.norm() <= 1.0e-14)
        return out;

    u.normalize();
    v.normalize();
    w.normalize();

    out.toGlobal.col(0) = u;
    out.toGlobal.col(1) = v;
    out.toGlobal.col(2) = w;
    out.toLocal = out.toGlobal.transpose();
    out.normal = w;
    out.valid = true;
    return out;
}

Eigen::Vector2d projectLocal(const TangentFrame& frame,
                             const Eigen::Vector3d& origin,
                             const Eigen::Vector3d& p) {
    Eigen::Vector3d d = p - origin;
    const double len = d.norm();
    if (len <= 1.0e-15)
        return Eigen::Vector2d::Zero();
    const double n2 = frame.normal.squaredNorm();
    if (n2 > 1.0e-30)
        d -= frame.normal * (d.dot(frame.normal) / n2);
    if (d.norm() <= 1.0e-15)
        return Eigen::Vector2d::Zero();
    d.normalize();
    const Eigen::Vector3d projected = origin + len * d;
    const Eigen::Vector3d local = frame.toLocal * (projected - origin);
    return Eigen::Vector2d(local.x(), local.y());
}

double lsqWeight(const Compressibility& flow, const Eigen::Vector3d& v) {
    const double amachm = 0.5 * (1.0 - flow.sbetam) * flow.mach;
    const double vn = v.norm();
    if (vn <= 1.0e-15)
        return 1.0;
    return 1.0 + amachm * (1.0 - v.dot(flow.compd) / vn);
}

// First row of the weighted quadratic least-squares inverse: the value of
// the fitted function at the local origin.
Eigen::VectorXd quadraticConstantRow(
    const std::array<Eigen::Vector2d, 16>& xy,
    const std::array<double, 16>& w, int n) {
    Eigen::MatrixXd A(n, 6);
    Eigen::VectorXd wv(n);
    for (int k = 0; k < n; ++k) {
        const double x = xy[static_cast<std::size_t>(k)].x();
        const double y = xy[static_cast<std::size_t>(k)].y();
        A(k, 0) = 1.0;
        A(k, 1) = x;
        A(k, 2) = y;
        A(k, 3) = 0.5 * x * x;
        A(k, 4) = x * y;
        A(k, 5) = 0.5 * y * y;
        wv(k) = w[static_cast<std::size_t>(k)];
    }
    const Eigen::MatrixXd W2 = wv.array().square().matrix().asDiagonal();
    const Eigen::Matrix<double, 6, 6> M = A.transpose() * W2 * A;
    if (std::abs(M.determinant()) <= 1.0e-30)
        return Eigen::VectorXd::Zero(n);
    const Eigen::Matrix<double, 6, Eigen::Dynamic> pinv =
        M.inverse() * A.transpose() * W2;
    return pinv.row(0).transpose();
}

bool cornerRow(const Network& net, int m, int n,
               const Compressibility& flow, std::array<double, 16>& b) {
    std::array<Eigen::Vector3d, 25> stencil;
    std::array<Eigen::Vector3d, 16> params;
    std::array<int, 16> panels;
    if (!buildStencil5(net, m, n, stencil) ||
        !buildParameterBlock4(net, m, n, params, panels))
        return false;

    const TangentFrame frame = makeTangentFrame(stencil.data(), 5, 5, 3, 3);
    if (!frame.valid)
        return false;
    const Eigen::Vector3d origin = stencil[12];

    std::array<Eigen::Vector2d, 16> xy;
    std::array<double, 16> w;
    for (int k = 0; k < 16; ++k) {
        const std::size_t kk = static_cast<std::size_t>(k);
        xy[kk] = projectLocal(frame, origin, params[kk]);
        w[kk] = lsqWeight(flow, params[kk] - origin);
        if (k == 5 || k == 6 || k == 9 || k == 10)
            w[kk] = kDeltf;
        if (k == 0 || k == 3 || k == 12 || k == 15)
            w[kk] = 0.0;
    }

    const Eigen::VectorXd row = quadraticConstantRow(xy, w, 16);
    for (int k = 0; k < 16; ++k)
        b[static_cast<std::size_t>(k)] = row(k);
    return true;
}

bool midColRow(const Network& net, int m, int n,
               const Compressibility& flow, std::array<double, 12>& b) {
    std::array<Eigen::Vector3d, 25> stencil;
    std::array<Eigen::Vector3d, 16> params;
    std::array<int, 16> panels;
    if (!buildStencil5(net, m, n, stencil) ||
        !buildParameterBlock4(net, m, n, params, panels))
        return false;

    std::array<Eigen::Vector3d, 9> ze;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            const Eigen::Vector3d a =
                stencil[static_cast<std::size_t>((i + 1) + 5 * 2)];
            const Eigen::Vector3d c =
                stencil[static_cast<std::size_t>((i + 1) + 5 * 3)];
            ze[static_cast<std::size_t>(i + 3 * j)] =
                a + 0.5 * static_cast<double>(j) * (c - a);
        }
    }

    const TangentFrame frame = makeTangentFrame(ze.data(), 3, 3, 2, 2);
    if (!frame.valid)
        return false;
    const Eigen::Vector3d origin = ze[4];

    std::array<Eigen::Vector2d, 16> xy;
    std::array<double, 16> w;
    int k = 0;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 4; ++i) {
            const Eigen::Vector3d p =
                params[static_cast<std::size_t>(i + 4 * (j + 1))];
            xy[static_cast<std::size_t>(k)] = projectLocal(frame, origin, p);
            w[static_cast<std::size_t>(k)] = lsqWeight(flow, p - origin);
            if (j == 1 && (i == 1 || i == 2))
                w[static_cast<std::size_t>(k)] = kDeltf;
            ++k;
        }
    }

    const Eigen::VectorXd row = quadraticConstantRow(xy, w, 12);
    for (int q = 0; q < 12; ++q)
        b[static_cast<std::size_t>(q)] = row(q);
    return true;
}

bool midRowRow(const Network& net, int m, int n,
               const Compressibility& flow, std::array<double, 12>& b) {
    std::array<Eigen::Vector3d, 25> stencil;
    std::array<Eigen::Vector3d, 16> params;
    std::array<int, 16> panels;
    if (!buildStencil5(net, m, n, stencil) ||
        !buildParameterBlock4(net, m, n, params, panels))
        return false;

    std::array<Eigen::Vector3d, 9> ze;
    for (int j = 0; j < 3; ++j) {
        for (int i = 0; i < 3; ++i) {
            const Eigen::Vector3d a =
                stencil[static_cast<std::size_t>(2 + 5 * (j + 1))];
            const Eigen::Vector3d c =
                stencil[static_cast<std::size_t>(3 + 5 * (j + 1))];
            ze[static_cast<std::size_t>(i + 3 * j)] =
                a + 0.5 * static_cast<double>(i) * (c - a);
        }
    }

    const TangentFrame frame = makeTangentFrame(ze.data(), 3, 3, 2, 2);
    if (!frame.valid)
        return false;
    const Eigen::Vector3d origin = ze[4];

    std::array<Eigen::Vector2d, 16> xy;
    std::array<double, 16> w;
    int k = 0;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 3; ++i) {
            const Eigen::Vector3d p =
                params[static_cast<std::size_t>((i + 1) + 4 * j)];
            xy[static_cast<std::size_t>(k)] = projectLocal(frame, origin, p);
            w[static_cast<std::size_t>(k)] = lsqWeight(flow, p - origin);
            if (i == 1 && (j == 1 || j == 2))
                w[static_cast<std::size_t>(k)] = kDeltf;
            ++k;
        }
    }

    const Eigen::VectorXd row = quadraticConstantRow(xy, w, 12);
    for (int q = 0; q < 12; ++q)
        b[static_cast<std::size_t>(q)] = row(q);
    return true;
}

} // namespace

DoubletOuterSpline buildDoubletOuterSpline(const Network& net,
                                           const Compressibility& flow,
                                           int panelI, int panelJ) {
    DoubletOuterSpline out;
    const int nm = net.nm;
    const int nn = net.nn;
    const int nPanelM = nm - 1;
    const int nPanelN = nn - 1;
    if (nPanelM < 5 || nPanelN < 5 || panelI < 0 || panelJ < 0 ||
        panelI >= nPanelM || panelJ >= nPanelN)
        return out;

    // Deep interior: the 5x5 block of panel-centre parameters exists.
    const int pRow = panelI + 1; // 1-based
    const int pCol = panelJ + 1;
    if (pRow < 3 || pRow > nPanelM - 2 || pCol < 3 ||
        pCol > nPanelN - 2)
        return out;

    const std::size_t nPoints =
        static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
    if (net.coordinates.size() < nPoints * 3)
        return out;

    for (int k = 0; k < 25; ++k)
        out.neighborPanel[static_cast<std::size_t>(k)] = -1;

    // Panel-centre indices of the 5x5 block.
    for (int j = 0; j < 5; ++j) {
        for (int i = 0; i < 5; ++i) {
            const int pr = pRow - 2 + i;
            const int pc = pCol - 2 + j;
            out.neighborPanel[static_cast<std::size_t>(i + 5 * j)] =
                (pr - 1) + nPanelM * (pc - 1);
        }
    }

    // Panel-centre value.
    out.coefficients(8, 12) = 1.0;

    // Four corner rows.
    std::array<std::array<double, 16>, 4> corner{};
    for (int lcp = 1; lcp <= 4; ++lcp) {
        const int mc = ((lcp - 1) % 2) + 1;
        const int nc = ((lcp - 1) / 2) + 1;
        if (!cornerRow(net, pRow + mc - 1, pCol + nc - 1, flow,
                       corner[static_cast<std::size_t>(lcp - 1)]))
            return DoubletOuterSpline{};

        const int lcpp = 2 * mc - 1 + std::abs(mc - nc);
        for (int i = 1; i <= 16; ++i) {
            const int mi = ((i - 1) % 4) + 1;
            const int nj = ((i - 1) / 4) + 1;
            const int col = (mc + mi - 1) + 5 * (nc + nj - 2); // 1-based
            out.coefficients(lcpp - 1, col - 1) =
                corner[static_cast<std::size_t>(lcp - 1)]
                      [static_cast<std::size_t>(i - 1)];
        }
    }

    // Four edge-midpoint rows.
    std::array<double, 12> colMid1{}, colMid2{}, rowMid1{}, rowMid3{};
    if (!midColRow(net, pRow, pCol, flow, colMid1))
        return DoubletOuterSpline{};
    if (!midColRow(net, pRow + 1, pCol, flow, colMid2))
        return DoubletOuterSpline{};
    if (!midRowRow(net, pRow, pCol, flow, rowMid1))
        return DoubletOuterSpline{};
    if (!midRowRow(net, pRow, pCol + 1, flow, rowMid3))
        return DoubletOuterSpline{};

    for (int i = 1; i <= 12; ++i) {
        const int m4 = ((i - 1) % 4) + 1;
        const int n4 = ((i - 1) / 4) + 1;
        const int m3 = ((i - 1) % 3) + 1;
        const int n3 = ((i - 1) / 3) + 1;

        const int k5 = m4 + 5 * n4;
        const int k6 = m3 + 5 * n3 + 1;
        const int k7 = m4 + 5 * n4 + 1;
        const int k8 = m3 + 5 * n3 - 4;

        out.coefficients(4, k5 - 1) = colMid1[static_cast<std::size_t>(i - 1)];
        out.coefficients(5, k6 - 1) = rowMid3[static_cast<std::size_t>(i - 1)];
        out.coefficients(6, k7 - 1) = colMid2[static_cast<std::size_t>(i - 1)];
        out.coefficients(7, k8 - 1) = rowMid1[static_cast<std::size_t>(i - 1)];
    }

    out.valid = true;
    return out;
}

} // namespace a502
