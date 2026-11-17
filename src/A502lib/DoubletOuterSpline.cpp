// DoubletOuterSpline.cpp
// ----------------------
// Doublet analysis outer splines: ports of blccal()/blcal().
//
//   * five-by-five grid stencil and tangent-plane frame (msrotm);
//   * 4x4 block of doublet parameters on the staggered parameter grid;
//   * weighted quadratic least-squares fit for panel corners (16 points)
//     and for the two adjacent edge midpoints (12 points);
//   * one-dimensional edgls() fit along network edges;
//   * network-level sparse B_D assembly for all panels, plus the dense 9x25
//     deep-interior form used by the first implementation.

#include "OuterSplines.h"

#include "AeroModel.h"
#include "EdgeSplines.h"
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
            const int mm = std::clamp(m + i - 2, 1, net.nm);
            const int nn = std::clamp(n + j - 2, 1, net.nn);
            out[static_cast<std::size_t>(i + 5 * j)] =
                networkPoint(net, mm, nn);
        }
    }
    return true;
}

bool buildParameterBlock4(const Network& net,
                          const NetworkDoubletParameters& parameters,
                          int m, int n,
                          std::array<Eigen::Vector3d, 16>& points,
                          std::array<int, 16>& ids) {
    const int nParamM = net.nm + 1;
    const int nParamN = net.nn + 1;
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            // blccal() clamps the (nm+1)x(nn+1) singularity-parameter
            // indices to the network edge rows/columns.
            const int mi = std::clamp(m + i - 1, 1, nParamM);
            const int nj = std::clamp(n + j - 1, 1, nParamN);
            const int id = parameters.parameterId(mi, nj);
            if (id < 0 || id >= parameters.count ||
                static_cast<std::size_t>(id) >= parameters.location.size())
                return false;
            const std::size_t k = static_cast<std::size_t>(i + 4 * j);
            points[k] = parameters.location[static_cast<std::size_t>(id)];
            ids[k] = id;
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
    if (std::abs(M.determinant()) > 1.0e-30) {
        const Eigen::Matrix<double, 6, Eigen::Dynamic> pinv =
            M.inverse() * A.transpose() * W2;
        return pinv.row(0).transpose();
    }

    // Rank-deficient stencil (e.g. the one-panel-wide wingtip network).
    // Fall back to the generalized inverse of the weighted basis matrix,
    // which is what lsqsf/sginvx use in panair.f90.
    const Eigen::MatrixXd weighted = wv.asDiagonal() * A;
    Eigen::JacobiSVD<Eigen::MatrixXd> svd(
        weighted, Eigen::ComputeThinU | Eigen::ComputeThinV);
    const Eigen::VectorXd s = svd.singularValues();
    if (s.size() == 0 || s(0) <= 0.0)
        return Eigen::VectorXd::Zero(n);

    const double tol = 1.0e-12 * s(0);
    Eigen::VectorXd row = Eigen::VectorXd::Zero(n);
    for (int k = 0; k < s.size(); ++k) {
        if (s(k) <= tol)
            continue;
        const double v = svd.matrixV()(0, k) / s(k);
        for (int j = 0; j < n; ++j)
            row(j) += wv(j) * svd.matrixU()(j, k) * v;
    }
    return row;
}

bool cornerRow(const Network& net,
               const NetworkDoubletParameters& parameters, int m, int n,
               const Compressibility& flow, std::array<double, 16>& b,
               std::array<int, 16>& ids) {
    std::array<Eigen::Vector3d, 25> stencil;
    std::array<Eigen::Vector3d, 16> params;
    if (!buildStencil5(net, m, n, stencil) ||
        !buildParameterBlock4(net, parameters, m, n, params, ids))
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

bool midColRow(const Network& net,
               const NetworkDoubletParameters& parameters, int m, int n,
               const Compressibility& flow, std::array<double, 12>& b,
               std::array<int, 12>& ids) {
    std::array<Eigen::Vector3d, 25> stencil;
    std::array<Eigen::Vector3d, 16> params;
    std::array<int, 16> blockIds;
    if (!buildStencil5(net, m, n, stencil) ||
        !buildParameterBlock4(net, parameters, m, n, params, blockIds))
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
            const std::size_t pk =
                static_cast<std::size_t>(i + 4 * (j + 1));
            const Eigen::Vector3d p = params[pk];
            xy[static_cast<std::size_t>(k)] = projectLocal(frame, origin, p);
            w[static_cast<std::size_t>(k)] = lsqWeight(flow, p - origin);
            ids[static_cast<std::size_t>(k)] = blockIds[pk];
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

bool midRowRow(const Network& net,
               const NetworkDoubletParameters& parameters, int m, int n,
               const Compressibility& flow, std::array<double, 12>& b,
               std::array<int, 12>& ids) {
    std::array<Eigen::Vector3d, 25> stencil;
    std::array<Eigen::Vector3d, 16> params;
    std::array<int, 16> blockIds;
    if (!buildStencil5(net, m, n, stencil) ||
        !buildParameterBlock4(net, parameters, m, n, params, blockIds))
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
            const std::size_t pk =
                static_cast<std::size_t>((i + 1) + 4 * j);
            const Eigen::Vector3d p = params[pk];
            xy[static_cast<std::size_t>(k)] = projectLocal(frame, origin, p);
            w[static_cast<std::size_t>(k)] = lsqWeight(flow, p - origin);
            ids[static_cast<std::size_t>(k)] = blockIds[pk];
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

// One-dimensional edge fit of blccal(): the value at edge grid point g is a
// weighted quadratic fit to the four neighbouring edge parameters.  Slots
// are clamped at the ends of the edge, exactly as zc()/zmc() are in Fortran;
// coincident parameter ids are combined by DoubletRow::add().
bool edgeRow(const Network& net,
             const NetworkDoubletParameters& parameters,
             const Compressibility& flow, int edge, int g, DoubletRow& row) {
    if (edge < 0 || edge > 3)
        return false;

    const int nEdge = (edge == 0 || edge == 2) ? parameters.nPanelN
                                               : parameters.nPanelM;
    if (g < 0 || g > nEdge)
        return false;

    auto edgePoint = [&](int index) -> Eigen::Vector3d {
        switch (edge) {
        case 0: return networkPoint(net, 1, index + 1);       // top row
        case 1: return networkPoint(net, index + 1, net.nn);  // right column
        case 2: return networkPoint(net, net.nm, index + 1);  // bottom row
        default: return networkPoint(net, index + 1, 1);      // left column
        }
    };

    const std::array<int, 4> startCorner{{0, 1, 3, 0}};
    const std::array<int, 4> endCorner{{1, 2, 2, 3}};

    auto parameterId = [&](int p) -> int {
        if (p <= 1)
            return parameters.cornerId[static_cast<std::size_t>(
                startCorner[static_cast<std::size_t>(edge)])];
        if (p >= nEdge + 2)
            return parameters.cornerId[static_cast<std::size_t>(
                endCorner[static_cast<std::size_t>(edge)])];
        return parameters.edgeMidBase[static_cast<std::size_t>(edge)] +
               (p - 2);
    };

    const Eigen::Vector3d origin = edgePoint(g);
    std::array<double, 4> d{};
    std::array<double, 4> w{};
    std::array<int, 4> ids{};

    for (int k = 0; k < 4; ++k) {
        const int p = std::clamp(g + k, 1, nEdge + 2);
        const int id = parameterId(p);
        if (id < 0 || id >= parameters.count ||
            static_cast<std::size_t>(id) >= parameters.location.size())
            return false;

        const int i0 = std::clamp(p - 2, 0, nEdge);
        const int i1 = std::clamp(p - 1, 0, nEdge);
        d[static_cast<std::size_t>(k)] =
            (edgePoint(i1) - edgePoint(i0)).norm();

        const Eigen::Vector3d location =
            parameters.location[static_cast<std::size_t>(id)];
        w[static_cast<std::size_t>(k)] = lsqWeight(flow, location - origin);
        ids[static_cast<std::size_t>(k)] = id;
    }

    const EdgeLeastSquares ls = edgls(d, w);
    for (int k = 0; k < 4; ++k)
        row.add(ids[static_cast<std::size_t>(k)],
                ls.a[static_cast<std::size_t>(k)]);
    return true;
}

} // namespace

void DoubletRow::add(int id, double value) {
    for (std::size_t k = 0; k < index.size(); ++k) {
        if (index[k] == id) {
            coefficient[k] += value;
            return;
        }
    }
    index.push_back(id);
    coefficient.push_back(value);
}

int NetworkDoubletParameters::parameterId(int m, int n) const {
    const int nm = nPanelM + 1;
    const int nn = nPanelN + 1;
    if (m < 1 || m > nm + 1 || n < 1 || n > nn + 1)
        return -1;

    if (m == 1 && n == 1)
        return cornerId[0];
    if (m == 1 && n == nn + 1)
        return cornerId[1];
    if (m == nm + 1 && n == nn + 1)
        return cornerId[2];
    if (m == nm + 1 && n == 1)
        return cornerId[3];

    if (m == 1)
        return edgeMidBase[0] + (n - 2);
    if (n == nn + 1)
        return edgeMidBase[1] + (m - 2);
    if (m == nm + 1)
        return edgeMidBase[2] + (n - 2);
    if (n == 1)
        return edgeMidBase[3] + (m - 2);

    return centerId(m - 2, n - 2);
}

NetworkDoubletParameters buildNetworkDoubletParameters(
    const Network& network) {
    NetworkDoubletParameters out;
    const int nm = network.nm;
    const int nn = network.nn;
    if (nm < 2 || nn < 2)
        return out;

    const std::size_t nPoints =
        static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
    if (network.coordinates.size() < nPoints * 3)
        return out;

    out.nPanelM = nm - 1;
    out.nPanelN = nn - 1;
    out.nCenters = out.nPanelM * out.nPanelN;

    // Panel centres.
    out.count = out.nCenters;
    out.location.resize(static_cast<std::size_t>(out.count));
    for (int pj = 0; pj < out.nPanelN; ++pj) {
        for (int pi = 0; pi < out.nPanelM; ++pi) {
            const PanelGeometry panel = buildPanelGeometry(
                network.coordinates, nm, nn, pi, pj);
            out.location[static_cast<std::size_t>(out.centerId(pi, pj))] =
                panel.cp[8];
        }
    }

    // Network corners.
    out.cornerId[0] = out.count;
    out.cornerId[1] = out.count + 1;
    out.cornerId[2] = out.count + 2;
    out.cornerId[3] = out.count + 3;
    out.count += 4;
    out.location.resize(static_cast<std::size_t>(out.count));
    out.location[static_cast<std::size_t>(out.cornerId[0])] =
        networkPoint(network, 1, 1);
    out.location[static_cast<std::size_t>(out.cornerId[1])] =
        networkPoint(network, 1, nn);
    out.location[static_cast<std::size_t>(out.cornerId[2])] =
        networkPoint(network, nm, nn);
    out.location[static_cast<std::size_t>(out.cornerId[3])] =
        networkPoint(network, nm, 1);

    // Edge midpoints, grouped by network edge.
    out.edgeMidBase[0] = out.count;
    out.count += out.nPanelN;
    out.edgeMidBase[1] = out.count;
    out.count += out.nPanelM;
    out.edgeMidBase[2] = out.count;
    out.count += out.nPanelN;
    out.edgeMidBase[3] = out.count;
    out.count += out.nPanelM;
    out.location.resize(static_cast<std::size_t>(out.count));

    for (int j = 0; j < out.nPanelN; ++j) {
        out.location[static_cast<std::size_t>(out.edgeMidBase[0] + j)] =
            0.5 * (networkPoint(network, 1, j + 1) +
                   networkPoint(network, 1, j + 2));
    }
    for (int i = 0; i < out.nPanelM; ++i) {
        out.location[static_cast<std::size_t>(out.edgeMidBase[1] + i)] =
            0.5 * (networkPoint(network, i + 1, nn) +
                   networkPoint(network, i + 2, nn));
    }
    for (int j = 0; j < out.nPanelN; ++j) {
        out.location[static_cast<std::size_t>(out.edgeMidBase[2] + j)] =
            0.5 * (networkPoint(network, nm, j + 1) +
                   networkPoint(network, nm, j + 2));
    }
    for (int i = 0; i < out.nPanelM; ++i) {
        out.location[static_cast<std::size_t>(out.edgeMidBase[3] + i)] =
            0.5 * (networkPoint(network, i + 1, 1) +
                   networkPoint(network, i + 2, 1));
    }

    return out;
}

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

    const NetworkDoubletParameters parameters =
        buildNetworkDoubletParameters(net);
    if (parameters.count <= 0)
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
    std::array<int, 16> cornerIds{};
    for (int lcp = 1; lcp <= 4; ++lcp) {
        const int mc = ((lcp - 1) % 2) + 1;
        const int nc = ((lcp - 1) / 2) + 1;
        if (!cornerRow(net, parameters, pRow + mc - 1, pCol + nc - 1, flow,
                       corner[static_cast<std::size_t>(lcp - 1)], cornerIds))
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
    std::array<int, 12> midIds{};
    if (!midColRow(net, parameters, pRow, pCol, flow, colMid1, midIds))
        return DoubletOuterSpline{};
    if (!midColRow(net, parameters, pRow + 1, pCol, flow, colMid2, midIds))
        return DoubletOuterSpline{};
    if (!midRowRow(net, parameters, pRow, pCol, flow, rowMid1, midIds))
        return DoubletOuterSpline{};
    if (!midRowRow(net, parameters, pRow, pCol + 1, flow, rowMid3, midIds))
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

DoubletPanelSpline buildDoubletPanelSpline(
    const Network& net, const Compressibility& flow, int panelI, int panelJ,
    const NetworkDoubletParameters& parameters) {
    DoubletPanelSpline out;
    const int nm = net.nm;
    const int nn = net.nn;
    const int nPanelM = nm - 1;
    const int nPanelN = nn - 1;
    if (nm < 2 || nn < 2 || panelI < 0 || panelJ < 0 ||
        panelI >= nPanelM || panelJ >= nPanelN ||
        parameters.nPanelM != nPanelM ||
        parameters.nPanelN != nPanelN || parameters.count <= 0 ||
        parameters.location.size() <
            static_cast<std::size_t>(parameters.count))
        return out;

    const std::size_t nPoints =
        static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
    if (net.coordinates.size() < nPoints * 3)
        return out;

    auto addUnit = [&](int row, int id) {
        out.rows[static_cast<std::size_t>(row)].add(id, 1.0);
    };

    auto addCornerPoint = [&](int row, int m, int n) -> bool {
        // Network corners are explicit doublet parameters.
        if (m == 0 && n == 0) {
            addUnit(row, parameters.cornerId[0]);
            return true;
        }
        if (m == 0 && n == nPanelN) {
            addUnit(row, parameters.cornerId[1]);
            return true;
        }
        if (m == nPanelM && n == nPanelN) {
            addUnit(row, parameters.cornerId[2]);
            return true;
        }
        if (m == nPanelM && n == 0) {
            addUnit(row, parameters.cornerId[3]);
            return true;
        }

        // Interior network-edge grid point: one-dimensional quadratic fit.
        if (m == 0 || m == nPanelM || n == 0 || n == nPanelN) {
            const int edge = (m == 0)   ? 0
                             : (m == nPanelM) ? 2
                             : (n == 0)  ? 3
                                         : 1;
            const int g = (edge == 0 || edge == 2) ? n : m;
            DoubletRow fit;
            if (!edgeRow(net, parameters, flow, edge, g, fit))
                return false;
            out.rows[static_cast<std::size_t>(row)] = fit;
            return true;
        }

        std::array<double, 16> b{};
        std::array<int, 16> ids{};
        if (!cornerRow(net, parameters, m + 1, n + 1, flow, b, ids))
            return false;
        for (int k = 0; k < 16; ++k)
            out.rows[static_cast<std::size_t>(row)].add(
                ids[static_cast<std::size_t>(k)],
                b[static_cast<std::size_t>(k)]);
        return true;
    };

    auto addMidColumn = [&](int row, int m, int n) -> bool {
        std::array<double, 12> b{};
        std::array<int, 12> ids{};
        if (!midColRow(net, parameters, m, n, flow, b, ids))
            return false;
        for (int k = 0; k < 12; ++k)
            out.rows[static_cast<std::size_t>(row)].add(
                ids[static_cast<std::size_t>(k)],
                b[static_cast<std::size_t>(k)]);
        return true;
    };

    auto addMidRow = [&](int row, int m, int n) -> bool {
        std::array<double, 12> b{};
        std::array<int, 12> ids{};
        if (!midRowRow(net, parameters, m, n, flow, b, ids))
            return false;
        for (int k = 0; k < 12; ++k)
            out.rows[static_cast<std::size_t>(row)].add(
                ids[static_cast<std::size_t>(k)],
                b[static_cast<std::size_t>(k)]);
        return true;
    };

    // P1..P4 (0-based network grid points).
    if (!addCornerPoint(0, panelI, panelJ) ||
        !addCornerPoint(1, panelI, panelJ + 1) ||
        !addCornerPoint(2, panelI + 1, panelJ + 1) ||
        !addCornerPoint(3, panelI + 1, panelJ))
        return DoubletPanelSpline{};

    // P5..P8.  A panel edge lying on a network edge has its midpoint as an
    // explicit doublet parameter; otherwise use the 12-point quadratic fit.
    if (panelI == 0)
        addUnit(4, parameters.edgeMidBase[0] + panelJ);
    else if (!addMidColumn(4, panelI + 1, panelJ + 1))
        return DoubletPanelSpline{};

    if (panelJ == nPanelN - 1)
        addUnit(5, parameters.edgeMidBase[1] + panelI);
    else if (!addMidRow(5, panelI + 1, panelJ + 2))
        return DoubletPanelSpline{};

    if (panelI == nPanelM - 1)
        addUnit(6, parameters.edgeMidBase[2] + panelJ);
    else if (!addMidColumn(6, panelI + 2, panelJ + 1))
        return DoubletPanelSpline{};

    if (panelJ == 0)
        addUnit(7, parameters.edgeMidBase[3] + panelI);
    else if (!addMidRow(7, panelI + 1, panelJ + 1))
        return DoubletPanelSpline{};

    // P9 (panel centre).
    addUnit(8, parameters.centerId(panelI, panelJ));

    for (const DoubletRow& row : out.rows) {
        if (row.index.empty() || row.index.size() != row.coefficient.size())
            return DoubletPanelSpline{};
    }

    out.valid = true;
    return out;
}

NetworkDoubletSpline buildNetworkDoubletSpline(
    const Network& net, const Compressibility& flow,
    const NetworkDoubletParameters& parameters) {
    NetworkDoubletSpline out;
    const int nPanelM = net.nm - 1;
    const int nPanelN = net.nn - 1;
    if (nPanelM <= 0 || nPanelN <= 0 ||
        parameters.nPanelM != nPanelM ||
        parameters.nPanelN != nPanelN)
        return out;

    out.nPanelM = nPanelM;
    out.nPanelN = nPanelN;
    out.nParameters = parameters.count;
    out.rows.reserve(static_cast<std::size_t>(9) *
                     static_cast<std::size_t>(nPanelM) *
                     static_cast<std::size_t>(nPanelN));

    for (int panelJ = 0; panelJ < nPanelN; ++panelJ) {
        for (int panelI = 0; panelI < nPanelM; ++panelI) {
            const DoubletPanelSpline panel = buildDoubletPanelSpline(
                net, flow, panelI, panelJ, parameters);
            if (!panel.valid)
                return NetworkDoubletSpline{};
            for (const DoubletRow& row : panel.rows)
                out.rows.push_back(row);
        }
    }

    out.valid = true;
    return out;
}

} // namespace a502
