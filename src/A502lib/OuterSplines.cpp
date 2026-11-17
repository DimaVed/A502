// OuterSplines.cpp
// ----------------
// Source outer spline for interior analysis panels.

#include "OuterSplines.h"

#include "AeroModel.h"
#include "PanelGeometry.h"
#include "PanelParameters.h"

#include <cmath>

namespace a502 {

SourceOuterSpline buildSourceOuterSpline(const Network& network,
                                         int panelI, int panelJ) {
    SourceOuterSpline out;

    const int nm = network.nm;
    const int nn = network.nn;
    if (nm < 2 || nn < 2 || panelI < 0 || panelJ < 0 ||
        panelI + 1 >= nm || panelJ + 1 >= nn)
        return out;

    const std::size_t nPoints =
        static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
    if (network.coordinates.size() < nPoints * 3)
        return out;

    const int nPanelM = nm - 1;

    const PanelGeometry panel = buildPanelGeometry(
        network.coordinates, nm, nn, panelI, panelJ);
    if (panel.degenerate)
        return out;

    const std::array<Eigen::Vector3d, 4> corners = {
        panel.cp[0], panel.cp[1], panel.cp[2], panel.cp[3]};
    const NearPlane np = buildNearPlane(corners);
    if (!np.valid)
        return out;

    const SubpanelGeometry& inner = panel.subpanels[4];
    const Eigen::Matrix3d innerFrame =
        inner.scaledValid ? inner.scaledToLocal : inner.toLocal;
    const Eigen::Vector3d center = panel.cp[8];

    // lproj(): project a point onto the near plane while preserving its
    // distance from the panel center.
    auto projectToNearPlane = [&](const Eigen::Vector3d& p)
        -> Eigen::Vector3d {
        Eigen::Vector3d d = p - center;
        const double len = d.norm();
        if (len <= 1.0e-15)
            return center;
        const double n2 = np.normal.squaredNorm();
        if (n2 > 1.0e-30)
            d -= np.normal * (d.dot(np.normal) / n2);
        if (d.norm() <= 1.0e-15)
            return center;
        d.normalize();
        return center + len * d;
    };

    Eigen::Matrix<double, 9, 3> A;
    Eigen::Matrix<double, 9, 1> w;

    int k = 0;
    for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
            const int ni = panelI + di;
            const int nj = panelJ + dj;
            if (ni < 0 || ni >= nPanelM || nj < 0 || nj >= nn - 1) {
                // sing() simply skips parameters outside the panel-center
                // grid.  Keep a zero-weight row so the 5x9 layout is stable.
                A(k, 0) = 1.0;
                A(k, 1) = 0.0;
                A(k, 2) = 0.0;
                w(k) = 0.0;
                out.neighborPanel[static_cast<std::size_t>(k)] = -1;
                ++k;
                continue;
            }

            const PanelGeometry neighbor = buildPanelGeometry(
                network.coordinates, nm, nn, ni, nj);
            if (neighbor.degenerate) {
                A(k, 0) = 1.0;
                A(k, 1) = 0.0;
                A(k, 2) = 0.0;
                w(k) = 0.0;
                out.neighborPanel[static_cast<std::size_t>(k)] = -1;
                ++k;
                continue;
            }

            const Eigen::Vector3d projected =
                projectToNearPlane(neighbor.cp[8]);
            const Eigen::Vector3d local =
                innerFrame * (projected - center);

            A(k, 0) = 1.0;
            A(k, 1) = local.x();
            A(k, 2) = local.y();
            w(k) = (di == 0 && dj == 0) ? 1.0e4 : 1.0;

            out.neighborPanel[static_cast<std::size_t>(k)] =
                ni + nPanelM * nj;
            ++k;
        }
    }

    const Eigen::Matrix<double, 9, 1> w2 = w.array().square().matrix();
    const Eigen::Matrix3d M = A.transpose() * w2.asDiagonal() * A;
    if (std::abs(M.determinant()) <= 1.0e-30)
        return out;

    const Eigen::Matrix<double, 3, 9> asts =
        M.inverse() * A.transpose() * w2.asDiagonal();

    // Panel source values at P1..P4,P9, in the inner-frame coordinates.
    for (int row = 0; row < 4; ++row) {
        const Eigen::Vector3d local = innerFrame * (panel.cp[row] - center);
        Eigen::Matrix<double, 1, 3> q;
        q << 1.0, local.x(), local.y();
        out.coefficients.row(row) = q * asts;
    }
    Eigen::Matrix<double, 1, 3> qc;
    qc << 1.0, 0.0, 0.0;
    out.coefficients.row(4) = qc * asts;

    out.valid = true;
    return out;
}

} // namespace a502
