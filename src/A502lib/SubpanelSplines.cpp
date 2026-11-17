// SubpanelSplines.cpp
// -------------------
// Source distribution on the eight triangular subpanels.

#include "SubpanelSplines.h"

#include <cmath>

namespace a502 {

namespace {

constexpr double kEps = 1.0e-14;

// Value of the panel source distribution at canonical point P_{index+1} as a
// linear combination of [sigma1, sigma2, sigma3, sigma4, sigma9].
Eigen::Matrix<double, 1, 5> sourceBasis(int cpIndex) {
    Eigen::Matrix<double, 1, 5> b =
        Eigen::Matrix<double, 1, 5>::Zero();
    switch (cpIndex) {
    case 0: b(0) = 1.0; break;                 // P1
    case 1: b(1) = 1.0; break;                 // P2
    case 2: b(2) = 1.0; break;                 // P3
    case 3: b(3) = 1.0; break;                 // P4
    case 4: b(0) = 0.5; b(1) = 0.5; break;     // P5 = (P1+P2)/2
    case 5: b(1) = 0.5; b(2) = 0.5; break;     // P6 = (P2+P3)/2
    case 6: b(2) = 0.5; b(3) = 0.5; break;     // P7 = (P3+P4)/2
    case 7: b(3) = 0.5; b(0) = 0.5; break;     // P8 = (P4+P1)/2
    case 8: b(4) = 1.0; break;                 // P9
    default: break;
    }
    return b;
}

} // namespace

std::array<SourceSubpanelSpline, 8>
buildSourceSubpanelSplines(const PanelGeometry& panel) {
    std::array<SourceSubpanelSpline, 8> out;

    if (panel.degenerate)
        return out;

    for (int s = 0; s < 8; ++s) {
        const SubpanelGeometry& sp = panel.subpanels[s];
        SourceSubpanelSpline& spline = out[s];
        if (sp.degenerate)
            continue;

        const Eigen::Matrix3d& frame =
            sp.scaledValid ? sp.scaledToLocal : sp.toLocal;

        Eigen::Matrix3d M;
        Eigen::Matrix<double, 3, 5> V;
        for (int k = 0; k < 3; ++k) {
            const Eigen::Vector3d local =
                frame * (sp.vertex[k] - sp.origin);
            M(k, 0) = 1.0;
            M(k, 1) = local.x();
            M(k, 2) = local.y();
            V.row(k) = sourceBasis(sp.cpIndex[k]);
        }

        const double det = M.determinant();
        if (std::abs(det) <= kEps)
            continue;

        spline.coefficients = M.inverse() * V;
        spline.valid = true;
    }

    return out;
}

Eigen::Vector3d sourceTaylorCoefficients(
    const SourceSubpanelSpline& spline, const PanelSourceValues& values) {
    return spline.coefficients * values;
}

double evaluateSource(const SourceSubpanelSpline& spline,
                      const Eigen::Vector2d& local,
                      const PanelSourceValues& values) {
    if (!spline.valid)
        return 0.0;
    const Eigen::Vector3d c = sourceTaylorCoefficients(spline, values);
    return c(0) + c(1) * local.x() + c(2) * local.y();
}

} // namespace a502
