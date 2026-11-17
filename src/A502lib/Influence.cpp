// Influence.cpp
// -------------
// Steady-flow port of influ() and the compressible diameter from geomc().

#include "Influence.h"

#include "PanelGeometry.h"

#include <algorithm>
#include <cmath>

namespace a502 {

namespace {

// compip(): inner product with the compressibility-direction component
// scaled by bs.
double compressibleInnerProduct(const Eigen::Vector3d& u,
                                const Eigen::Vector3d& v,
                                const Eigen::Vector3d& direction,
                                double scale) {
    const double c2 = direction.squaredNorm();
    if (c2 <= 1.0e-30)
        return scale * u.dot(v);
    return (1.0 - scale) * u.dot(direction) * v.dot(direction) / c2 +
           scale * u.dot(v);
}

} // namespace

double compressedPanelDiameter(const PanelGeometry& panel,
                               const Compressibility& flow) {
    double diameter = 0.0;
    for (int i = 0; i < 4; ++i) {
        const Eigen::Vector3d d = panel.cp[i] - panel.cp[8];
        const double ww = compressibleInnerProduct(
            d, d, flow.compd, flow.abetms);
        diameter = std::max(diameter, std::sqrt(std::max(ww, 0.0)));
    }
    return 2.0 * diameter;
}

int influenceLevel(const PanelGeometry& panel, const Compressibility& flow,
                   const Eigen::Vector3d& fieldPoint) {
    if (panel.degenerate)
        return 6;

    constexpr double delt1 = 0.042;
    constexpr double delt2 = 0.125;
    constexpr double delt3 = 0.2;
    constexpr double delt4 = 0.4;
    constexpr double delt5 = 0.75;

    const Eigen::Vector3d pc = fieldPoint - panel.cp[8];
    double pcm = compressibleInnerProduct(
        pc, pc, flow.compd, flow.betams);
    if (pcm <= 0.0)
        return 6; // influ() leaves iflu at its near-field default

    pcm = 1.0 / std::sqrt(pcm);
    const double pcn = compressibleInnerProduct(
        pc, pc, flow.compd, flow.abetms);
    const double eps = 0.5 * compressedPanelDiameter(panel, flow) *
                       pcm * pcm * std::sqrt(std::max(pcn, 0.0));

    if (eps < delt1) return 1;
    if (eps < delt2) return 2;
    if (eps < delt3) return 3;
    if (eps < delt4) return 4;
    if (eps < delt5) return 5;
    return 6;
}

} // namespace a502
