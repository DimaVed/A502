// FlowSolution.cpp
// ----------------
// Surface velocity/Cp reconstruction from the solved singularity vector.

#include "FlowSolution.h"

#include "AeroModel.h"
#include "ControlPoints.h"

#include <cmath>

namespace a502 {

double linearPressureCoefficient(const Eigen::Vector3d& velocity,
                                 const Eigen::Vector3d& freestream) {
    const double f2 = freestream.squaredNorm();
    if (f2 <= 1.0e-30)
        return 0.0;
    return -2.0 * (velocity.dot(freestream) - f2) / f2;
}

SurfaceFlowResult evaluateSurfaceFlow(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout,
    const GlobalParameterValues& values, bool includeEdgeAndCorner) {
    SurfaceFlowResult out;
    if (!values.valid || networks.empty() ||
        layout.network.size() != networks.size() ||
        values.source.size() != layout.totalSource ||
        values.doublet.size() != layout.totalDoublet)
        return out;

    std::vector<Eigen::Vector3d> points;
    for (const Network& net : networks) {
        const std::vector<ControlPoint> cps = buildNetworkControlPoints(net);
        for (const ControlPoint& cp : cps) {
            if (!includeEdgeAndCorner &&
                cp.kind != ControlPointKind::PanelCenter)
                continue;
            points.push_back(cp.position);
        }
    }
    if (points.empty())
        return out;

    const auto globalAic =
        buildGlobalNearFieldAic(networks, flow, layout, points);
    if (!globalAic.valid)
        return SurfaceFlowResult{};

    out.velocity.resize(points.size());
    out.cp.resize(points.size());
    for (std::size_t k = 0; k < points.size(); ++k) {
        const Eigen::Index row = 4 * static_cast<Eigen::Index>(k);
        const Eigen::Vector3d perturbation =
            globalAic.source.block(row + 1, 0, 3, globalAic.nSource) *
                values.source +
            globalAic.doublet.block(row + 1, 0, 3, globalAic.nDoublet) *
                values.doublet;
        const Eigen::Vector3d velocity = flow.compd + perturbation;
        out.velocity[k] = velocity;
        out.cp[k] = linearPressureCoefficient(velocity, flow.compd);
    }

    out.valid = true;
    return out;
}

} // namespace a502
