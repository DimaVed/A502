// FlowSolution.cpp
// ----------------
// Surface velocity/Cp reconstruction from the solved singularity vector.

#include "FlowSolution.h"

#include "AeroModel.h"
#include "ControlPoints.h"

#include <algorithm>
#include <cmath>

namespace a502 {

double pressureCoefficient(
    const Eigen::Vector3d& perturbationVelocity,
    const Eigen::Vector3d& freestream, double mach, CpFormula formula,
    const CpMaterialRatios& material) {
    const double f2 = freestream.squaredNorm();
    if (f2 <= 1.0e-30)
        return 0.0;

    // pv = perturbationVelocity / |freestream|; pvf is its component
    // along the freestream, pvs the squared magnitude.
    const double pvf = freestream.dot(perturbationVelocity) / f2;
    const double pvs = perturbationVelocity.squaredNorm() / f2;

    switch (formula) {
    case CpFormula::Linear:
        return material.cpFactor * (-2.0 * pvf);
    case CpFormula::SlenderBody:
        return material.cpFactor *
               (-2.0 * pvf - pvs + pvf * pvf);
    case CpFormula::SecondOrder: {
        const double mach2 = mach * mach;
        return material.cpFactor *
               (-2.0 * pvf - pvs + mach2 * pvf * pvf);
    }
    case CpFormula::Isentropic: {
        const double fac =
            material.velocity * material.velocity / material.temperature;
        double cp = 1.0 - fac * (1.0 + 2.0 * pvf + pvs);
        if (mach > 0.0) {
            const double mach2 = mach * mach;
            const double base = std::max(0.0, 1.0 + 0.2 * mach2 * cp);
            cp = (1.4285714285714286 / mach2) *
                 (material.pressure * std::pow(base, 3.5) - 1.0);
        }
        return cp;
    }
    }
    return 0.0;
}

double linearPressureCoefficient(const Eigen::Vector3d& velocity,
                                 const Eigen::Vector3d& freestream) {
    return pressureCoefficient(velocity - freestream, freestream, 0.0,
                               CpFormula::Linear);
}

SurfaceFlowResult evaluateSurfaceFlow(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout,
    const GlobalParameterValues& values, bool includeEdgeAndCorner,
    CpFormula formula, const CpMaterialRatios& material) {
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
        out.cp[k] = pressureCoefficient(perturbation, flow.compd,
                                        flow.mach, formula, material);
    }

    out.valid = true;
    return out;
}

} // namespace a502
