#pragma once

// FlowSolution.h
// --------------
// Surface velocity and pressure coefficient from a solved singularity
// vector.  The present Cp rule is the linear (first-order) rule of
// cpcalx(); the second-order and isentropic rules will be added with the
// material-property input.

#define EIGEN_NO_DEBUG

#include <vector>

#include <Eigen/Dense>

#include "GlobalBoundarySystem.h"

namespace a502 {

struct SurfaceFlowResult {
    std::vector<Eigen::Vector3d> velocity;
    std::vector<double> cp;
    bool valid{false};
};

// Linear pressure coefficient for the total velocity:
//   cp = -2 (V - Vinf) . Vinf / |Vinf|^2.
double linearPressureCoefficient(const Eigen::Vector3d& velocity,
                                 const Eigen::Vector3d& freestream);

// Evaluates potential/velocity rows at all control points of all networks
// and applies the linear Cp rule.
SurfaceFlowResult evaluateSurfaceFlow(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout,
    const GlobalParameterValues& values, bool includeEdgeAndCorner = true);

} // namespace a502
