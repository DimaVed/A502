#pragma once

// GlobalBoundarySystem.h
// ----------------------
// Assemble the first global system from:
//   * one standard normal-mass-flux row per panel-centre control point;
//   * full-edge doublet matching rows.
// Inactive source/doublet types are removed from the unknown vector using
// the kt mapping stored in GlobalSingularityLayout.

#define EIGEN_NO_DEBUG

#include <vector>

#include <Eigen/Dense>

#include "GlobalInfluence.h"

namespace a502 {

struct Network;

struct GlobalBoundarySystem {
    int nBcRows{0};
    int nMatchingRows{0};
    int nRows{0};
    int nActiveSource{0};
    int nActiveDoublet{0};
    int nColumns{0};
    Eigen::MatrixXd matrix;
    Eigen::VectorXd rhs;
    bool valid{false};
};

GlobalBoundarySystem buildGlobalNormalFluxSystem(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout, double cu, double cl,
    double matchingTolerance);

} // namespace a502
