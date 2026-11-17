#pragma once

// BoundaryConditions.h
// --------------------
// Core boundary-condition row assembly from aical() in panair.f90.
//
// The full A502 option processing (bconcl/bcopt) is not ported yet.  This
// module covers the nct=1 (normal mass flux only) and nct=4 (potential only)
// branches for a panel-centre control point, given the user coefficients
// cu,cl,du,dl and right-hand side bet.

#define EIGEN_NO_DEBUG

#include <vector>

#include <Eigen/Dense>

#include "GlobalInfluence.h"

namespace a502 {

struct Network;

struct BoundaryConditionCoefficients {
    int nct{1}; // 1 = normal mass flux only, 4 = potential only
    double cu{0.0}, cl{0.0}; // upper/lower normal mass flux
    double du{0.0}, dl{0.0}; // upper/lower potential
    double bet{0.0};         // user right-hand side
};

// Dependence of the source and doublet values at the panel-centre control
// point on the local source/doublet parameters of that panel.
struct LocalSingularityTerms {
    std::vector<int> sourceId;
    Eigen::VectorXd sourceValue;
    std::vector<int> doubletId;
    Eigen::VectorXd doubletValue;
    bool valid{false};
};

LocalSingularityTerms buildPanelCenterLocalTerms(
    const Network& network, const Compressibility& flow, int panelI,
    int panelJ, const NetworkDoubletParameters& doubletParameters);

struct BoundaryRow {
    Eigen::VectorXd source;
    Eigen::VectorXd doublet;
    double rhs{0.0};
    bool valid{false};
};

BoundaryRow buildBoundaryRow(
    const GlobalInfluence& influence, const Eigen::Vector3d& normal,
    const Compressibility& flow,
    const BoundaryConditionCoefficients& coefficients,
    const LocalSingularityTerms& local);

} // namespace a502
