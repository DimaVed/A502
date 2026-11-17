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

#include "ControlPoints.h"
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

// Standard nropt=2 right-hand side: minus the freestream normal component.
BoundaryConditionCoefficients standardNormalMassFlux(
    double cu, double cl, const Eigen::Vector3d& normal,
    const Eigen::Vector3d& freestream);

struct BoundarySystemRows {
    int nRows{0};
    int nSource{0};
    int nDoublet{0};
    Eigen::MatrixXd source;
    Eigen::MatrixXd doublet;
    Eigen::VectorXd rhs;
    std::vector<int> controlPointIndex; // index into the input list
    bool valid{false};
};

// Assemble one normal-mass-flux boundary row per panel-centre control point
// of a network.  Edge/corner control points require the matching/abutment
// equations and are not included here.
BoundarySystemRows buildNetworkNormalFluxRows(
    const Network& network, const Compressibility& flow,
    const NetworkSingularityLayout& layout,
    const std::vector<ControlPoint>& controlPoints, double cu, double cl);

} // namespace a502
