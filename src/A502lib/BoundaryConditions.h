#pragma once

// BoundaryConditions.h
// --------------------
// Core boundary-condition row assembly from aical() in panair.f90.
//
// The full A502 option processing (bconcl/bcopt) is not ported yet.  This
// module covers the nct=1 (normal mass flux only), nct=2 (general) and nct=4
// (potential only) branches for a panel-centre control point, given the user
// coefficients cu,cl,tu,tl,du,dl and right-hand side bet.

#define EIGEN_NO_DEBUG

#include <vector>

#include <Eigen/Dense>

#include "ControlPoints.h"
#include "GlobalInfluence.h"

namespace a502 {

struct Network;

struct BoundaryConditionCoefficients {
    int nct{1}; // 1 = normal mass flux, 2 = general, 4 = potential
    double cu{0.0}, cl{0.0}; // upper/lower normal mass flux
    Eigen::Vector3d tu{0.0, 0.0, 0.0};
    Eigen::Vector3d tl{0.0, 0.0, 0.0};
    double du{0.0}, dl{0.0}; // upper/lower potential
    double bet{0.0};         // user right-hand side
};

// Dependence of the source and doublet values at the panel-centre control
// point on the local source/doublet parameters of that panel.  doubletValue
// has rows [mu, dmu/dx, dmu/dy, dmu/dz] in global coordinates.  hasSource /
// hasDoublet follow the network singularity types (A502 nts/ntd): terms of an
// absent type are not evaluated and their index lists stay empty.
struct LocalSingularityTerms {
    std::vector<int> sourceId;
    Eigen::VectorXd sourceValue;
    std::vector<int> doubletId;
    Eigen::Matrix<double, 4, Eigen::Dynamic> doubletValue;
    bool valid{false};
};

LocalSingularityTerms buildPanelPointLocalTerms(
    const Network& network, const Compressibility& flow, int panelI,
    int panelJ, const Eigen::Vector3d& fieldPoint,
    const NetworkDoubletParameters& doubletParameters,
    bool hasSource = true, bool hasDoublet = true);

LocalSingularityTerms buildPanelCenterLocalTerms(
    const Network& network, const Compressibility& flow, int panelI,
    int panelJ, const NetworkDoubletParameters& doubletParameters,
    bool hasSource = true, bool hasDoublet = true);

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

// Port of the bcopt() right-hand-side dispatch.  reference is the A502
// zmctr used by nropt=3 and 8 (zero when the condition is not tied to a
// region surface).  nropt=0 produces no right-hand side; 1 is zero;
// 4 (Fortran-coded cbet) and 9 (input value) are not ported.
struct BoundaryRightHandSide {
    double value{0.0};
    bool supported{false};
    bool hasValue{false};
};

BoundaryRightHandSide buildBoundaryRightHandSide(
    int nropt, const Compressibility& flow,
    const Eigen::Vector3d& normal,
    const Eigen::Vector3d& controlPoint,
    const Eigen::Vector3d& reference = Eigen::Vector3d::Zero());

// Result of the bcopt() left-hand-side dispatch (A502 nlopt).  supported is
// false for options that need input coefficients (nlopt=1), a Fortran-coded
// routine (10), material data (18..20), or are not ported (21..25).
// hasCondition is false only for nlopt=0, which means "no boundary
// condition".  Velocity options (11..14) use the default $compressibility
// axis-angle constants al1=al2=0, enx1=enx2=1; the corresponding input
// options are not parsed yet.
struct LeftHandBoundaryOption {
    BoundaryConditionCoefficients coefficients;
    bool supported{false};
    bool hasCondition{false};
};

LeftHandBoundaryOption buildLeftHandBoundaryOption(
    int nlopt, const Compressibility& flow,
    const Eigen::Vector3d& normal);

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
