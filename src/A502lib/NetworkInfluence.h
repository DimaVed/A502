#pragma once

// NetworkInfluence.h
// ------------------
// Single-network influence row: the potential and velocity at a field point
// induced by all source and doublet parameters of one network.
//
// Source parameters are one per panel centre; the source outer spline
// (SourceOuterSpline) maps the 3x3 neighboring parameters to sigma at P1..P4,
// P9.  Doublet parameters use the staggered (nm+1)x(nn+1) layout of gcpcal.
// The panel near-field PIC then maps the canonical values at P1..P9 to the
// field point.  This is the row-assembly step used by the AIC matrix; the
// present implementation evaluates every panel with the near-field kernel.

#define EIGEN_NO_DEBUG

#include <vector>

#include <Eigen/Dense>

#include "Compressibility.h"
#include "OuterSplines.h"

namespace a502 {

struct Network;

struct NetworkSourceParameters {
    int nPanelM{0};
    int nPanelN{0};
    int count{0};
    std::vector<Eigen::Vector3d> location; // panel centres

    int panelId(int panelI, int panelJ) const {
        return panelI + nPanelM * panelJ;
    }
};

struct NetworkSingularityTypes {
    int sourceType{0};  // A502 nts
    int doubletType{0}; // A502 ntd
};

// Network-type default from the contrl/kt dispatch table.  Unknown kt codes
// keep both singularity types active.
NetworkSingularityTypes networkSingularityTypes(int kt);

struct NetworkSingularityLayout {
    NetworkSourceParameters source;
    NetworkDoubletParameters doublet;
    NetworkSingularityTypes types;
    bool hasSource{true};
    bool hasDoublet{true};
};

NetworkSingularityLayout buildNetworkSingularityLayout(
    const Network& network);

struct NetworkInfluence {
    // rows: potential, global vx,vy,vz
    Eigen::Matrix<double, 4, Eigen::Dynamic> source;
    Eigen::Matrix<double, 4, Eigen::Dynamic> doublet;
    bool valid{false};
};

// The network must have at least two points per direction and valid
// coordinates.  The field point may be anywhere; all panels are evaluated
// with the subsonic near-field kernel.
NetworkInfluence evaluateNetworkNearFieldInfluence(
    const Network& network, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const NetworkSingularityLayout& layout);

// Raw AIC rows at a list of field points, four rows per point (potential,
// global vx,vy,vz).  Boundary conditions select/combine these rows later.
struct NetworkAic {
    int nControlPoints{0};
    int nSource{0};
    int nDoublet{0};
    Eigen::MatrixXd source;
    Eigen::MatrixXd doublet;
    bool valid{false};
};

NetworkAic buildNetworkNearFieldAic(
    const Network& network, const Compressibility& flow,
    const NetworkSingularityLayout& layout,
    const std::vector<Eigen::Vector3d>& fieldPoints);

} // namespace a502
