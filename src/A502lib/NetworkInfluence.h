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

// Default network type and boundary-condition options from the per-kt
// dispatch in A502 inbc():
//   nts / ntd  singularity types (sourceType / doubletType),
//   ipot       potential/velocity flag,
//   nlopt1/2   left-hand-side coefficient options of the first/second
//              boundary condition at a control point,
//   nropt1/2   right-hand-side options.
// Deck-specific overrides read from the kt card (thin-wake thickness,
// type-18 source flag/wake matching, type-19 user options) are not parsed
// yet.  Unknown kt codes fall back to the kt=1 defaults; kt=13 returns its
// supersonic settings for any Mach number.
struct NetworkTypeDefaults {
    int sourceType{1};
    int doubletType{12};
    int ipot{2};
    int nlopt1{5};
    int nropt1{3};
    int nlopt2{7};
    int nropt2{2};
};

NetworkTypeDefaults networkTypeDefaults(int kt, double mach = 0.0);

struct NetworkSingularityLayout {
    NetworkSourceParameters source;
    NetworkDoubletParameters doublet;
    NetworkSingularityTypes types;
    bool hasSource{true};
    bool hasDoublet{true};
};

NetworkSingularityLayout buildNetworkSingularityLayout(
    const Network& network);

// Flow symmetry planes from $symmetry (A502 misym/mjsym): +1 symmetric,
// -1 antisymmetric, 0 absent.  misym is the x-z plane (y reflected), mjsym
// the x-y plane (z reflected).  Image panels share the same singularity
// parameters; their influence is added with the A502 sign pattern.
struct SymmetryOptions {
    int misym{0};
    int mjsym{0};
};

struct NetworkInfluence {
    // rows: potential, global vx,vy,vz
    Eigen::Matrix<double, 4, Eigen::Dynamic> source;
    Eigen::Matrix<double, 4, Eigen::Dynamic> doublet;
    bool valid{false};
};

// The network must have at least two points per direction and valid
// coordinates.  The field point may be anywhere; all panels are evaluated
// with the near-field PIC.  With symmetry the image contributions are
// evaluated at the reflected field points and accumulated with the A502
// sign factors.  Networks touching a symmetry plane are not yet given the
// special half-panel treatment.
NetworkInfluence evaluateNetworkNearFieldInfluence(
    const Network& network, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const NetworkSingularityLayout& layout,
    const SymmetryOptions& symmetry = {});

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
    const std::vector<Eigen::Vector3d>& fieldPoints,
    const SymmetryOptions& symmetry = {});

} // namespace a502
