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
    double matchingTolerance, bool includeEdgeAndCorner = false);

// Diagonal row/column equilibration of the global matrix:
//   A_scaled = diag(row) * A * diag(column),
//   b_scaled = diag(row) * b,
//   x        = diag(column) * y.
// Rows and columns whose largest entry is below tolerance times the largest
// entry of the matrix keep a unit scale, so zero rows/columns stay visible.
struct GlobalSystemScaling {
    Eigen::VectorXd row;
    Eigen::VectorXd column;
    bool valid{false};
};

GlobalSystemScaling buildGlobalSystemScaling(
    const GlobalBoundarySystem& system, double tolerance = 1.0e-12);

struct GlobalSystemDiagnostics {
    int nRows{0};
    int nColumns{0};
    int rank{0};
    int dependentRows{0}; // nRows - rank
    int zeroRows{0};
    int zeroColumns{0};
    double maxAbs{0.0};
    // 2-norm condition number estimate; infinity when rank < min(rows, cols).
    double conditionNumber{0.0};
    bool valid{false};
};

// Rank, conditioning and zero-row/column diagnostics.  Uses a thin SVD and is
// intended for tests and failure analysis, not for the inner solve loop.
GlobalSystemDiagnostics analyzeGlobalBoundarySystem(
    const GlobalBoundarySystem& system, double tolerance = 1.0e-12);

struct GlobalSystemSolution {
    Eigen::VectorXd unknowns; // active source block, then active doublet
    double residualNorm{0.0};
    double scaledResidualNorm{0.0};
    int rank{0};
    int nullity{0};
    // True when a singular square system was solved with the minimum-norm
    // gauge instead of a plain LU solve.
    bool gauged{false};
    bool valid{false};
};

struct GlobalSolveOptions {
    // Closed doublet-only configurations have a constant-doublet null mode:
    // the induced velocity is unchanged when the same constant is added to
    // every doublet parameter.  With this flag a singular square system is
    // solved by a complete orthogonal decomposition, which returns the
    // minimum-norm least-squares solution in the equilibrated unknown space.
    // This is a numerical gauge, not an A502 boundary condition; the
    // resulting surface velocities/Cp do not depend on the choice of gauge.
    // Inconsistent systems are still rejected.
    bool allowRankDeficient{false};
};

// Direct solve for a square system.  The matrix is row/column equilibrated
// before the factorization; the return vector is in the original unknown
// space.  Under/over-determined systems are left for the later
// network-matching/unknown-elimination stage.
GlobalSystemSolution solveGlobalBoundarySystem(
    const GlobalBoundarySystem& system,
    const GlobalSolveOptions& options = {});

struct GlobalParameterValues {
    Eigen::VectorXd source;  // totalSource, zero for inactive networks
    Eigen::VectorXd doublet; // totalDoublet
    bool valid{false};
};

// Scatter the active unknown vector back into full per-network source and
// doublet parameter vectors.
GlobalParameterValues expandGlobalSolution(
    const GlobalSystemSolution& solution,
    const GlobalSingularityLayout& layout);

} // namespace a502
