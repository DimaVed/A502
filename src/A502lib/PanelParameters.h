#pragma once

// PanelParameters.h
// -----------------
// Utilities for the bilinear panel surface used by nrpthp()/surfit():
//   * closest (s,t) parameters of a point to the four-corner bilinear patch;
//   * near-plane transformation and skewness parameters c1,c2.
//
// Canonical points map to
//   P1(1,1)  P2(-1,1)  P3(-1,-1)  P4(1,-1)
//   P5(0,1)  P6(-1,0)  P7(0,-1)   P8(1,0)  P9(0,0).

#define EIGEN_NO_DEBUG

#include <array>

#include <Eigen/Dense>

namespace a502 {

struct PanelParameters {
    double s{0.0};
    double t{0.0};
    bool   valid{false};
};

Eigen::Vector3d evalBilinearPanel(
    const std::array<Eigen::Vector3d, 4>& corners, double s, double t);

PanelParameters closestPanelParameters(
    const std::array<Eigen::Vector3d, 4>& corners, const Eigen::Vector3d& p);

struct NearPlane {
    Eigen::Matrix3d toNear{Eigen::Matrix3d::Identity()};
    // Third column of surfit's art matrix: v1 x v2 / sqrt(|v1 x v2|).
    Eigen::Vector3d normal{Eigen::Vector3d::UnitZ()};
    bool            valid{false};
};

// surfit(): near-plane coordinate transform relative to the panel center.
NearPlane buildNearPlane(const std::array<Eigen::Vector3d, 4>& corners);

// ccaln()/panmom(): area moments of the panel polygon projected onto the
// local x-y plane:
//   c(i,j) = integral x^(i-1) y^(j-1) dA,  i,j = 1..degree+1.
// ics is the skipped (collapsed) edge, 0 for a quadrilateral.  The vertex
// order must be counterclockwise for a positive area.
Eigen::MatrixXd panelAreaMoments(
    const std::array<Eigen::Vector3d, 4>& vertices, int nside, int ics,
    int degree);

// qcof(): near-plane quadratic doublet distribution used by the far-field
// expansion.  The rows are the coefficients of
//   [1, x, y, x^2/2, x*y, y^2/2]
// in the inner-panel frame, and the columns are the canonical doublet
// values at P1..P9.
struct NearPlaneQuadratic {
    Eigen::Matrix<double, 6, 9> coefficients{
        Eigen::Matrix<double, 6, 9>::Zero()};
    bool valid{false};
};

NearPlaneQuadratic buildNearPlaneQuadratic(
    const std::array<Eigen::Vector3d, 9>& cp,
    const Eigen::Matrix3d& innerToLocal);

} // namespace a502
