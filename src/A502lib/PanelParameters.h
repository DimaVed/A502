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
    bool            valid{false};
};

// surfit(): near-plane coordinate transform relative to the panel center.
NearPlane buildNearPlane(const std::array<Eigen::Vector3d, 4>& corners);

} // namespace a502
