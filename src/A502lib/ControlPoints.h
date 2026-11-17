#pragma once

// ControlPoints.h
// ---------------
// First port of the A502 control point generator (contrl / gcpcal / zcadj).
//
// This version covers the control points of a single network:
//   * one control point near each panel center;
//   * one near the midpoint of every panel edge that lies on a network edge;
//   * one near every network corner.
//
// Control points are withdrawn from subpanel vertices/edges (xxadj/zcadj in
// panair.f90) to avoid the infinite self-induced velocity on the singularity
// surface.  Abutment matching points and extra control points at abutment
// intersections will be added together with the network matching module.
//
// A network edge control point belongs to a specific panel and edge so the
// future AIC code can associate it with the correct boundary condition and
// doublet parameter.

#define EIGEN_NO_DEBUG

#include <Eigen/Dense>

#include <vector>

namespace a502 {

struct Network; // defined in AeroModel.h

enum class ControlPointKind {
    PanelCenter,
    EdgeMidpoint,
    Corner
};

struct ControlPoint {
    ControlPointKind kind{ControlPointKind::PanelCenter};
    Eigen::Vector3d  position{Eigen::Vector3d::Zero()};
    Eigen::Vector3d  normal{Eigen::Vector3d::UnitZ()};

    int network{0};
    int panelI{0};      // row index of the owning panel
    int panelJ{0};      // column index of the owning panel
    int panelIndex{0};  // i + (nm-1)*j within the network
    int edge{0};        // 1..4 for EdgeMidpoint, 0 otherwise
    int corner{0};      // 1..4 for Corner, 0 otherwise
};

// Builds the control points of one network.  The input coordinates follow the
// A502 column-major order used by AeroModel::Network.
std::vector<ControlPoint> buildNetworkControlPoints(const Network& network);

} // namespace a502
