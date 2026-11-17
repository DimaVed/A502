#pragma once

// PanelGeometry.h
// ---------------
// Geometric precursor of the higher-order panel method.
//
// A $points network is an nm x nn grid of corner points.  Every cell
// (panel) is described by the nine canonical points used by PAN AIR:
//
//      P1----P5----P2
//      |  \   |   /  |
//      P8----P9----P6
//      |  /   |   \  |
//      P4----P7----P3
//
// and is split into eight planar triangular subpanels (gensub in panair.f90):
//
//   1: P1 P5 P8     5: P5 P9 P8
//   2: P2 P6 P5     6: P6 P9 P5
//   3: P3 P7 P6     7: P7 P9 P6
//   4: P4 P8 P7     8: P8 P9 P7
//
// This module currently stores the global geometry, unit normals, areas and
// an orthonormal local frame for each subpanel.  The compressibility-scaled
// local frame of refloc() will be added together with the influence
// coefficient module; the frames here are its geometric counterpart and are
// already useful for visualization and area/normal checks.

#define EIGEN_NO_DEBUG

#include <array>
#include <cstddef>
#include <vector>

#include <Eigen/Dense>

namespace a502 {

struct SubpanelGeometry {
    std::array<Eigen::Vector3d, 3> vertex{};
    Eigen::Vector3d                centroid{Eigen::Vector3d::Zero()};
    Eigen::Vector3d                normal{Eigen::Vector3d::UnitZ()};
    double                         area{0.0};

    // Orthonormal frame: local = toLocal * (global - origin),
    // global = origin + toGlobal * local.
    Eigen::Vector3d origin{Eigen::Vector3d::Zero()};
    Eigen::Matrix3d toLocal{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d toGlobal{Eigen::Matrix3d::Identity()};

    // Compressibility-scaled frame from refloc() (set by
    // applyCompressibleFrames() once the freestream is known).
    Eigen::Matrix3d scaledToLocal{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d scaledToGlobal{Eigen::Matrix3d::Identity()};
    double          jacobian{1.0};
    int             inclination{1}; // +1 subinclined, -1 superinclined
    bool            scaledValid{false};

    bool degenerate{false};
};

struct PanelGeometry {
    int                                  i{0}; // row index in the network
    int                                  j{0}; // column index in the network
    std::array<Eigen::Vector3d, 9>       cp{}; // canonical points P1..P9
    Eigen::Vector3d                      center{Eigen::Vector3d::Zero()};
    Eigen::Vector3d                      normal{Eigen::Vector3d::UnitZ()};
    double                               area{0.0};

    std::array<SubpanelGeometry, 8>      subpanels{};

    int  collapsedEdge{0}; // 0 = quadrilateral; 1..4 = collapsed edge
    bool degenerate{false};
};

// Builds the geometry of panel (i,j) of an nm x nn network.  Coordinates are
// stored as x,y,z triplets in A502 order: point(m,n) = coordinates[(n*nm+m)*3].
// Returns a degenerate panel when indices are out of range or the cell has no
// usable plane.
PanelGeometry buildPanelGeometry(const std::vector<double>& coordinates,
                                 int nm, int nn, int i, int j);

} // namespace a502
