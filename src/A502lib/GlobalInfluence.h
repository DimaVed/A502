#pragma once

// GlobalInfluence.h
// -----------------
// Multi-network singularity layout and raw influence rows for the full
// model.  Each network keeps its local parameter ordering; offsets provide
// the global column indices used by the eventual AIC assembly.

#define EIGEN_NO_DEBUG

#include <vector>

#include <Eigen/Dense>

#include "NetworkInfluence.h"

namespace a502 {

struct Network; // defined in AeroModel.h

struct GlobalSingularityLayout {
    std::vector<NetworkSingularityLayout> network;
    std::vector<int> sourceOffset;
    std::vector<int> doubletOffset;
    int totalSource{0};
    int totalDoublet{0};
    // Columns of the active unknown vector, after dropping singularity
    // types that are not present for a network (kt mapping).  A value of -1
    // means that the corresponding type is inactive for that network.
    std::vector<int> activeSourceOffset;
    std::vector<int> activeDoubletOffset;
    int activeSource{0};
    int activeDoublet{0};
};

GlobalSingularityLayout buildGlobalSingularityLayout(
    const std::vector<Network>& networks);

struct GlobalInfluence {
    // rows: potential, global vx,vy,vz
    Eigen::Matrix<double, 4, Eigen::Dynamic> source;
    Eigen::Matrix<double, 4, Eigen::Dynamic> doublet;
    bool valid{false};
};

GlobalInfluence evaluateGlobalNearFieldInfluence(
    const std::vector<Network>& networks, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const GlobalSingularityLayout& layout,
    const SymmetryOptions& symmetry = {});

GlobalInfluence evaluateGlobalPicInfluence(
    const std::vector<Network>& networks, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const GlobalSingularityLayout& layout,
    const SymmetryOptions& symmetry = {});

struct GlobalAic {
    int nControlPoints{0};
    int nSource{0};
    int nDoublet{0};
    Eigen::MatrixXd source;
    Eigen::MatrixXd doublet;
    bool valid{false};
};

GlobalAic buildGlobalNearFieldAic(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout,
    const std::vector<Eigen::Vector3d>& fieldPoints,
    const SymmetryOptions& symmetry = {});

GlobalAic buildGlobalPicAic(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout,
    const std::vector<Eigen::Vector3d>& fieldPoints,
    const SymmetryOptions& symmetry = {});

} // namespace a502
