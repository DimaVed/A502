#pragma once

// Matching.h
// ----------
// First, simple network-edge matching: detect full edges that coincide
// geometrically (same or reversed point order) and equate the corresponding
// doublet parameters.  This covers the common wing/wing or wing/wake full-
// edge abutments; partial overlaps and multi-network intersection points
// still need the general abutment machinery from panair.f90.

#define EIGEN_NO_DEBUG

#include <vector>

#include <Eigen/Dense>

#include "GlobalInfluence.h"

namespace a502 {

struct Network;

struct NetworkEdge {
    int network{-1};
    int edge{0}; // 1..4, as in A502
    std::vector<Eigen::Vector3d> points;
};

std::vector<NetworkEdge> buildNetworkEdges(
    const std::vector<Network>& networks);

struct EdgeMatch {
    int edgeA{0}; // index into the NetworkEdge list
    int edgeB{0};
    bool reversed{false};
};

std::vector<EdgeMatch> findFullEdgeMatches(
    const std::vector<NetworkEdge>& edges, double tolerance);

struct MatchingRows {
    int nRows{0};
    Eigen::MatrixXd doublet; // nRows x totalDoublet
    Eigen::VectorXd rhs;
    bool valid{false};
};

MatchingRows buildEdgeDoubletMatchingRows(
    const std::vector<Network>& networks,
    const GlobalSingularityLayout& layout,
    const std::vector<EdgeMatch>& matches);

} // namespace a502
