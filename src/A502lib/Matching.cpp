// Matching.cpp
// ------------
// Simple full-edge coincidence detection and doublet-parameter matching.

#include "Matching.h"

#include "AeroModel.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <set>
#include <utility>

namespace a502 {

namespace {

Eigen::Vector3d networkPoint(const Network& net, int m, int n) {
    // 1-based m,n.
    const std::size_t g =
        static_cast<std::size_t>(n - 1) * static_cast<std::size_t>(net.nm) +
        static_cast<std::size_t>(m - 1);
    return Eigen::Vector3d(net.coordinates[g * 3],
                           net.coordinates[g * 3 + 1],
                           net.coordinates[g * 3 + 2]);
}

} // namespace

std::vector<NetworkEdge> buildNetworkEdges(
    const std::vector<Network>& networks) {
    std::vector<NetworkEdge> out;
    for (std::size_t k = 0; k < networks.size(); ++k) {
        const Network& net = networks[k];
        const int nm = net.nm;
        const int nn = net.nn;
        if (nm < 2 || nn < 2)
            continue;
        const std::size_t nPoints =
            static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
        if (net.coordinates.size() < nPoints * 3)
            continue;

        for (int edge = 1; edge <= 4; ++edge) {
            NetworkEdge e;
            e.network = static_cast<int>(k);
            e.edge = edge;
            if (edge == 1) {
                for (int n = 1; n <= nn; ++n)
                    e.points.push_back(networkPoint(net, 1, n));
            } else if (edge == 2) {
                for (int m = 1; m <= nm; ++m)
                    e.points.push_back(networkPoint(net, m, nn));
            } else if (edge == 3) {
                for (int n = 1; n <= nn; ++n)
                    e.points.push_back(networkPoint(net, nm, n));
            } else {
                for (int m = 1; m <= nm; ++m)
                    e.points.push_back(networkPoint(net, m, 1));
            }
            out.push_back(std::move(e));
        }
    }
    return out;
}

std::vector<EdgeMatch> findFullEdgeMatches(
    const std::vector<NetworkEdge>& edges, double tolerance) {
    std::vector<EdgeMatch> out;
    const double tol = std::max(tolerance, 0.0);
    for (std::size_t i = 0; i < edges.size(); ++i) {
        for (std::size_t j = i + 1; j < edges.size(); ++j) {
            const NetworkEdge& a = edges[i];
            const NetworkEdge& b = edges[j];
            if (a.points.size() != b.points.size() ||
                a.points.empty())
                continue;

            bool same = true;
            bool reversed = true;
            const std::size_t n = a.points.size();
            for (std::size_t k = 0; k < n; ++k) {
                if ((a.points[k] - b.points[k]).norm() > tol)
                    same = false;
                if ((a.points[k] - b.points[n - 1 - k]).norm() > tol)
                    reversed = false;
                if (!same && !reversed)
                    break;
            }
            if (!same && !reversed)
                continue;

            EdgeMatch match;
            match.edgeA = static_cast<int>(i);
            match.edgeB = static_cast<int>(j);
            match.reversed = !same && reversed;
            out.push_back(match);
        }
    }
    return out;
}

MatchingRows buildEdgeDoubletMatchingRows(
    const std::vector<Network>& networks,
    const GlobalSingularityLayout& layout,
    const std::vector<EdgeMatch>& matches) {
    MatchingRows out;
    const std::vector<NetworkEdge> edges = buildNetworkEdges(networks);
    if (edges.empty() || layout.network.size() != networks.size())
        return out;

    const std::array<int, 4> startCorner{{0, 1, 3, 0}};
    const std::array<int, 4> endCorner{{1, 2, 2, 3}};
    std::vector<Eigen::VectorXd> rows;
    std::set<std::pair<int, int>> cornerPairs;

    auto addRow = [&](int globalA, int globalB) {
        if (globalA < 0 || globalB < 0 ||
            globalA >= layout.totalDoublet ||
            globalB >= layout.totalDoublet)
            return false;
        Eigen::VectorXd row =
            Eigen::VectorXd::Zero(layout.totalDoublet);
        row(globalA) += 1.0;
        row(globalB) -= 1.0;
        rows.push_back(std::move(row));
        return true;
    };

    for (const EdgeMatch& match : matches) {
        if (match.edgeA < 0 || match.edgeB < 0 ||
            match.edgeA >= static_cast<int>(edges.size()) ||
            match.edgeB >= static_cast<int>(edges.size()))
            return MatchingRows{};
        const NetworkEdge& a = edges[static_cast<std::size_t>(match.edgeA)];
        const NetworkEdge& b = edges[static_cast<std::size_t>(match.edgeB)];
        if (a.network < 0 || b.network < 0 ||
            a.network >= static_cast<int>(layout.network.size()) ||
            b.network >= static_cast<int>(layout.network.size()) ||
            a.points.size() != b.points.size() || a.points.size() < 2)
            return MatchingRows{};

        // A502 abtecd() removes source-alone networks (ntd = 0) from the
        // abutment search when xsrcab is set, so their edges never produce
        // doublet matching equations.
        if (!layout.network[static_cast<std::size_t>(a.network)].hasDoublet ||
            !layout.network[static_cast<std::size_t>(b.network)].hasDoublet)
            continue;

        const NetworkDoubletParameters& da =
            layout.network[static_cast<std::size_t>(a.network)].doublet;
        const NetworkDoubletParameters& db =
            layout.network[static_cast<std::size_t>(b.network)].doublet;
        const int offsetA =
            layout.doubletOffset[static_cast<std::size_t>(a.network)];
        const int offsetB =
            layout.doubletOffset[static_cast<std::size_t>(b.network)];

        const int nIntervals = static_cast<int>(a.points.size()) - 1;
        const int baseA =
            da.edgeMidBase[static_cast<std::size_t>(a.edge - 1)];
        const int baseB =
            db.edgeMidBase[static_cast<std::size_t>(b.edge - 1)];
        for (int k = 0; k < nIntervals; ++k) {
            const int kb = match.reversed ? (nIntervals - 1 - k) : k;
            if (!addRow(offsetA + baseA + k, offsetB + baseB + kb))
                return MatchingRows{};
        }

        // Endpoint corners.
        const int startA =
            da.cornerId[static_cast<std::size_t>(
                startCorner[static_cast<std::size_t>(a.edge - 1)])];
        const int endA =
            da.cornerId[static_cast<std::size_t>(
                endCorner[static_cast<std::size_t>(a.edge - 1)])];
        const int startB =
            db.cornerId[static_cast<std::size_t>(
                startCorner[static_cast<std::size_t>(b.edge - 1)])];
        const int endB =
            db.cornerId[static_cast<std::size_t>(
                endCorner[static_cast<std::size_t>(b.edge - 1)])];

        const int a0 = offsetA + startA;
        const int a1 = offsetA + endA;
        const int b0 = offsetB + startB;
        const int b1 = offsetB + endB;
        const std::pair<int, int> p0 =
            match.reversed ? std::make_pair(a0, b1)
                           : std::make_pair(a0, b0);
        const std::pair<int, int> p1 =
            match.reversed ? std::make_pair(a1, b0)
                           : std::make_pair(a1, b1);
        if (cornerPairs.insert(p0).second && !addRow(p0.first, p0.second))
            return MatchingRows{};
        if (cornerPairs.insert(p1).second && !addRow(p1.first, p1.second))
            return MatchingRows{};
    }

    // Matches may all be ineligible (e.g. source-alone networks); this is a
    // valid empty matching set, not an assembly failure.
    out.nRows = static_cast<int>(rows.size());
    out.doublet = Eigen::MatrixXd::Zero(out.nRows, layout.totalDoublet);
    out.rhs = Eigen::VectorXd::Zero(out.nRows);
    for (int i = 0; i < out.nRows; ++i)
        out.doublet.row(i) = rows[static_cast<std::size_t>(i)].transpose();
    out.valid = true;
    return out;
}

} // namespace a502
