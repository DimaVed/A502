// Matching.cpp
// ------------
// Full-edge coincidence detection and doublet-parameter matching using
// equivalence classes.

#include "Matching.h"

#include "AeroModel.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
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

std::vector<EdgeMatch> findEdgeMatches(
    const std::vector<NetworkEdge>& edges, double tolerance) {
    std::vector<EdgeMatch> out;
    const double tol = std::max(tolerance, 0.0);
    for (std::size_t i = 0; i < edges.size(); ++i) {
        for (std::size_t j = i + 1; j < edges.size(); ++j) {
            const NetworkEdge& a = edges[i];
            const NetworkEdge& b = edges[j];
            const Eigen::Index na =
                static_cast<Eigen::Index>(a.points.size());
            const Eigen::Index nb =
                static_cast<Eigen::Index>(b.points.size());
            if (na < 2 || nb < 2)
                continue;

            int bestCount = 0;
            int bestEndA = 0;
            int bestEndB = 0;
            bool bestReversed = false;
            for (int orientation = 0; orientation < 2; ++orientation) {
                const bool reversed = orientation == 1;
                // dp[j] = length of the common chain ending at a[i] and at
                // the (possibly reversed) b[j].
                std::vector<int> previous(static_cast<std::size_t>(nb), 0);
                std::vector<int> current(static_cast<std::size_t>(nb), 0);
                for (Eigen::Index ai = 0; ai < na; ++ai) {
                    std::fill(current.begin(), current.end(), 0);
                    for (Eigen::Index bj = 0; bj < nb; ++bj) {
                        const Eigen::Index brow = reversed ? nb - 1 - bj : bj;
                        if ((a.points[static_cast<std::size_t>(ai)] -
                             b.points[static_cast<std::size_t>(brow)])
                                .norm() > tol)
                            continue;
                        const int diag =
                            (ai > 0 && bj > 0)
                                ? previous[static_cast<std::size_t>(bj - 1)]
                                : 0;
                        current[static_cast<std::size_t>(bj)] = diag + 1;
                        if (diag + 1 > bestCount) {
                            bestCount = diag + 1;
                            bestEndA = static_cast<int>(ai);
                            bestEndB = static_cast<int>(bj);
                            bestReversed = reversed;
                        }
                    }
                    previous.swap(current);
                }
            }
            if (bestCount < 2)
                continue;

            EdgeMatch match;
            match.edgeA = static_cast<int>(i);
            match.edgeB = static_cast<int>(j);
            match.reversed = bestReversed;
            match.count = bestCount;
            match.startA = bestEndA - bestCount + 1;
            match.startB = bestReversed
                               ? static_cast<int>(nb) - 1 - bestEndB
                               : bestEndB - bestCount + 1;
            match.partial =
                !(match.count == static_cast<int>(na) &&
                  match.count == static_cast<int>(nb));
            out.push_back(match);
        }
    }
    return out;
}

std::vector<EdgeMatch> findFullEdgeMatches(
    const std::vector<NetworkEdge>& edges, double tolerance) {
    std::vector<EdgeMatch> out;
    for (const EdgeMatch& match : findEdgeMatches(edges, tolerance))
        if (!match.partial)
            out.push_back(match);
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

    // Equivalence classes of matched doublet parameters, as in A502's
    // abutment processing.  A constraint joining two parameters that are
    // already in the same class is implied by the class representative
    // rows and is not emitted again (otherwise corners shared by three or
    // more networks would add dependent rows).
    std::vector<int> parent(
        static_cast<std::size_t>(std::max(layout.totalDoublet, 0)));
    std::iota(parent.begin(), parent.end(), 0);
    auto findRoot = [&](int id) {
        int root = id;
        while (parent[static_cast<std::size_t>(root)] != root)
            root = parent[static_cast<std::size_t>(root)];
        while (parent[static_cast<std::size_t>(id)] != id) {
            const int next = parent[static_cast<std::size_t>(id)];
            parent[static_cast<std::size_t>(id)] = root;
            id = next;
        }
        return root;
    };

    auto addRow = [&](int globalA, int globalB) {
        if (globalA < 0 || globalB < 0 ||
            globalA >= layout.totalDoublet ||
            globalB >= layout.totalDoublet)
            return false;
        const int rootA = findRoot(globalA);
        const int rootB = findRoot(globalB);
        if (rootA == rootB)
            return true; // already implied by the equivalence class
        Eigen::VectorXd row =
            Eigen::VectorXd::Zero(layout.totalDoublet);
        row(globalA) += 1.0;
        row(globalB) -= 1.0;
        rows.push_back(std::move(row));
        parent[static_cast<std::size_t>(rootB)] = rootA;
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
            b.network >= static_cast<int>(layout.network.size()))
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

        // EdgeMatch objects created without start/count describe a complete
        // edge coincidence (legacy full-edge form).
        int startIndexA = match.startA;
        int startIndexB = match.startB;
        int count = match.count;
        if (count <= 0) {
            if (a.points.size() != b.points.size())
                return MatchingRows{};
            startIndexA = 0;
            startIndexB = 0;
            count = static_cast<int>(a.points.size());
        }
        if (count < 2 || startIndexA < 0 || startIndexB < 0 ||
            startIndexA + count > static_cast<int>(a.points.size()) ||
            startIndexB + count > static_cast<int>(b.points.size()))
            return MatchingRows{};

        const int nIntervals = count - 1;
        const int baseA =
            da.edgeMidBase[static_cast<std::size_t>(a.edge - 1)];
        const int baseB =
            db.edgeMidBase[static_cast<std::size_t>(b.edge - 1)];
        for (int k = 0; k < nIntervals; ++k) {
            const int indexA = startIndexA + k;
            const int indexB = match.reversed
                                   ? startIndexB + nIntervals - 1 - k
                                   : startIndexB + k;
            if (!addRow(offsetA + baseA + indexA,
                        offsetB + baseB + indexB))
                return MatchingRows{};
        }

        // Network corners at the shared chain endpoints.  Interior edge
        // points of a partial match are not doublet parameters.
        auto cornerAt = [&](const NetworkEdge& edge,
                            const NetworkDoubletParameters& parameters,
                            int pointIndex) -> int {
            if (edge.edge < 1 || edge.edge > 4)
                return -1;
            const int last =
                static_cast<int>(edge.points.size()) - 1;
            if (pointIndex == 0)
                return parameters.cornerId[static_cast<std::size_t>(
                    startCorner[static_cast<std::size_t>(edge.edge - 1)])];
            if (pointIndex == last)
                return parameters.cornerId[static_cast<std::size_t>(
                    endCorner[static_cast<std::size_t>(edge.edge - 1)])];
            return -1;
        };
        auto cornerPair = [&](int indexA, int indexB) {
            const int cornerA = cornerAt(a, da, indexA);
            const int cornerB = cornerAt(b, db, indexB);
            if (cornerA < 0 || cornerB < 0)
                return true;
            return addRow(offsetA + cornerA, offsetB + cornerB);
        };
        const int startBIndex =
            match.reversed ? startIndexB + count - 1 : startIndexB;
        const int endBIndex =
            match.reversed ? startIndexB : startIndexB + count - 1;
        if (!cornerPair(startIndexA, startBIndex))
            return MatchingRows{};
        if (!cornerPair(startIndexA + count - 1, endBIndex))
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
