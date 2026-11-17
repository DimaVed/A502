// GlobalBoundarySystem.cpp
// ------------------------
// First global boundary-condition system assembly.

#include "GlobalBoundarySystem.h"

#include "AeroModel.h"
#include "BoundaryConditions.h"
#include "ControlPoints.h"
#include "Matching.h"
#include "NetworkInfluence.h"

#include <cmath>

namespace a502 {

GlobalBoundarySystem buildGlobalNormalFluxSystem(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout, double cu, double cl,
    double matchingTolerance, bool includeEdgeAndCorner) {
    GlobalBoundarySystem out;
    if (networks.empty() || layout.network.size() != networks.size())
        return out;

    std::vector<BoundaryRow> bcRows;
    std::vector<double> bcRhs;

    for (std::size_t k = 0; k < networks.size(); ++k) {
        const Network& net = networks[k];
        const NetworkSingularityLayout& netLayout = layout.network[k];
        const std::vector<ControlPoint> controlPoints =
            buildNetworkControlPoints(net);

        for (const ControlPoint& cp : controlPoints) {
            if (!includeEdgeAndCorner &&
                cp.kind != ControlPointKind::PanelCenter)
                continue;

            const GlobalInfluence global =
                evaluateGlobalNearFieldInfluence(
                    networks, flow, cp.position, layout);
            if (!global.valid)
                return GlobalBoundarySystem{};

            const LocalSingularityTerms local =
                buildPanelPointLocalTerms(
                    net, flow, cp.panelI, cp.panelJ, cp.position,
                    netLayout.doublet);
            if (!local.valid)
                return GlobalBoundarySystem{};

            const BoundaryConditionCoefficients bc =
                standardNormalMassFlux(cu, cl, cp.normal, flow.compd);
            BoundaryRow row = buildBoundaryRow(
                global, cp.normal, flow, bc, local);
            if (!row.valid)
                return GlobalBoundarySystem{};

            bcRows.push_back(std::move(row));
            bcRhs.push_back(bc.bet);
        }
    }

    if (bcRows.empty())
        return out;

    const std::vector<NetworkEdge> edges = buildNetworkEdges(networks);
    const std::vector<EdgeMatch> matches =
        findFullEdgeMatches(edges, matchingTolerance);
    MatchingRows matching;
    if (!matches.empty()) {
        matching = buildEdgeDoubletMatchingRows(networks, layout, matches);
        if (!matching.valid)
            return GlobalBoundarySystem{};
    }

    out.nBcRows = static_cast<int>(bcRows.size());
    out.nMatchingRows = matching.valid ? matching.nRows : 0;
    out.nRows = out.nBcRows + out.nMatchingRows;
    out.nActiveSource = layout.activeSource;
    out.nActiveDoublet = layout.activeDoublet;
    out.nColumns = out.nActiveSource + out.nActiveDoublet;
    out.matrix = Eigen::MatrixXd::Zero(out.nRows, out.nColumns);
    out.rhs = Eigen::VectorXd::Zero(out.nRows);

    std::vector<int> sourceMap(
        static_cast<std::size_t>(layout.totalSource), -1);
    std::vector<int> doubletMap(
        static_cast<std::size_t>(layout.totalDoublet), -1);
    for (std::size_t k = 0; k < networks.size(); ++k) {
        const NetworkSingularityLayout& netLayout = layout.network[k];
        if (netLayout.hasSource) {
            for (int i = 0; i < netLayout.source.count; ++i)
                sourceMap[static_cast<std::size_t>(
                    layout.sourceOffset[k] + i)] =
                    layout.activeSourceOffset[k] + i;
        }
        if (netLayout.hasDoublet) {
            for (int i = 0; i < netLayout.doublet.count; ++i)
                doubletMap[static_cast<std::size_t>(
                    layout.doubletOffset[k] + i)] =
                    layout.activeDoubletOffset[k] + i;
        }
    }

    for (int r = 0; r < out.nBcRows; ++r) {
        const BoundaryRow& row = bcRows[static_cast<std::size_t>(r)];
        if (row.source.size() != layout.totalSource ||
            row.doublet.size() != layout.totalDoublet)
            return GlobalBoundarySystem{};

        for (int j = 0; j < layout.totalSource; ++j) {
            const int col = sourceMap[static_cast<std::size_t>(j)];
            if (col >= 0)
                out.matrix(r, col) = row.source(j);
        }
        for (int j = 0; j < layout.totalDoublet; ++j) {
            const int col = doubletMap[static_cast<std::size_t>(j)];
            if (col >= 0)
                out.matrix(r, out.nActiveSource + col) = row.doublet(j);
        }
        out.rhs(r) = bcRhs[static_cast<std::size_t>(r)];
    }

    for (int r = 0; r < out.nMatchingRows; ++r) {
        const int outRow = out.nBcRows + r;
        for (int j = 0; j < layout.totalDoublet; ++j) {
            const double c = matching.doublet(r, j);
            if (c == 0.0)
                continue;
            const int col = doubletMap[static_cast<std::size_t>(j)];
            if (col < 0)
                return GlobalBoundarySystem{};
            out.matrix(outRow, out.nActiveSource + col) += c;
        }
    }

    out.valid = true;
    return out;
}

GlobalSystemSolution solveGlobalBoundarySystem(
    const GlobalBoundarySystem& system) {
    GlobalSystemSolution out;
    if (!system.valid || system.matrix.rows() == 0 ||
        system.matrix.rows() != system.matrix.cols() ||
        system.rhs.size() != system.matrix.rows())
        return out;

    Eigen::FullPivLU<Eigen::MatrixXd> lu(system.matrix);
    if (!lu.isInvertible())
        return out;

    out.unknowns = lu.solve(system.rhs);
    // A couple of iterative-refinement steps keep the residual small for
    // the near-surface near-field matrices.
    for (int it = 0; it < 3; ++it) {
        const Eigen::VectorXd residual =
            system.rhs - system.matrix * out.unknowns;
        if (residual.norm() <= 1.0e-13 * system.rhs.norm())
            break;
        out.unknowns += lu.solve(residual);
    }
    out.residualNorm =
        (system.matrix * out.unknowns - system.rhs).norm();
    out.valid = true;
    return out;
}

} // namespace a502
