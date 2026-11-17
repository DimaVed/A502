// GlobalBoundarySystem.cpp
// ------------------------
// First global boundary-condition system assembly.

#include "GlobalBoundarySystem.h"

#include "AeroModel.h"
#include "BoundaryConditions.h"
#include "ControlPoints.h"
#include "Matching.h"
#include "NetworkInfluence.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace a502 {

namespace {

bool isUsableGlobalSystem(const GlobalBoundarySystem& system) {
    return system.valid && system.matrix.rows() > 0 &&
           system.matrix.cols() > 0 &&
           system.rhs.size() == system.matrix.rows() &&
           system.matrix.allFinite() && system.rhs.allFinite();
}

} // namespace

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
                    netLayout.doublet, netLayout.hasSource,
                    netLayout.hasDoublet);
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

GlobalSystemScaling buildGlobalSystemScaling(
    const GlobalBoundarySystem& system, double tolerance) {
    GlobalSystemScaling out;
    if (!isUsableGlobalSystem(system))
        return out;

    const Eigen::Index rows = system.matrix.rows();
    const Eigen::Index cols = system.matrix.cols();
    out.row = Eigen::VectorXd::Ones(rows);
    out.column = Eigen::VectorXd::Ones(cols);

    const double tol = std::max(tolerance, 0.0);
    const double maxAbs = system.matrix.cwiseAbs().maxCoeff();
    const double threshold = maxAbs > 0.0 ? tol * maxAbs : 0.0;

    const Eigen::VectorXd rowMax =
        system.matrix.cwiseAbs().rowwise().maxCoeff();
    for (Eigen::Index i = 0; i < rows; ++i) {
        if (rowMax(i) > threshold)
            out.row(i) = 1.0 / rowMax(i);
    }

    const Eigen::MatrixXd rowScaled =
        out.row.asDiagonal() * system.matrix;
    const double scaledMax = rowScaled.cwiseAbs().maxCoeff();
    const double columnThreshold =
        scaledMax > 0.0 ? tol * scaledMax : 0.0;
    const Eigen::VectorXd columnMax =
        rowScaled.cwiseAbs().colwise().maxCoeff();
    for (Eigen::Index j = 0; j < cols; ++j) {
        if (columnMax(j) > columnThreshold)
            out.column(j) = 1.0 / columnMax(j);
    }

    out.valid = true;
    return out;
}

GlobalSystemDiagnostics analyzeGlobalBoundarySystem(
    const GlobalBoundarySystem& system, double tolerance) {
    GlobalSystemDiagnostics out;
    if (!isUsableGlobalSystem(system))
        return out;

    out.nRows = static_cast<int>(system.matrix.rows());
    out.nColumns = static_cast<int>(system.matrix.cols());
    out.maxAbs = system.matrix.cwiseAbs().maxCoeff();

    const double tol = std::max(tolerance, 0.0);
    const double threshold = out.maxAbs > 0.0 ? tol * out.maxAbs : 0.0;

    const Eigen::VectorXd rowMax =
        system.matrix.cwiseAbs().rowwise().maxCoeff();
    const Eigen::VectorXd columnMax =
        system.matrix.cwiseAbs().colwise().maxCoeff();
    for (Eigen::Index i = 0; i < rowMax.size(); ++i)
        if (rowMax(i) <= threshold)
            ++out.zeroRows;
    for (Eigen::Index j = 0; j < columnMax.size(); ++j)
        if (columnMax(j) <= threshold)
            ++out.zeroColumns;

    const Eigen::JacobiSVD<Eigen::MatrixXd> svd(system.matrix);
    const Eigen::VectorXd values = svd.singularValues();
    if (values.size() == 0 || !(values(0) > 0.0)) {
        out.dependentRows = out.nRows;
        out.conditionNumber = std::numeric_limits<double>::infinity();
        out.valid = true;
        return out;
    }

    const double rankThreshold = tol * values(0);
    for (Eigen::Index k = 0; k < values.size(); ++k)
        if (values(k) > rankThreshold)
            ++out.rank;
    out.dependentRows = std::max(0, out.nRows - out.rank);

    const Eigen::Index minDimension =
        std::min(system.matrix.rows(), system.matrix.cols());
    if (out.rank == minDimension && values(minDimension - 1) > 0.0)
        out.conditionNumber = values(0) / values(minDimension - 1);
    else
        out.conditionNumber = std::numeric_limits<double>::infinity();

    out.valid = true;
    return out;
}

GlobalEquationAudit auditGlobalEquationCounts(
    const std::vector<Network>& networks,
    const GlobalSingularityLayout& layout,
    double matchingTolerance) {
    GlobalEquationAudit out;
    if (networks.empty() || layout.network.size() != networks.size())
        return out;

    out.networks.reserve(networks.size());
    for (std::size_t k = 0; k < networks.size(); ++k) {
        const Network& net = networks[k];
        const NetworkSingularityLayout& netLayout = layout.network[k];

        NetworkEquationAudit entry;
        entry.network = static_cast<int>(k);
        entry.kt = net.kt;
        entry.hasSource = netLayout.hasSource;
        entry.hasDoublet = netLayout.hasDoublet;
        entry.nSourceParameters = netLayout.source.count;
        entry.nDoubletParameters = netLayout.doublet.count;
        entry.activeSource =
            netLayout.hasSource ? netLayout.source.count : 0;
        entry.activeDoublet =
            netLayout.hasDoublet ? netLayout.doublet.count : 0;

        const std::vector<ControlPoint> controlPoints =
            buildNetworkControlPoints(net);
        entry.nControlPoints = static_cast<int>(controlPoints.size());
        for (const ControlPoint& cp : controlPoints) {
            switch (cp.kind) {
            case ControlPointKind::PanelCenter:
                ++entry.nPanelCenters;
                break;
            case ControlPointKind::EdgeMidpoint:
                ++entry.nEdgePoints;
                break;
            case ControlPointKind::Corner:
                ++entry.nCornerPoints;
                break;
            }
        }

        out.totalActiveUnknowns +=
            entry.activeSource + entry.activeDoublet;
        out.totalControlPoints += entry.nControlPoints;
        out.totalCenterRows += entry.nPanelCenters;
        out.totalAllPointRows += entry.nControlPoints;
        out.networks.push_back(entry);
    }

    const std::vector<NetworkEdge> edges = buildNetworkEdges(networks);
    const std::vector<EdgeMatch> matches =
        findFullEdgeMatches(edges, matchingTolerance);
    if (!matches.empty()) {
        const MatchingRows rows =
            buildEdgeDoubletMatchingRows(networks, layout, matches);
        if (!rows.valid)
            return GlobalEquationAudit{};
        out.nMatchingRows = rows.nRows;
    }

    out.centerDeficit = out.totalActiveUnknowns -
                        out.totalCenterRows - out.nMatchingRows;
    out.allPointDeficit = out.totalActiveUnknowns -
                          out.totalAllPointRows - out.nMatchingRows;
    out.valid = true;
    return out;
}

GlobalSystemSolution solveGlobalBoundarySystem(
    const GlobalBoundarySystem& system,
    const GlobalSolveOptions& options) {
    GlobalSystemSolution out;
    if (!isUsableGlobalSystem(system))
        return out;
    const bool square = system.matrix.rows() == system.matrix.cols();
    if (!square && !options.allowRankDeficient)
        return out;

    const GlobalSystemScaling scaling =
        buildGlobalSystemScaling(system);
    if (!scaling.valid)
        return out;

    const Eigen::MatrixXd scaledMatrix =
        scaling.row.asDiagonal() * system.matrix *
        scaling.column.asDiagonal();
    const Eigen::VectorXd scaledRhs =
        scaling.row.asDiagonal() * system.rhs;

    Eigen::FullPivLU<Eigen::MatrixXd> lu(scaledMatrix);
    if (square && lu.isInvertible()) {
        out.rank = static_cast<int>(lu.rank());

        // Iterative refinement is done entirely in the equilibrated system:
        // forming the residual in the original units and mapping the
        // correction back through large column scales would amplify
        // round-off.
        Eigen::VectorXd scaledUnknowns = lu.solve(scaledRhs);
        for (int it = 0; it < 3; ++it) {
            const Eigen::VectorXd residual =
                scaledRhs - scaledMatrix * scaledUnknowns;
            if (residual.norm() <= 1.0e-13 * scaledRhs.norm())
                break;
            scaledUnknowns += lu.solve(residual);
        }

        out.unknowns =
            scaling.column.asDiagonal() * scaledUnknowns;
        out.residualNorm =
            (system.matrix * out.unknowns - system.rhs).norm();
        out.scaledResidualNorm =
            (scaledMatrix * scaledUnknowns - scaledRhs).norm();
        out.valid = true;
        return out;
    }

    if (!options.allowRankDeficient)
        return out;

    // Rank-deficient (singular) square system: solve the consistent
    // least-squares problem and use the minimum-norm solution as the gauge.
    // For closed doublet-only networks this fixes the arbitrary additive
    // constant of the doublet potential.
    const Eigen::Index n = scaledMatrix.cols();
    const double scaledMax = scaledMatrix.cwiseAbs().maxCoeff();
    Eigen::CompleteOrthogonalDecomposition<Eigen::MatrixXd> cod;
    cod.setThreshold(1.0e-12 * (scaledMax > 0.0 ? scaledMax : 1.0));
    cod.compute(scaledMatrix);
    const int rank = static_cast<int>(cod.rank());
    if (rank <= 0)
        return out;

    Eigen::VectorXd scaledUnknowns = cod.solve(scaledRhs);
    for (int it = 0; it < 3; ++it) {
        const Eigen::VectorXd residual =
            scaledRhs - scaledMatrix * scaledUnknowns;
        if (residual.norm() <= 1.0e-13 * scaledRhs.norm())
            break;
        scaledUnknowns += cod.solve(residual);
    }

    out.scaledResidualNorm =
        (scaledMatrix * scaledUnknowns - scaledRhs).norm();
    // Reject an inconsistent system instead of silently returning a
    // least-squares compromise.
    const double consistencyScale =
        std::max(1.0, scaledRhs.norm());
    if (!(out.scaledResidualNorm <= 1.0e-8 * consistencyScale))
        return out;

    out.unknowns =
        scaling.column.asDiagonal() * scaledUnknowns;
    out.residualNorm =
        (system.matrix * out.unknowns - system.rhs).norm();
    out.rank = rank;
    out.nullity = static_cast<int>(n) - rank;
    out.gauged = rank < static_cast<int>(n);
    out.valid = true;
    return out;
}

GlobalParameterValues expandGlobalSolution(
    const GlobalSystemSolution& solution,
    const GlobalSingularityLayout& layout) {
    GlobalParameterValues out;
    if (!solution.valid ||
        solution.unknowns.size() !=
            layout.activeSource + layout.activeDoublet)
        return out;

    out.source = Eigen::VectorXd::Zero(layout.totalSource);
    out.doublet = Eigen::VectorXd::Zero(layout.totalDoublet);
    for (std::size_t k = 0; k < layout.network.size(); ++k) {
        const NetworkSingularityLayout& netLayout = layout.network[k];
        if (netLayout.hasSource) {
            const int count = netLayout.source.count;
            out.source.segment(layout.sourceOffset[k], count) =
                solution.unknowns.segment(
                    layout.activeSourceOffset[k], count);
        }
        if (netLayout.hasDoublet) {
            const int count = netLayout.doublet.count;
            out.doublet.segment(layout.doubletOffset[k], count) =
                solution.unknowns.segment(
                    layout.activeSource + layout.activeDoubletOffset[k],
                    count);
        }
    }

    out.valid = true;
    return out;
}

} // namespace a502
