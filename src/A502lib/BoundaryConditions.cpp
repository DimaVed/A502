// BoundaryConditions.cpp
// ----------------------
// Port of the nct=1/nct=4 core of aical().

#include "BoundaryConditions.h"

#include "AeroModel.h"
#include "OuterSplines.h"
#include "PanelGeometry.h"
#include "SubpanelSplines.h"

#include <cmath>

namespace a502 {

namespace {

Eigen::Vector3d compressibilityScale(double bs,
                                     const Eigen::Vector3d& direction,
                                     const Eigen::Vector3d& v) {
    const double c2 = direction.squaredNorm();
    if (c2 <= 1.0e-30)
        return v;
    const double f = (bs - 1.0) * v.dot(direction) / c2;
    return v + f * direction;
}

} // namespace

LocalSingularityTerms buildPanelCenterLocalTerms(
    const Network& network, const Compressibility& flow, int panelI,
    int panelJ, const NetworkDoubletParameters& doubletParameters) {
    LocalSingularityTerms out;
    const int nPanelM = network.nm - 1;
    const int nPanelN = network.nn - 1;
    if (nPanelM <= 0 || nPanelN <= 0 || panelI < 0 || panelJ < 0 ||
        panelI >= nPanelM || panelJ >= nPanelN)
        return out;

    const SourceOuterSpline source =
        buildSourceOuterSpline(network, panelI, panelJ);
    if (!source.valid)
        return out;

    // P9 row of the source outer spline.
    for (int k = 0; k < 9; ++k) {
        const int id = source.neighborPanel[static_cast<std::size_t>(k)];
        if (id < 0)
            continue;
        out.sourceId.push_back(id);
        out.sourceValue.conservativeResize(out.sourceValue.size() + 1);
        out.sourceValue(out.sourceValue.size() - 1) =
            source.coefficients(4, k);
    }

    const DoubletPanelSpline doublet = buildDoubletPanelSpline(
        network, flow, panelI, panelJ, doubletParameters);
    if (!doublet.valid)
        return LocalSingularityTerms{};

    PanelGeometry panel = buildPanelGeometry(
        network.coordinates, network.nm, network.nn, panelI, panelJ);
    if (panel.degenerate)
        return LocalSingularityTerms{};
    applyCompressibleFrames(panel, flow);

    const auto doubletSubpanels = buildDoubletSubpanelSplines(panel, flow);
    const DoubletSubpanelSpline& inner = doubletSubpanels[4];
    if (!inner.valid)
        return LocalSingularityTerms{};
    const SubpanelGeometry& sp = panel.subpanels[4];
    const Eigen::Matrix3d toGlobal =
        sp.scaledValid ? sp.scaledToGlobal : sp.toGlobal;

    // Gradient of mu at P9 with respect to the nine canonical values.
    Eigen::Matrix<double, 3, 9> canonicalGradient;
    for (int j = 0; j < 9; ++j) {
        canonicalGradient.col(j) =
            toGlobal * Eigen::Vector3d(inner.coefficients(1, j),
                                       inner.coefficients(2, j), 0.0);
    }

    // Collect all local doublet parameter ids referenced by any row.
    auto findOrAdd = [&](int id) -> int {
        for (std::size_t k = 0; k < out.doubletId.size(); ++k) {
            if (out.doubletId[k] == id)
                return static_cast<int>(k);
        }
        out.doubletId.push_back(id);
        return static_cast<int>(out.doubletId.size() - 1);
    };
    for (int row = 0; row < 9; ++row) {
        for (int id : doublet.rows[static_cast<std::size_t>(row)].index)
            findOrAdd(id);
    }

    out.doubletValue = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(
        4, static_cast<Eigen::Index>(out.doubletId.size()));
    for (int row = 0; row < 9; ++row) {
        const DoubletRow& sourceRow =
            doublet.rows[static_cast<std::size_t>(row)];
        for (std::size_t k = 0; k < sourceRow.index.size(); ++k) {
            const int col = findOrAdd(sourceRow.index[k]);
            const double c = sourceRow.coefficient[k];
            out.doubletValue.col(col).segment<3>(1) +=
                c * canonicalGradient.col(row);
            if (row == 8)
                out.doubletValue(0, col) += c;
        }
    }

    out.valid = true;
    return out;
}

BoundaryRow buildBoundaryRow(
    const GlobalInfluence& influence, const Eigen::Vector3d& normal,
    const Compressibility& flow,
    const BoundaryConditionCoefficients& coefficients,
    const LocalSingularityTerms& local) {
    BoundaryRow out;
    if (!influence.valid || !local.valid)
        return out;
    if (coefficients.nct != 1 && coefficients.nct != 2 &&
        coefficients.nct != 4)
        return out;

    const double nnorm = normal.norm();
    if (!(nnorm > 0.0))
        return out;
    const Eigen::Vector3d n = normal / nnorm;

    const Eigen::Index nSource = influence.source.cols();
    const Eigen::Index nDoublet = influence.doublet.cols();
    out.source = Eigen::VectorXd::Zero(nSource);
    out.doublet = Eigen::VectorXd::Zero(nDoublet);
    out.rhs = coefficients.bet;

    if (coefficients.nct == 1) {
        const double ca = coefficients.cu + coefficients.cl;
        const double cd = 0.5 * (coefficients.cu - coefficients.cl);
        const Eigen::Vector3d tn =
            compressibilityScale(flow.betams, flow.compd, n);

        for (Eigen::Index j = 0; j < nSource; ++j) {
            const double v = tn.x() * influence.source(1, j) +
                             tn.y() * influence.source(2, j) +
                             tn.z() * influence.source(3, j);
            out.source(j) = ca * v;
        }
        for (Eigen::Index j = 0; j < nDoublet; ++j) {
            const double v = tn.x() * influence.doublet(1, j) +
                             tn.y() * influence.doublet(2, j) +
                             tn.z() * influence.doublet(3, j);
            out.doublet(j) = ca * v;
        }

        for (std::size_t k = 0; k < local.sourceId.size(); ++k) {
            const int id = local.sourceId[k];
            if (id < 0 || id >= nSource)
                return BoundaryRow{};
            out.source(id) +=
                cd * local.sourceValue(static_cast<Eigen::Index>(k));
        }
    } else if (coefficients.nct == 4) {
        const double da = coefficients.du + coefficients.dl;
        const double dd = 0.5 * (coefficients.du - coefficients.dl);

        out.source = da * influence.source.row(0).transpose();
        out.doublet = da * influence.doublet.row(0).transpose();

        for (std::size_t k = 0; k < local.doubletId.size(); ++k) {
            const int id = local.doubletId[k];
            if (id < 0 || id >= nDoublet)
                return BoundaryRow{};
            out.doublet(id) +=
                dd * local.doubletValue(0, static_cast<Eigen::Index>(k));
        }
    } else {
        // General (nct=2) coefficients: aical's ta/td projection.
        const double ca = coefficients.cu + coefficients.cl;
        const double cd = 0.5 * (coefficients.cu - coefficients.cl);
        const double da = coefficients.du + coefficients.dl;
        const double dd = 0.5 * (coefficients.du - coefficients.dl);
        const Eigen::Vector3d tn =
            compressibilityScale(flow.betams, flow.compd, n);

        auto projectTangential =
            [&](const Eigen::Vector3d& v) -> Eigen::Vector3d {
            const double v2 = v.squaredNorm();
            const double vn = v.dot(n);
            if (v2 > vn * vn + 1.0e-30)
                return std::sqrt(1.0 / (1.0 - vn * vn / v2)) *
                       (v - vn * n);
            return Eigen::Vector3d::Zero();
        };

        const Eigen::Vector3d ta =
            projectTangential(coefficients.tu + coefficients.tl) +
            ca * tn;
        const Eigen::Vector3d td =
            projectTangential(0.5 * (coefficients.tu - coefficients.tl)) +
            cd * tn;

        for (Eigen::Index j = 0; j < nSource; ++j) {
            out.source(j) =
                da * influence.source(0, j) +
                ta.x() * influence.source(1, j) +
                ta.y() * influence.source(2, j) +
                ta.z() * influence.source(3, j);
        }
        for (Eigen::Index j = 0; j < nDoublet; ++j) {
            out.doublet(j) =
                da * influence.doublet(0, j) +
                ta.x() * influence.doublet(1, j) +
                ta.y() * influence.doublet(2, j) +
                ta.z() * influence.doublet(3, j);
        }

        const double factor = tn.dot(n);
        if (std::abs(factor) <= 1.0e-14)
            return BoundaryRow{};
        const double sourceFactor = td.dot(n) / factor;
        for (std::size_t k = 0; k < local.sourceId.size(); ++k) {
            const int id = local.sourceId[k];
            if (id < 0 || id >= nSource)
                return BoundaryRow{};
            out.source(id) +=
                sourceFactor * local.sourceValue(static_cast<Eigen::Index>(k));
        }

        for (std::size_t k = 0; k < local.doubletId.size(); ++k) {
            const int id = local.doubletId[k];
            if (id < 0 || id >= nDoublet)
                return BoundaryRow{};
            const Eigen::Index kk = static_cast<Eigen::Index>(k);
            out.doublet(id) +=
                dd * local.doubletValue(0, kk) +
                td.dot(local.doubletValue.col(kk).segment<3>(1));
        }
    }

    out.valid = true;
    return out;
}

BoundaryConditionCoefficients standardNormalMassFlux(
    double cu, double cl, const Eigen::Vector3d& normal,
    const Eigen::Vector3d& freestream) {
    BoundaryConditionCoefficients out;
    out.nct = 1;
    out.cu = cu;
    out.cl = cl;
    const double nnorm = normal.norm();
    out.bet = (nnorm > 0.0)
                  ? -normal.dot(freestream) / nnorm
                  : 0.0;
    return out;
}

BoundarySystemRows buildNetworkNormalFluxRows(
    const Network& network, const Compressibility& flow,
    const NetworkSingularityLayout& layout,
    const std::vector<ControlPoint>& controlPoints, double cu, double cl) {
    BoundarySystemRows out;

    std::vector<BoundaryRow> rows;
    std::vector<int> indices;
    rows.reserve(controlPoints.size());
    indices.reserve(controlPoints.size());

    for (std::size_t k = 0; k < controlPoints.size(); ++k) {
        const ControlPoint& cp = controlPoints[k];
        if (cp.kind != ControlPointKind::PanelCenter)
            continue;

        const NetworkInfluence raw = evaluateNetworkNearFieldInfluence(
            network, flow, cp.position, layout);
        if (!raw.valid)
            return BoundarySystemRows{};

        GlobalInfluence global;
        global.source = raw.source;
        global.doublet = raw.doublet;
        global.valid = true;

        const LocalSingularityTerms local = buildPanelCenterLocalTerms(
            network, flow, cp.panelI, cp.panelJ, layout.doublet);
        if (!local.valid)
            return BoundarySystemRows{};

        const BoundaryConditionCoefficients bc =
            standardNormalMassFlux(cu, cl, cp.normal, flow.compd);
        BoundaryRow row = buildBoundaryRow(
            global, cp.normal, flow, bc, local);
        if (!row.valid)
            return BoundarySystemRows{};

        rows.push_back(std::move(row));
        indices.push_back(static_cast<int>(k));
    }

    if (rows.empty())
        return out;

    out.nRows = static_cast<int>(rows.size());
    out.nSource = layout.source.count;
    out.nDoublet = layout.doublet.count;
    out.source = Eigen::MatrixXd::Zero(out.nRows, out.nSource);
    out.doublet = Eigen::MatrixXd::Zero(out.nRows, out.nDoublet);
    out.rhs = Eigen::VectorXd::Zero(out.nRows);
    out.controlPointIndex = std::move(indices);

    for (int r = 0; r < out.nRows; ++r) {
        const BoundaryRow& row = rows[static_cast<std::size_t>(r)];
        if (row.source.size() != out.nSource ||
            row.doublet.size() != out.nDoublet)
            return BoundarySystemRows{};
        out.source.row(r) = row.source.transpose();
        out.doublet.row(r) = row.doublet.transpose();
        out.rhs(r) = row.rhs;
    }

    out.valid = true;
    return out;
}

} // namespace a502
