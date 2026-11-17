// BoundaryConditions.cpp
// ----------------------
// Port of the nct=1/nct=4 core of aical().

#include "BoundaryConditions.h"

#include "AeroModel.h"
#include "OuterSplines.h"

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

    const DoubletRow& centerRow = doublet.rows[8];
    out.doubletId = centerRow.index;
    out.doubletValue.resize(
        static_cast<Eigen::Index>(centerRow.coefficient.size()));
    for (std::size_t k = 0; k < centerRow.coefficient.size(); ++k)
        out.doubletValue(static_cast<Eigen::Index>(k)) =
            centerRow.coefficient[k];

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
    if (coefficients.nct != 1 && coefficients.nct != 4)
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
    } else {
        const double da = coefficients.du + coefficients.dl;
        const double dd = 0.5 * (coefficients.du - coefficients.dl);

        out.source = da * influence.source.row(0).transpose();
        out.doublet = da * influence.doublet.row(0).transpose();

        for (std::size_t k = 0; k < local.doubletId.size(); ++k) {
            const int id = local.doubletId[k];
            if (id < 0 || id >= nDoublet)
                return BoundaryRow{};
            out.doublet(id) +=
                dd * local.doubletValue(static_cast<Eigen::Index>(k));
        }
    }

    out.valid = true;
    return out;
}

} // namespace a502
