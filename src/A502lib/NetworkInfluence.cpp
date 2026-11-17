// NetworkInfluence.cpp
// --------------------
// Assemble one AIC row for a single network from the source/doublet outer
// splines and the panel near-field PIC.

#include "NetworkInfluence.h"

#include "AeroModel.h"
#include "NearField.h"
#include "PanelGeometry.h"

namespace a502 {

NetworkSingularityLayout buildNetworkSingularityLayout(
    const Network& network) {
    NetworkSingularityLayout out;
    const int nm = network.nm;
    const int nn = network.nn;
    if (nm < 2 || nn < 2)
        return out;

    const std::size_t nPoints =
        static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
    if (network.coordinates.size() < nPoints * 3)
        return out;

    out.source.nPanelM = nm - 1;
    out.source.nPanelN = nn - 1;
    out.source.count = out.source.nPanelM * out.source.nPanelN;
    out.source.location.resize(
        static_cast<std::size_t>(out.source.count));
    for (int j = 0; j < out.source.nPanelN; ++j) {
        for (int i = 0; i < out.source.nPanelM; ++i) {
            const PanelGeometry panel = buildPanelGeometry(
                network.coordinates, nm, nn, i, j);
            out.source.location[static_cast<std::size_t>(
                out.source.panelId(i, j))] = panel.cp[8];
        }
    }

    out.doublet = buildNetworkDoubletParameters(network);
    if (out.doublet.count <= 0)
        return NetworkSingularityLayout{};
    return out;
}

NetworkInfluence evaluateNetworkNearFieldInfluence(
    const Network& network, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const NetworkSingularityLayout& layout) {
    NetworkInfluence out;
    const int nm = network.nm;
    const int nn = network.nn;
    const int nPanelM = nm - 1;
    const int nPanelN = nn - 1;
    if (nm < 2 || nn < 2 ||
        layout.source.nPanelM != nPanelM ||
        layout.source.nPanelN != nPanelN ||
        layout.doublet.nPanelM != nPanelM ||
        layout.doublet.nPanelN != nPanelN)
        return out;

    const std::size_t nPoints =
        static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
    if (network.coordinates.size() < nPoints * 3)
        return out;

    out.source = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(
        4, layout.source.count);
    out.doublet = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(
        4, layout.doublet.count);

    for (int panelJ = 0; panelJ < nPanelN; ++panelJ) {
        for (int panelI = 0; panelI < nPanelM; ++panelI) {
            PanelGeometry panel = buildPanelGeometry(
                network.coordinates, nm, nn, panelI, panelJ);
            if (panel.degenerate)
                continue;
            applyCompressibleFrames(panel, flow);

            const PanelNearFieldInfluence pic =
                subsonicPanelNearFieldInfluence(panel, flow, fieldPoint, 6);
            if (!pic.valid)
                return NetworkInfluence{};

            const SourceOuterSpline sourceSpline =
                buildSourceOuterSpline(network, panelI, panelJ);
            if (!sourceSpline.valid)
                return NetworkInfluence{};

            const Eigen::Matrix<double, 4, 9> sourceBlock =
                pic.source * sourceSpline.coefficients;
            for (int k = 0; k < 9; ++k) {
                const int id =
                    sourceSpline.neighborPanel[static_cast<std::size_t>(k)];
                if (id >= 0 && id < layout.source.count)
                    out.source.col(id) += sourceBlock.col(k);
            }

            const DoubletPanelSpline doubletSpline = buildDoubletPanelSpline(
                network, flow, panelI, panelJ, layout.doublet);
            if (!doubletSpline.valid)
                return NetworkInfluence{};

            for (int row = 0; row < 9; ++row) {
                const DoubletRow& splineRow =
                    doubletSpline.rows[static_cast<std::size_t>(row)];
                if (splineRow.index.size() != splineRow.coefficient.size())
                    return NetworkInfluence{};
                for (std::size_t k = 0; k < splineRow.index.size(); ++k) {
                    const int id = splineRow.index[k];
                    if (id < 0 || id >= layout.doublet.count)
                        return NetworkInfluence{};
                    out.doublet.col(id) +=
                        splineRow.coefficient[k] * pic.doublet.col(row);
                }
            }
        }
    }

    out.valid = true;
    return out;
}

NetworkAic buildNetworkNearFieldAic(
    const Network& network, const Compressibility& flow,
    const NetworkSingularityLayout& layout,
    const std::vector<Eigen::Vector3d>& fieldPoints) {
    NetworkAic out;
    if (fieldPoints.empty())
        return out;

    out.nControlPoints = static_cast<int>(fieldPoints.size());
    out.nSource = layout.source.count;
    out.nDoublet = layout.doublet.count;
    out.source = Eigen::MatrixXd::Zero(4 * out.nControlPoints, out.nSource);
    out.doublet = Eigen::MatrixXd::Zero(4 * out.nControlPoints, out.nDoublet);

    for (int k = 0; k < out.nControlPoints; ++k) {
        const NetworkInfluence row = evaluateNetworkNearFieldInfluence(
            network, flow, fieldPoints[static_cast<std::size_t>(k)],
            layout);
        if (!row.valid)
            return NetworkAic{};

        out.source.block(4 * k, 0, 4, out.nSource) = row.source;
        out.doublet.block(4 * k, 0, 4, out.nDoublet) = row.doublet;
    }

    out.valid = true;
    return out;
}

} // namespace a502
