// NetworkInfluence.cpp
// --------------------
// Assemble one AIC row for a single network from the source/doublet outer
// splines and the panel near-field PIC.

#include "NetworkInfluence.h"

#include "AeroModel.h"
#include "NearField.h"
#include "PanelGeometry.h"
#include "Pic.h"

namespace a502 {

NetworkTypeDefaults networkTypeDefaults(int kt, double mach) {
    // Values from inbc(), labels 610..810 in panair.f90.
    NetworkTypeDefaults out;
    switch (kt) {
    case 1: // indirect mass flux b.c.
        out = {1, 12, 2, 5, 3, 7, 2};
        break;
    case 2: // thin wake (thickness flag not parsed; default tkw = 0)
        out = {0, 12, 0, 0, 0, 4, 3};
        break;
    case 3: // superinclined panels
        out = {1, 12, 0, 3, 2, 7, 2};
        break;
    case 4: // inlets
        out = {1, 12, 0, 2, 1, 7, 2};
        break;
    case 5: // bases
        if (mach > 1.0)
            out = {1, 12, 2, 2, 2, 6, 2};
        else
            out = {1, 12, 2, 6, 9, 7, 2};
        break;
    case 6: // off-body sources
        out = {1, 0, 0, 5, 2, 0, 2};
        break;
    case 7: // source alone, mass flux b.c.
        out = {1, 0, 0, 2, 3, 0, 2};
        break;
    case 8: // direct mass flux b.c., composite networks
        out = {1, 12, 2, 2, 3, 7, 2};
        break;
    case 9: // inlet b.c. with total mass flux
        out = {1, 12, 2, 2, 7, 7, 2};
        break;
    case 11: // indirect velocity b.c., composite networks
        out = {1, 12, 2, 14, 3, 7, 2};
        break;
    case 12: // thin wings, velocity b.c., doublets alone
        out = {0, 12, 0, 0, 3, 13, 3};
        break;
    case 13: // superinclined (A502 requires Mach > 1)
        out = {1, 12, 0, 12, 2, 7, 2};
        break;
    case 14: // inlets, velocity b.c.
        out = {1, 12, 0, 11, 1, 7, 2};
        break;
    case 15: // direct velocity b.c.
        out = {1, 12, 2, 11, 3, 7, 2};
        break;
    case 17: // source alone, velocity b.c.
        out = {1, 0, 0, 11, 3, 0, 2};
        break;
    case 18: // type 18 wake (source option and wake matching not parsed)
        out = {0, 18, 2, 0, 0, 15, 2};
        break;
    case 19: // type 6 design wake
        out = {0, 6, 0, 0, 0, 18, 2};
        break;
    case 20: // type 20 wake
        out = {0, 20, 2, 0, 0, 6, 2};
        break;
    case 21: // doublets alone
        out = {0, 12, 4, 0, 3, 7, 4};
        break;
    case 22:
        // Not accepted by the bundled inbc() dispatch, but other A502
        // versions describe kt=22 as doublets alone; keep that mapping.
        out = {0, 12, 4, 0, 3, 7, 4};
        break;
    default:
        // Unknown/arbitrary network type: keep the kt=1 composite defaults.
        break;
    }
    return out;
}

NetworkSingularityTypes networkSingularityTypes(int kt) {
    const NetworkTypeDefaults defaults = networkTypeDefaults(kt);
    NetworkSingularityTypes out;
    out.sourceType = defaults.sourceType;
    out.doubletType = defaults.doubletType;
    return out;
}

NetworkSingularityLayout buildNetworkSingularityLayout(
    const Network& network) {
    NetworkSingularityLayout out;
    const int nm = network.nm;
    const int nn = network.nn;
    if (nm < 2 || nn < 2)
        return out;

    out.types = networkSingularityTypes(network.kt);
    out.hasSource = out.types.sourceType != 0;
    out.hasDoublet = out.types.doubletType != 0;

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

namespace {

NetworkInfluence evaluateNetworkInfluenceSingle(
    const Network& network, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const NetworkSingularityLayout& layout, bool influenceDispatch) {
    NetworkInfluence out;
    const int nm = network.nm;
    const int nn = network.nn;
    const int nPanelM = nm - 1;
    const int nPanelN = nn - 1;
    if (nm < 2 || nn < 2 || (!layout.hasSource && !layout.hasDoublet))
        return out;
    if (layout.hasSource &&
        (layout.source.nPanelM != nPanelM ||
         layout.source.nPanelN != nPanelN))
        return out;
    if (layout.hasDoublet &&
        (layout.doublet.nPanelM != nPanelM ||
         layout.doublet.nPanelN != nPanelN))
        return out;

    const std::size_t nPoints =
        static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
    if (network.coordinates.size() < nPoints * 3)
        return out;

    // The raw row keeps the full per-network column layout so that global
    // offsets stay stable; columns of an inactive singularity type are zero.
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

            PanelNearFieldInfluence pic;
            if (influenceDispatch) {
                const int its =
                    (layout.hasSource ? 1 : 0) |
                    (layout.hasDoublet ? 2 : 0);
                pic = panelPicInfluence(panel, flow, fieldPoint, its, 6);
            } else {
                pic = panelNearFieldInfluence(panel, flow, fieldPoint, 6);
            }
            if (!pic.valid)
                return NetworkInfluence{};

            if (layout.hasSource) {
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
            }

            if (!layout.hasDoublet)
                continue;

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

} // namespace

namespace {

NetworkInfluence evaluateNetworkInfluenceImpl(
    const Network& network, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const NetworkSingularityLayout& layout,
    const SymmetryOptions& symmetry, bool influenceDispatch) {
    NetworkInfluence out;
    if (symmetry.misym == 0 && symmetry.mjsym == 0)
        return evaluateNetworkInfluenceSingle(
            network, flow, fieldPoint, layout, influenceDispatch);

    const int nImagesI = (symmetry.misym != 0) ? 2 : 1;
    const int nImagesJ = (symmetry.mjsym != 0) ? 2 : 1;
    out.source = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(
        4, layout.source.count);
    out.doublet = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(
        4, layout.doublet.count);

    for (int isymm = 1; isymm <= nImagesI; ++isymm) {
        const double sgni = (isymm == 1) ? 1.0 : -1.0;
        for (int jsymm = 1; jsymm <= nImagesJ; ++jsymm) {
            const double sgnj = (jsymm == 1) ? 1.0 : -1.0;
            double sgnk = 1.0;
            if (isymm * symmetry.misym == -2)
                sgnk = -sgnk;
            if (jsymm * symmetry.mjsym == -2)
                sgnk = -sgnk;

            const Eigen::Vector3d imagePoint(
                fieldPoint.x(), sgni * fieldPoint.y(),
                sgnj * fieldPoint.z());
            const NetworkInfluence row =
                evaluateNetworkInfluenceSingle(
                    network, flow, imagePoint, layout, influenceDispatch);
            if (!row.valid)
                return NetworkInfluence{};

            // Same sign pattern as the isymm/jsymm loop in the PIC caller:
            // potential/vx use sgnk, vy uses sgnk*sgni, vz uses sgnk*sgnj.
            const Eigen::Matrix<double, 4, 1> signs(
                sgnk, sgnk, sgnk * sgni, sgnk * sgnj);
            for (Eigen::Index j = 0; j < row.source.cols(); ++j)
                out.source.col(j) += signs.cwiseProduct(row.source.col(j));
            for (Eigen::Index j = 0; j < row.doublet.cols(); ++j)
                out.doublet.col(j) += signs.cwiseProduct(row.doublet.col(j));
        }
    }

    out.valid = true;
    return out;
}

NetworkAic buildNetworkAicImpl(
    const Network& network, const Compressibility& flow,
    const NetworkSingularityLayout& layout,
    const std::vector<Eigen::Vector3d>& fieldPoints,
    const SymmetryOptions& symmetry, bool influenceDispatch) {
    NetworkAic out;
    if (fieldPoints.empty())
        return out;

    out.nControlPoints = static_cast<int>(fieldPoints.size());
    out.nSource = layout.source.count;
    out.nDoublet = layout.doublet.count;
    out.source = Eigen::MatrixXd::Zero(4 * out.nControlPoints, out.nSource);
    out.doublet = Eigen::MatrixXd::Zero(4 * out.nControlPoints, out.nDoublet);

    for (int k = 0; k < out.nControlPoints; ++k) {
        const NetworkInfluence row = evaluateNetworkInfluenceImpl(
            network, flow, fieldPoints[static_cast<std::size_t>(k)],
            layout, symmetry, influenceDispatch);
        if (!row.valid)
            return NetworkAic{};

        out.source.block(4 * k, 0, 4, out.nSource) = row.source;
        out.doublet.block(4 * k, 0, 4, out.nDoublet) = row.doublet;
    }

    out.valid = true;
    return out;
}

} // namespace

NetworkInfluence evaluateNetworkNearFieldInfluence(
    const Network& network, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const NetworkSingularityLayout& layout,
    const SymmetryOptions& symmetry) {
    return evaluateNetworkInfluenceImpl(network, flow, fieldPoint, layout,
                                        symmetry, false);
}

NetworkInfluence evaluateNetworkPicInfluence(
    const Network& network, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const NetworkSingularityLayout& layout,
    const SymmetryOptions& symmetry) {
    return evaluateNetworkInfluenceImpl(network, flow, fieldPoint, layout,
                                        symmetry, true);
}

NetworkAic buildNetworkNearFieldAic(
    const Network& network, const Compressibility& flow,
    const NetworkSingularityLayout& layout,
    const std::vector<Eigen::Vector3d>& fieldPoints,
    const SymmetryOptions& symmetry) {
    return buildNetworkAicImpl(network, flow, layout, fieldPoints, symmetry,
                               false);
}

NetworkAic buildNetworkPicAic(
    const Network& network, const Compressibility& flow,
    const NetworkSingularityLayout& layout,
    const std::vector<Eigen::Vector3d>& fieldPoints,
    const SymmetryOptions& symmetry) {
    return buildNetworkAicImpl(network, flow, layout, fieldPoints, symmetry,
                               true);
}

} // namespace a502
