// GlobalInfluence.cpp
// -------------------
// Multi-network singularity offsets and raw influence rows.

#include "GlobalInfluence.h"

#include "AeroModel.h"

namespace a502 {

GlobalSingularityLayout buildGlobalSingularityLayout(
    const std::vector<Network>& networks) {
    GlobalSingularityLayout out;
    if (networks.empty())
        return out;

    out.network.reserve(networks.size());
    out.sourceOffset.reserve(networks.size());
    out.doubletOffset.reserve(networks.size());
    out.activeSourceOffset.reserve(networks.size());
    out.activeDoubletOffset.reserve(networks.size());

    for (const Network& net : networks) {
        NetworkSingularityLayout layout =
            buildNetworkSingularityLayout(net);
        if (layout.source.count <= 0 || layout.doublet.count <= 0)
            return GlobalSingularityLayout{};

        out.sourceOffset.push_back(out.totalSource);
        out.doubletOffset.push_back(out.totalDoublet);
        out.totalSource += layout.source.count;
        out.totalDoublet += layout.doublet.count;

        if (layout.hasSource) {
            out.activeSourceOffset.push_back(out.activeSource);
            out.activeSource += layout.source.count;
        } else {
            out.activeSourceOffset.push_back(-1);
        }
        if (layout.hasDoublet) {
            out.activeDoubletOffset.push_back(out.activeDoublet);
            out.activeDoublet += layout.doublet.count;
        } else {
            out.activeDoubletOffset.push_back(-1);
        }

        out.network.push_back(std::move(layout));
    }

    return out;
}

GlobalInfluence evaluateGlobalNearFieldInfluence(
    const std::vector<Network>& networks, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint,
    const GlobalSingularityLayout& layout) {
    GlobalInfluence out;
    if (networks.empty() || layout.network.size() != networks.size() ||
        layout.sourceOffset.size() != networks.size() ||
        layout.doubletOffset.size() != networks.size())
        return out;

    out.source = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(
        4, layout.totalSource);
    out.doublet = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(
        4, layout.totalDoublet);

    for (std::size_t k = 0; k < networks.size(); ++k) {
        const NetworkInfluence row = evaluateNetworkNearFieldInfluence(
            networks[k], flow, fieldPoint, layout.network[k]);
        if (!row.valid)
            return GlobalInfluence{};

        const int sourceOffset =
            layout.sourceOffset[k];
        const int doubletOffset =
            layout.doubletOffset[k];
        out.source.block(0, sourceOffset, 4, row.source.cols()) =
            row.source;
        out.doublet.block(0, doubletOffset, 4, row.doublet.cols()) =
            row.doublet;
    }

    out.valid = true;
    return out;
}

GlobalAic buildGlobalNearFieldAic(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout,
    const std::vector<Eigen::Vector3d>& fieldPoints) {
    GlobalAic out;
    if (fieldPoints.empty() || layout.totalSource <= 0 ||
        layout.totalDoublet <= 0)
        return out;

    out.nControlPoints = static_cast<int>(fieldPoints.size());
    out.nSource = layout.totalSource;
    out.nDoublet = layout.totalDoublet;
    out.source = Eigen::MatrixXd::Zero(4 * out.nControlPoints, out.nSource);
    out.doublet = Eigen::MatrixXd::Zero(4 * out.nControlPoints, out.nDoublet);

    for (int k = 0; k < out.nControlPoints; ++k) {
        const GlobalInfluence row = evaluateGlobalNearFieldInfluence(
            networks, flow, fieldPoints[static_cast<std::size_t>(k)],
            layout);
        if (!row.valid)
            return GlobalAic{};

        out.source.block(4 * k, 0, 4, out.nSource) = row.source;
        out.doublet.block(4 * k, 0, 4, out.nDoublet) = row.doublet;
    }

    out.valid = true;
    return out;
}

} // namespace a502
