// Pic.cpp
// -------
// Influence-level dispatch of the panel influence coefficients.

#include "Pic.h"

#include "FarField.h"
#include "Influence.h"
#include "PanelGeometry.h"
#include "PanelParameters.h"
#include "SubpanelSplines.h"

namespace a502 {

PanelNearFieldInfluence panelIntermediateInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int its, int nf) {
    PanelNearFieldInfluence out;
    if (panel.degenerate || nf != 6 || its < 1 || its > 3)
        return out;

    const SubpanelGeometry& inner = panel.subpanels[4];
    if (!inner.scaledValid)
        return out;

    const Eigen::Matrix3d& frame = inner.scaledToLocal;
    const Eigen::Vector3d origin = panel.cp[8];

    // qffcal() passes the four panel corners (pf) and the collapsed-edge
    // index to nftpic() as a single four-sided subpanel.
    std::array<Eigen::Vector3d, 4> vertices{};
    for (int k = 0; k < 4; ++k)
        vertices[static_cast<std::size_t>(k)] =
            frame *
            (panel.cp[static_cast<std::size_t>(k)] - origin);
    const Eigen::Vector3d localField = frame * (fieldPoint - origin);

    SubsonicSubpanelInfluence infl;
    if (flow.mach > 1.0 && inner.inclination < 0) {
        const double sgnx =
            (flow.compd.dot(inner.normal) >= 0.0) ? 1.0 : -1.0;
        infl = supersonicSuperinclinedSubpanelInfluence(
            vertices, 4, panel.collapsedEdge, its, localField,
            inner.jacobian, nf, sgnx);
    } else if (flow.mach > 1.0 && inner.inclination > 0) {
        infl = supersonicSubpanelInfluence(
            vertices, 4, panel.collapsedEdge, its, localField,
            inner.jacobian, nf);
    } else if (flow.mach <= 1.0) {
        infl = subsonicSubpanelInfluence(
            vertices, 4, panel.collapsedEdge, its, localField,
            inner.jacobian, nf);
    }
    if (!infl.valid)
        return out;

    // qffcal() zeroes the inactive block before calling nftpic(); the
    // kernel may still fill its scratch output.
    if (its == 1)
        infl.doublet.setZero();
    if (its == 2)
        infl.source.setZero();

    const auto sourceSplines = buildSourceSubpanelSplines(panel);
    if (!sourceSplines[4].valid)
        return out;
    const NearPlaneQuadratic qa =
        buildNearPlaneQuadratic(panel.cp, frame);
    if (!qa.valid)
        return out;

    // The caller of qffcal() applies asts to the source Taylor coefficients
    // and qa to the doublet Taylor coefficients.
    out.source = infl.source.block<4, 3>(0, 0) *
                 sourceSplines[4].coefficients;
    out.doublet = infl.doublet.block<4, 6>(0, 0) * qa.coefficients;

    for (int col = 0; col < 5; ++col) {
        const Eigen::Vector3d v =
            inner.scaledToGlobal *
            Eigen::Vector3d(out.source(1, col), out.source(2, col),
                            out.source(3, col));
        out.source(1, col) = v.x();
        out.source(2, col) = v.y();
        out.source(3, col) = v.z();
    }
    for (int col = 0; col < 9; ++col) {
        const Eigen::Vector3d v =
            inner.scaledToGlobal *
            Eigen::Vector3d(out.doublet(1, col), out.doublet(2, col),
                            out.doublet(3, col));
        out.doublet(1, col) = v.x();
        out.doublet(2, col) = v.y();
        out.doublet(3, col) = v.z();
    }

    out.valid = true;
    return out;
}

PanelNearFieldInfluence panelPicInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int its, int nf, int level) {
    PanelNearFieldInfluence out;
    if (its < 1 || its > 3 || nf != 6)
        return out;

    int iflun = level;
    if (iflun < 0)
        iflun = influenceLevel(panel, flow, fieldPoint);
    if (iflun <= 0) {
        out.valid = true;
        return out;
    }

    if (iflun <= 3) {
        const PanelFarFieldInfluence far = farFieldPanelInfluence(
            panel, flow, fieldPoint, iflun, nf, its);
        if (!far.valid)
            return out;
        out.source = far.source;
        out.doublet = far.doublet;
        out.valid = true;
        return out;
    }

    if (iflun == 4) {
        const PanelNearFieldInfluence intermediate =
            panelIntermediateInfluence(panel, flow, fieldPoint, its, nf);
        if (intermediate.valid)
            return intermediate;
        // Fall through to the near field if the single-panel frame is not
        // available.
    }

    PanelNearFieldInfluence near =
        panelNearFieldInfluence(panel, flow, fieldPoint, nf);
    if (near.valid) {
        if (its == 1)
            near.doublet.setZero();
        if (its == 2)
            near.source.setZero();
    }
    return near;
}

} // namespace a502
