// ControlPoints.cpp
// -----------------
// Control point withdrawal from panel centers, network edge midpoints and
// network corners (contrl, gcpcal, zcadj, xxadj in panair.f90).

#include "ControlPoints.h"

#include "AeroModel.h"
#include "PanelGeometry.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace a502 {

namespace {

constexpr double kCenterDelta = 0.01;  // contrl: data delta
constexpr double kEdgeDelta   = 0.125; // zcadj: deltam
constexpr double kCornerDelta = 0.25;  // zcadj: deltac

// xxadj: withdraw p0 a fraction al of the distance along the angle bisector
// of the triangle (p0,p1,p2).
Eigen::Vector3d xxadj(const Eigen::Vector3d& p0,
                      const Eigen::Vector3d& p1,
                      const Eigen::Vector3d& p2,
                      double al) {
    const Eigen::Vector3d v1 = p1 - p0;
    const Eigen::Vector3d v2 = p2 - p0;
    const double m1 = v1.norm();
    const double m2 = v2.norm();
    if (m1 + m2 <= 1.0e-15)
        return p0;
    return p0 + al * (m2 * v1 + m1 * v2) / (m1 + m2);
}

int oppositeEdge(int edge) {
    // zcadj: isp3 = mod(is+2,4)+1 for is=1..4.
    static const int opp[5] = {0, 4, 1, 2, 3};
    return opp[edge];
}

int cornerOfEdge(int edge) {
    // edge 1 = P1-P2, 2 = P2-P3, 3 = P3-P4, 4 = P4-P1.
    return edge;
}

} // namespace

std::vector<ControlPoint> buildNetworkControlPoints(const Network& network) {
    std::vector<ControlPoint> out;
    const int nm = network.nm;
    const int nn = network.nn;
    if (nm < 2 || nn < 2)
        return out;

    const std::size_t nPoints =
        static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
    if (network.coordinates.size() < nPoints * 3)
        return out;

    const int nPanelsM = nm - 1;
    const int nPanelsN = nn - 1;
    std::vector<PanelGeometry> panels(
        static_cast<std::size_t>(nPanelsM) * nPanelsN);

    double scale = 1.0;
    for (std::size_t i = 0; i + 2 < network.coordinates.size(); i += 3) {
        scale = std::max(scale,
                         std::abs(network.coordinates[i]) +
                             std::abs(network.coordinates[i + 1]) +
                             std::abs(network.coordinates[i + 2]));
    }
    const double dedupTol = 1.0e-10 * scale;

    auto addUnique = [&](const ControlPoint& cp) {
        for (const auto& p : out) {
            if ((p.position - cp.position).norm() <= dedupTol)
                return; // coincident control points (collapsed edge/corner)
        }
        out.push_back(cp);
    };

    auto panelAt = [&](int i, int j) -> PanelGeometry& {
        return panels[static_cast<std::size_t>(i) +
                      static_cast<std::size_t>(nPanelsM) * j];
    };
    const auto panelAtC = [&](int i, int j) -> const PanelGeometry& {
        return panels[static_cast<std::size_t>(i) +
                      static_cast<std::size_t>(nPanelsM) * j];
    };

    for (int j = 0; j < nPanelsN; ++j)
        for (int i = 0; i < nPanelsM; ++i)
            panelAt(i, j) = buildPanelGeometry(network.coordinates, nm, nn, i, j);

    constexpr int kFirstPanel = 0; // panelIndex = i + nPanelsM*j
    auto fillCommon = [&](ControlPoint& cp, int i, int j) {
        cp.network = 0; // AeroModel overwrites this when appending
        cp.panelI = i;
        cp.panelJ = j;
        cp.panelIndex = i + nPanelsM * j + kFirstPanel;
    };

    // ---- panel centers ---------------------------------------------------
    for (int j = 0; j < nPanelsN; ++j) {
        for (int i = 0; i < nPanelsM; ++i) {
            const PanelGeometry& g = panelAtC(i, j);
            if (g.degenerate)
                continue;
            ControlPoint cp;
            cp.kind = ControlPointKind::PanelCenter;
            cp.position = xxadj(g.cp[8], g.cp[4], g.cp[7], kCenterDelta);
            cp.normal = g.normal;
            fillCommon(cp, i, j);
            addUnique(cp);
        }
    }

    // ---- edge midpoints on network edges --------------------------------
    auto addEdge = [&](int i, int j, int edge) {
        const PanelGeometry& g = panelAtC(i, j);
        if (g.degenerate)
            return;

        const int midIndex = 3 + edge; // edge 1 -> P5 ... edge 4 -> P8
        const int opp = oppositeEdge(edge);
        const Eigen::Vector3d mid = g.cp[midIndex];
        const Eigen::Vector3d center = g.cp[8];
        const Eigen::Vector3d oppMid = g.cp[3 + opp];
        const Eigen::Vector3d corner = g.cp[cornerOfEdge(edge) - 1];

        const double d4 = (mid - center).norm();
        const double d7 = (oppMid - center).norm();
        const double di = (mid - corner).norm();
        const double dm = std::min({d4, d7, di});
        double delta = (d4 > 1.0e-15) ? kEdgeDelta * dm / d4 : kEdgeDelta;
        if (g.collapsedEdge == edge)
            delta = kCornerDelta;

        ControlPoint cp;
        cp.kind = ControlPointKind::EdgeMidpoint;
        cp.position = xxadj(mid, center, oppMid, delta);
        // The withdrawn edge midpoint lies in the adjacent inner subpanel,
        // whose normal equals the panel normal (en(5) in geomc).
        cp.normal = g.normal;
        cp.edge = edge;
        fillCommon(cp, i, j);
        addUnique(cp);
    };

    for (int j = 0; j < nPanelsN; ++j) {
        for (int i = 0; i < nPanelsM; ++i) {
            if (i == 0)
                addEdge(i, j, 1); // P1-P2, first row
            if (j == nPanelsN - 1)
                addEdge(i, j, 2); // P2-P3, last column
            if (i == nPanelsM - 1)
                addEdge(i, j, 3); // P3-P4, last row
            if (j == 0)
                addEdge(i, j, 4); // P4-P1, first column
        }
    }

    // ---- network corners -------------------------------------------------
    auto addCorner = [&](int i, int j, int corner) {
        const PanelGeometry& g = panelAtC(i, j);
        if (g.degenerate)
            return;

        // Adjacent edge midpoints: corner 1 -> edges 1,4; 2 -> 2,1;
        // 3 -> 3,2; 4 -> 4,3.
        static const int edgeA[4] = {1, 2, 3, 4};
        static const int edgeB[4] = {4, 1, 2, 3};
        const Eigen::Vector3d p0 = g.cp[corner - 1];
        const Eigen::Vector3d p1 = g.cp[3 + edgeA[corner - 1]];
        const Eigen::Vector3d p2 = g.cp[3 + edgeB[corner - 1]];

        ControlPoint cp;
        cp.kind = ControlPointKind::Corner;
        cp.position = xxadj(p0, p1, p2, kCornerDelta);
        cp.normal = g.subpanels[corner - 1].degenerate
                        ? g.normal
                        : g.subpanels[corner - 1].normal;
        cp.corner = corner;
        fillCommon(cp, i, j);
        addUnique(cp);
    };

    addCorner(0, 0, 1);
    addCorner(0, nPanelsN - 1, 2);
    addCorner(nPanelsM - 1, nPanelsN - 1, 3);
    addCorner(nPanelsM - 1, 0, 4);

    return out;
}

} // namespace a502
