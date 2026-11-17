// PanelGeometry.cpp
// -----------------
// Canonical nine-point panel geometry and the eight triangular subpanels,
// following gensub / geomc / surfit / norcal in panair.f90.

#include "PanelGeometry.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace a502 {

namespace {

constexpr double kEps = 1.0e-12;

void makeFrame(SubpanelGeometry& s) {
    s.origin = s.vertex[0];
    s.centroid = (s.vertex[0] + s.vertex[1] + s.vertex[2]) / 3.0;

    Eigen::Vector3d z = s.normal;
    if (z.norm() <= kEps)
        z = Eigen::Vector3d::UnitZ();
    z.normalize();

    Eigen::Vector3d x = s.vertex[1] - s.vertex[0];
    x -= z * x.dot(z);
    if (x.norm() <= kEps) {
        x = s.vertex[2] - s.vertex[0];
        x -= z * x.dot(z);
    }
    if (x.norm() <= kEps)
        x = z.unitOrthogonal();
    x.normalize();

    const Eigen::Vector3d y = z.cross(x);

    s.toLocal.setZero();
    s.toLocal.row(0) = x.transpose();
    s.toLocal.row(1) = y.transpose();
    s.toLocal.row(2) = z.transpose();
    s.toGlobal = s.toLocal.transpose(); // orthonormal frame
}

} // namespace

PanelGeometry buildPanelGeometry(const std::vector<double>& coordinates,
                                 int nm, int nn, int i, int j) {
    PanelGeometry p;
    p.i = i;
    p.j = j;

    const std::size_t nPoints =
        static_cast<std::size_t>(std::max(nm, 0)) *
        static_cast<std::size_t>(std::max(nn, 0));
    if (nm < 2 || nn < 2 || i < 0 || j < 0 || i + 1 >= nm || j + 1 >= nn ||
        coordinates.size() < nPoints * 3) {
        p.degenerate = true;
        return p;
    }

    auto pointAt = [&](int m, int n) -> Eigen::Vector3d {
        const std::size_t g =
            static_cast<std::size_t>(n) * static_cast<std::size_t>(nm) +
            static_cast<std::size_t>(m);
        const std::size_t k = g * 3;
        return Eigen::Vector3d(coordinates[k], coordinates[k + 1],
                               coordinates[k + 2]);
    };

    // ---- canonical points P1..P9 -----------------------------------------
    p.cp[0] = pointAt(i, j);         // P1
    p.cp[1] = pointAt(i, j + 1);     // P2
    p.cp[2] = pointAt(i + 1, j + 1); // P3
    p.cp[3] = pointAt(i + 1, j);     // P4
    p.cp[4] = 0.5 * (p.cp[0] + p.cp[1]); // P5
    p.cp[5] = 0.5 * (p.cp[1] + p.cp[2]); // P6
    p.cp[6] = 0.5 * (p.cp[2] + p.cp[3]); // P7
    p.cp[7] = 0.5 * (p.cp[3] + p.cp[0]); // P8
    p.cp[8] = 0.25 * (p.cp[0] + p.cp[1] + p.cp[2] + p.cp[3]); // P9
    p.center = p.cp[8];

    double scale = 0.0;
    for (int k = 0; k < 4; ++k)
        for (int l = k + 1; l < 4; ++l)
            scale = std::max(scale, (p.cp[k] - p.cp[l]).norm());
    if (scale <= kEps) {
        p.degenerate = true;
        return p;
    }

    // ---- collapsed edges (iscal) ----------------------------------------
    const double collapseTol = 1.0e-13 * scale;
    for (int e = 0; e < 4; ++e) {
        const int e1 = e;
        const int e2 = (e + 1) % 4;
        if ((p.cp[e1] - p.cp[e2]).norm() <= collapseTol) {
            p.collapsedEdge = e + 1;
            break;
        }
    }

    // ---- panel normal: surfit's v1 x v2 --------------------------------
    const Eigen::Vector3d v1 =
        0.25 * (p.cp[0] + p.cp[3] - p.cp[1] - p.cp[2]);
    const Eigen::Vector3d v2 =
        0.25 * (p.cp[0] + p.cp[1] - p.cp[2] - p.cp[3]);
    Eigen::Vector3d nvec = v1.cross(v2);
    if (nvec.norm() <= kEps * scale * scale)
        nvec = (p.cp[1] - p.cp[0]).cross(p.cp[3] - p.cp[0]);
    if (nvec.norm() <= kEps * scale * scale) {
        p.degenerate = true;
        return p;
    }
    p.normal = nvec.normalized();

    // ---- eight subpanels (gensub order) --------------------------------
    static const int tri[8][3] = {
        {0, 4, 7}, // 1: P1 P5 P8
        {1, 5, 4}, // 2: P2 P6 P5
        {2, 6, 5}, // 3: P3 P7 P6
        {3, 7, 6}, // 4: P4 P8 P7
        {4, 8, 7}, // 5: P5 P9 P8
        {5, 8, 4}, // 6: P6 P9 P5
        {6, 8, 5}, // 7: P7 P9 P6
        {7, 8, 6}  // 8: P8 P9 P7
    };

    for (int s = 0; s < 8; ++s) {
        SubpanelGeometry& sp = p.subpanels[s];
        sp.cpIndex = {tri[s][0], tri[s][1], tri[s][2]};
        sp.vertex[0] = p.cp[tri[s][0]];
        sp.vertex[1] = p.cp[tri[s][1]];
        sp.vertex[2] = p.cp[tri[s][2]];

        Eigen::Vector3d cross =
            (sp.vertex[1] - sp.vertex[0])
                .cross(sp.vertex[2] - sp.vertex[0]);
        const double mag = cross.norm();
        if (mag <= kEps * scale * scale) {
            sp.degenerate = true;
            sp.area = 0.0;
            sp.normal = p.normal;
        } else {
            sp.normal = cross / mag;
            if (sp.normal.dot(p.normal) < 0.0)
                sp.normal = -sp.normal;
            sp.area = 0.5 * mag;
        }
        p.area += sp.area;
        makeFrame(sp);
    }

    if (p.area <= kEps * scale * scale)
        p.degenerate = true;

    return p;
}

} // namespace a502
