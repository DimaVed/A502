// PanelParameters.cpp
// -------------------
// Bilinear panel parameterization and near-plane transform.

#include "PanelParameters.h"

#include <algorithm>
#include <cmath>

namespace a502 {

namespace {
constexpr double kEps = 1.0e-14;
}

Eigen::Vector3d evalBilinearPanel(
    const std::array<Eigen::Vector3d, 4>& c, double s, double t) {
    const Eigen::Vector3d qc = 0.25 * (c[0] + c[1] + c[2] + c[3]);
    const Eigen::Vector3d ps = 0.25 * (c[0] - c[1] - c[2] + c[3]);
    const Eigen::Vector3d pt = 0.25 * (c[0] + c[1] - c[2] - c[3]);
    const Eigen::Vector3d pst = 0.25 * (c[0] - c[1] + c[2] - c[3]);
    return qc + s * ps + t * pt + s * t * pst;
}

PanelParameters closestPanelParameters(
    const std::array<Eigen::Vector3d, 4>& c, const Eigen::Vector3d& p) {
    PanelParameters out;

    const Eigen::Vector3d qc = 0.25 * (c[0] + c[1] + c[2] + c[3]);
    const Eigen::Vector3d ps = 0.25 * (c[0] - c[1] - c[2] + c[3]);
    const Eigen::Vector3d pt = 0.25 * (c[0] + c[1] - c[2] - c[3]);
    const Eigen::Vector3d pst = 0.25 * (c[0] - c[1] + c[2] - c[3]);

    double scale = 1.0;
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j)
            scale = std::max(scale, (c[i] - c[j]).norm());
    const double tol = 1.0e-12 * scale;

    double s = 0.0;
    double t = 0.0;
    bool converged = false;

    for (int iter = 0; iter < 40; ++iter) {
        const Eigen::Vector3d Q = qc + s * ps + t * pt + s * t * pst;
        const Eigen::Vector3d Qs = ps + t * pst;
        const Eigen::Vector3d Qt = pt + s * pst;
        const Eigen::Vector3d r = Q - p;

        const double f1 = r.dot(Qs);
        const double f2 = r.dot(Qt);
        if (std::max(std::abs(f1), std::abs(f2)) <= tol) {
            converged = true;
            break;
        }

        const double h11 = Qs.dot(Qs);
        const double h12 = Qs.dot(Qt) + r.dot(pst);
        const double h22 = Qt.dot(Qt);
        const double det = h11 * h22 - h12 * h12;
        if (std::abs(det) <= kEps * scale * scale)
            break;

        const double ds = -(h22 * f1 - h12 * f2) / det;
        const double dt = -(-h12 * f1 + h11 * f2) / det;

        double sNew = s + ds;
        double tNew = t + dt;
        if (!std::isfinite(sNew) || !std::isfinite(tNew))
            break;

        sNew = std::max(-1.0, std::min(1.0, sNew));
        tNew = std::max(-1.0, std::min(1.0, tNew));

        if (std::abs(sNew - s) + std::abs(tNew - t) <= 1.0e-14) {
            s = sNew;
            t = tNew;
            converged = true;
            break;
        }
        s = sNew;
        t = tNew;
    }

    out.s = s;
    out.t = t;
    out.valid = converged;
    return out;
}

NearPlane buildNearPlane(const std::array<Eigen::Vector3d, 4>& c) {
    NearPlane out;

    const Eigen::Vector3d v1 = 0.25 * (c[0] + c[3] - c[1] - c[2]);
    const Eigen::Vector3d v2 = 0.25 * (c[0] + c[1] - c[2] - c[3]);
    const Eigen::Vector3d v3 = v1.cross(v2);
    const double mag = v3.norm();
    if (mag <= kEps)
        return out;

    const double artm = std::sqrt(mag); // surfit uses sqrt(|v1 x v2|)
    if (artm <= kEps)
        return out;

    Eigen::Matrix3d art;
    art.col(0) = v1;
    art.col(1) = v2;
    art.col(2) = v3 / artm;

    const double det = art.determinant();
    if (std::abs(det) <= kEps)
        return out;

    out.toNear = art.inverse();
    out.valid = true;
    return out;
}

} // namespace a502
