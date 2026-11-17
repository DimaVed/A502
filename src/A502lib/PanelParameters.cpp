// PanelParameters.cpp
// -------------------
// Bilinear panel parameterization and near-plane transform.

#include "PanelParameters.h"

#include <algorithm>
#include <cmath>
#include <vector>

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

namespace {

// lchvar(): apply x = x0 + dx*tau to the polynomial coefficients
// v[k] = coefficient of tau^(k-1), in place.  base is a 0-based offset into
// the flat buffer; inc matches the Fortran stride parameter.
void lchvar(std::vector<double>& v, std::size_t base, double x0, double dx,
            int idim, int inc) {
    int last = 1 + idim;
    for (int idbk = 1; idbk <= idim; ++idbk) {
        const int id = 1 + idim - idbk;
        const int lstinc = last + inc;
        for (int j = 1; j <= id; ++j) {
            const std::size_t dst =
                base + static_cast<std::size_t>(lstinc + 1 - j) - 1;
            const std::size_t srcDx =
                base + static_cast<std::size_t>(last + 1 - j) - 1;
            const std::size_t srcX0 =
                base + static_cast<std::size_t>(last - j) - 1;
            v[dst] = dx * v[srcDx] + x0 * v[srcX0];
        }
        last = lstinc;
    }
}

} // namespace

Eigen::MatrixXd panelAreaMoments(
    const std::array<Eigen::Vector3d, 4>& vertices, int nside, int ics,
    int degree) {
    Eigen::MatrixXd c = Eigen::MatrixXd::Zero(degree + 1, degree + 1);
    if (nside < 3 || nside > 4 || degree < 0 || degree > 8)
        return c;

    const int np1 = degree + 1;
    const int np2 = degree + 2;
    std::vector<double> e(
        static_cast<std::size_t>(np2) * static_cast<std::size_t>(np1),
        0.0);

    for (int ic = 1; ic <= nside; ++ic) {
        if (ic == ics)
            continue;
        int icp1 = ic % nside + 1;
        if (icp1 == ics)
            icp1 = ics % nside + 1;

        const Eigen::Vector3d& p1 =
            vertices[static_cast<std::size_t>(ic - 1)];
        const Eigen::Vector3d& p2 =
            vertices[static_cast<std::size_t>(icp1 - 1)];
        const double delx = p2.x() - p1.x();
        const double dely = p2.y() - p1.y();

        for (int i = 1; i <= np2; ++i)
            e[static_cast<std::size_t>(i - 1)] = dely / i;

        lchvar(e, 0, p1.y(), dely, np1, np1);
        for (int j = 1; j <= np1; ++j) {
            const std::size_t base =
                static_cast<std::size_t>(j - 1) *
                static_cast<std::size_t>(np2);
            lchvar(e, base, p1.x(), delx, np2 - j, 0);
        }

        for (int i = 1; i <= np1; ++i) {
            const double aiinv = 1.0 / i;
            const int jmax = np2 - i;
            for (int j = 1; j <= jmax; ++j) {
                const std::size_t idx =
                    static_cast<std::size_t>(j - 1) *
                        static_cast<std::size_t>(np2) +
                    static_cast<std::size_t>(i + 1) - 1;
                c(i - 1, j - 1) += aiinv * e[idx];
            }
        }
    }
    return c;
}

NearPlaneQuadratic buildNearPlaneQuadratic(
    const std::array<Eigen::Vector3d, 9>& cp,
    const Eigen::Matrix3d& innerToLocal) {
    NearPlaneQuadratic out;

    // qcof(): pc are the bilinear shape vectors of the four corners.
    Eigen::Matrix<double, 3, 3> pc;
    for (int l = 0; l < 3; ++l) {
        pc(l, 0) =
            0.25 * (cp[0][l] - cp[1][l] - cp[2][l] + cp[3][l]);
        pc(l, 1) =
            0.25 * (cp[0][l] + cp[1][l] - cp[2][l] - cp[3][l]);
        pc(l, 2) =
            0.25 * (cp[0][l] - cp[1][l] + cp[2][l] - cp[3][l]);
    }
    const Eigen::Matrix3d z = innerToLocal * pc;

    const double kses = z(0, 0);
    const double etas = z(1, 0);
    const double kset = z(0, 1);
    const double etat = z(1, 1);
    const double ksest = z(0, 2);
    const double etast = z(1, 2);

    const double det = kses * etat - kset * etas;
    if (std::abs(det) <= kEps)
        return out;
    const double al0 = 1.0 / det;
    const double al0s = al0 * al0;
    const double deta =
        al0 * (2.0 * etas * etat * ksest -
               (etat * kses + etas * kset) * etast);
    const double dkse =
        -al0 * (2.0 * kses * kset * etast -
                (etat * kses + etas * kset) * ksest);

    Eigen::Matrix<double, 6, 6> c =
        Eigen::Matrix<double, 6, 6>::Zero();
    c(0, 0) = 1.0;
    c(1, 1) = al0 * etat;
    c(1, 2) = -al0 * etas;
    c(2, 1) = -al0 * kset;
    c(2, 2) = al0 * kses;
    c(3, 1) = al0s * etat * (etast + deta);
    c(3, 2) = al0s * etas * (etast - deta);
    c(3, 3) = al0s * etat * etat;
    c(3, 4) = -2.0 * al0s * etas * etat;
    c(3, 5) = al0s * etas * etas;
    c(4, 1) = -al0s * (kset * etast + etat * dkse);
    c(4, 2) = -al0s * (etas * ksest - kses * deta);
    c(4, 3) = -al0s * kset * etat;
    c(4, 4) = al0s * (kset * etas + etat * kses);
    c(4, 5) = -al0s * kses * etas;
    c(5, 1) = al0s * kset * (ksest + dkse);
    c(5, 2) = al0s * kses * (ksest - dkse);
    c(5, 3) = al0s * kset * kset;
    c(5, 4) = -2.0 * al0s * kses * kset;
    c(5, 5) = al0s * kses * kses;

    for (int k = 0; k < 6; ++k) {
        out.coefficients(k, 0) = 0.25 * c(k, 4);
        out.coefficients(k, 1) = -0.25 * c(k, 4);
        out.coefficients(k, 2) = 0.25 * c(k, 4);
        out.coefficients(k, 3) = -0.25 * c(k, 4);
        out.coefficients(k, 4) = 0.5 * c(k, 2) + c(k, 5);
        out.coefficients(k, 5) = -0.5 * c(k, 1) + c(k, 3);
        out.coefficients(k, 6) = -0.5 * c(k, 2) + c(k, 5);
        out.coefficients(k, 7) = 0.5 * c(k, 1) + c(k, 3);
        out.coefficients(k, 8) =
            c(k, 0) - 2.0 * c(k, 3) - 2.0 * c(k, 5);
    }

    out.valid = true;
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
    out.normal = art.col(2);
    out.valid = true;
    return out;
}

} // namespace a502
