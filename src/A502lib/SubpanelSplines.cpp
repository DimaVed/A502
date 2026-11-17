// SubpanelSplines.cpp
// -------------------
// Source distribution on the eight triangular subpanels.

#include "SubpanelSplines.h"

#include "Compressibility.h"
#include "PanelParameters.h"
#include "TriangleInterpolation.h"

#include <array>
#include <cmath>

namespace a502 {

namespace {

constexpr double kEps = 1.0e-14;

// Value of the panel source distribution at canonical point P_{index+1} as a
// linear combination of [sigma1, sigma2, sigma3, sigma4, sigma9].
Eigen::Matrix<double, 1, 5> sourceBasis(int cpIndex) {
    Eigen::Matrix<double, 1, 5> b =
        Eigen::Matrix<double, 1, 5>::Zero();
    switch (cpIndex) {
    case 0: b(0) = 1.0; break;                 // P1
    case 1: b(1) = 1.0; break;                 // P2
    case 2: b(2) = 1.0; break;                 // P3
    case 3: b(3) = 1.0; break;                 // P4
    case 4: b(0) = 0.5; b(1) = 0.5; break;     // P5 = (P1+P2)/2
    case 5: b(1) = 0.5; b(2) = 0.5; break;     // P6 = (P2+P3)/2
    case 6: b(2) = 0.5; b(3) = 0.5; break;     // P7 = (P3+P4)/2
    case 7: b(3) = 0.5; b(0) = 0.5; break;     // P8 = (P4+P1)/2
    case 8: b(4) = 1.0; break;                 // P9
    default: break;
    }
    return b;
}

} // namespace

std::array<SourceSubpanelSpline, 8>
buildSourceSubpanelSplines(const PanelGeometry& panel) {
    std::array<SourceSubpanelSpline, 8> out;

    if (panel.degenerate)
        return out;

    for (int s = 0; s < 8; ++s) {
        const SubpanelGeometry& sp = panel.subpanels[s];
        SourceSubpanelSpline& spline = out[s];
        if (sp.degenerate)
            continue;

        const Eigen::Matrix3d& frame =
            sp.scaledValid ? sp.scaledToLocal : sp.toLocal;
        // subpqr uses the panel center as the origin of the inner
        // subpanels (ic > 4) and the first triangle corner otherwise.
        const Eigen::Vector3d origin =
            (s >= 4) ? panel.cp[8] : sp.origin;

        Eigen::Matrix3d M;
        Eigen::Matrix<double, 3, 5> V;
        for (int k = 0; k < 3; ++k) {
            const Eigen::Vector3d local =
                frame * (sp.vertex[k] - origin);
            M(k, 0) = 1.0;
            M(k, 1) = local.x();
            M(k, 2) = local.y();
            V.row(k) = sourceBasis(sp.cpIndex[k]);
        }

        const double det = M.determinant();
        if (std::abs(det) <= kEps)
            continue;

        spline.coefficients = M.inverse() * V;
        spline.origin = origin;
        spline.valid = true;
    }

    return out;
}

Eigen::Vector3d sourceTaylorCoefficients(
    const SourceSubpanelSpline& spline, const PanelSourceValues& values) {
    return spline.coefficients * values;
}

double evaluateSource(const SourceSubpanelSpline& spline,
                      const Eigen::Vector2d& local,
                      const PanelSourceValues& values) {
    if (!spline.valid)
        return 0.0;
    const Eigen::Vector3d c = sourceTaylorCoefficients(spline, values);
    return c(0) + c(1) * local.x() + c(2) * local.y();
}

Eigen::Matrix<double, 9, 4> buildAlam(const PanelGeometry& panel,
                                      const Compressibility& flow) {
    Eigen::Matrix<double, 9, 4> alam =
        Eigen::Matrix<double, 9, 4>::Zero();
    if (panel.degenerate)
        return alam;

    const std::array<Eigen::Vector3d, 4> corners = {
        panel.cp[0], panel.cp[1], panel.cp[2], panel.cp[3]};

    const NearPlane np = buildNearPlane(corners);
    if (!np.valid)
        return alam;

    const Eigen::Vector3d w = np.toNear * (panel.cp[0] - panel.cp[8]);
    const double c1 = w.x() - 1.0;
    const double c2 = w.y() - 1.0;

    const double fc[4] = {1.0 + c1 + c2, 1.0 + c1 - c2,
                          1.0 - c1 - c2, 1.0 - c1 + c2};
    const double fs[4] = {1.0 + c1, 1.0 - c2, 1.0 - c1, 1.0 + c2};

    for (int ic = 1; ic <= 4; ++ic) {
        const int icp1 = (ic % 4) + 1;
        const int icp2 = (icp1 % 4) + 1;
        const int icp3 = (icp2 % 4) + 1;
        const int icp4 = ic + 4;
        const int icp6 = icp2 + 4;

        const double wt1 = 1.0 / (fs[ic - 1] + fs[icp3 - 1]);
        const double wt2 = 1.0 / (fs[ic - 1] + fs[icp1 - 1]);

        alam(8, ic - 1) += wt1 * fc[ic - 1];
        alam(8, icp1 - 1) += wt2 * fc[icp1 - 1];

        alam(icp4 - 1, ic - 1) += wt1 * (fs[ic - 1] - 0.75 * fc[ic - 1]);
        alam(icp4 - 1, icp1 - 1) +=
            wt2 * (fs[ic - 1] - 0.75 * fc[icp1 - 1]);

        alam(icp6 - 1, ic - 1) -= 0.25 * wt1 * fc[ic - 1];
        alam(icp6 - 1, icp1 - 1) -= 0.25 * wt2 * fc[icp1 - 1];

        alam(ic - 1, ic - 1) += 0.25 * wt1;
        alam(ic - 1, icp1 - 1) -= 0.25 * wt2;
        alam(icp1 - 1, ic - 1) -= 0.25 * wt1;
        alam(icp1 - 1, icp1 - 1) += 0.25 * wt2;
    }

    std::array<Eigen::Vector3d, 4> scaledCorners;
    for (int i = 0; i < 4; ++i)
        scaledCorners[i] = flow.ggcp * panel.cp[i];

    for (int ic = 1; ic <= 4; ++ic) {
        const int icp4 = ic + 4;
        const int icp7 = ((ic + 2) % 4) + 5;
        const Eigen::Vector3d z =
            0.5 * (panel.cp[icp4 - 1] + panel.cp[icp7 - 1]);

        const PanelParameters par = closestPanelParameters(
            scaledCorners, flow.ggcp * z);
        if (!par.valid)
            continue;

        const BiquadraticBasis b = bqbfun(par.s, par.t);
        for (int i = 0; i < 9; ++i)
            alam(i, ic - 1) = 2.0 * b.phi(i);
        alam(icp4 - 1, ic - 1) -= 0.5;
        alam(icp7 - 1, ic - 1) -= 0.5;
    }

    return alam;
}

std::array<DoubletSubpanelSpline, 8>
buildDoubletSubpanelSplines(const PanelGeometry& panel,
                            const Compressibility& flow) {
    std::array<DoubletSubpanelSpline, 8> out;
    if (panel.degenerate)
        return out;

    const Eigen::Matrix<double, 9, 4> alam = buildAlam(panel, flow);

    // Panel corner coordinates in the inner subpanel frame, origin P9.
    const SubpanelGeometry& inner = panel.subpanels[4];
    const Eigen::Matrix3d innerFrame =
        inner.scaledValid ? inner.scaledToLocal : inner.toLocal;
    Eigen::Matrix<double, 3, 4> p;
    for (int i = 0; i < 4; ++i)
        p.col(i) = innerFrame * (panel.cp[i] - panel.cp[8]);

    for (int s = 0; s < 8; ++s) {
        const int ic = s + 1; // 1..8
        const int icp1 = (ic % 4) + 1;
        const int icp2 = (icp1 % 4) + 1;
        const int icp3 = (icp2 % 4) + 1;
        const int icp0 = (icp3 % 4) + 1;

        int ic0 = 0, ic1 = 0, ic2 = 0;
        int ice1 = 0, ice2 = 0;
        int icf1 = 0, icf2 = 0, icf3 = 0, icf4 = 0, icf5 = 0, icf6 = 0;

        Eigen::Matrix<double, 3, 3> e;
        Eigen::Matrix<double, 6, 9> f =
            Eigen::Matrix<double, 6, 9>::Zero();

        const bool innerTriangle = (ic > 4);
        if (!innerTriangle) {
            // psddq6 skips collapsed outer subpanels.
            if (panel.collapsedEdge == icp0 ||
                panel.collapsedEdge == icp3)
                continue;

            ic0 = icp0;
            ic1 = icp0 + 4;
            ic2 = icp3 + 4;
            ice1 = icp1;
            ice2 = icp3;
            icf1 = icp3 + 4;
            icf2 = icp0;
            icf3 = icp3;
            icf4 = icp0 + 4;
            icf5 = icp0;
            icf6 = icp1;

            e.row(0) << 1.0, p(0, ic0 - 1), p(1, ic0 - 1);
            e.row(1) << 1.0,
                0.5 * (p(0, ic0 - 1) + p(0, ice1 - 1)),
                0.5 * (p(1, ic0 - 1) + p(1, ice1 - 1));
            e.row(2) << 1.0,
                0.5 * (p(0, ic0 - 1) + p(0, ice2 - 1)),
                0.5 * (p(1, ic0 - 1) + p(1, ice2 - 1));
        } else {
            ic0 = 9;
            ic1 = icp3 + 4;
            ic2 = icp0 + 4;
            ice1 = icp3;
            ice2 = icp1;
            icf1 = 9;
            icf2 = icp0 + 4;
            icf3 = icp2 + 4;
            icf4 = 9;
            icf5 = icp3 + 4;
            icf6 = icp1 + 4;

            const double cx =
                0.25 * (p(0, 0) + p(0, 1) + p(0, 2) + p(0, 3));
            const double cy =
                0.25 * (p(1, 0) + p(1, 1) + p(1, 2) + p(1, 3));
            e.row(0) << 1.0, cx, cy;
            e.row(1) << 1.0,
                0.5 * (p(0, icp0 - 1) + p(0, ice1 - 1)),
                0.5 * (p(1, icp0 - 1) + p(1, ice1 - 1));
            e.row(2) << 1.0,
                0.5 * (p(0, icp0 - 1) + p(0, ice2 - 1)),
                0.5 * (p(1, icp0 - 1) + p(1, ice2 - 1));
        }

        f(0, ic0 - 1) = 1.0;
        f(1, ic1 - 1) = 1.0;
        f(2, ic2 - 1) = 1.0;
        f.row(3) = alam.col(icp0 - 1).transpose();
        f(4, icf1 - 1) = 1.0;
        f(4, icf2 - 1) = 0.25;
        f(4, icf3 - 1) = -0.25;
        f(5, icf4 - 1) = 1.0;
        f(5, icf5 - 1) = 0.25;
        f(5, icf6 - 1) = -0.25;

        const SubpanelGeometry& sp = panel.subpanels[s];
        const Eigen::Matrix3d frame =
            sp.scaledValid ? sp.scaledToLocal : sp.toLocal;
        const Eigen::Vector3d origin =
            innerTriangle ? panel.cp[8] : panel.cp[ic0 - 1];

        Eigen::Matrix<double, 3, 3> pp;
        pp.col(0).setZero();
        pp.col(1) = frame * (panel.cp[ic1 - 1] - origin);
        pp.col(2) = frame * (panel.cp[ic2 - 1] - origin);

        const TriangleQuadraticMap quad = buildQuadraticTriangleMap(pp);
        if (!quad.valid)
            continue;

        // e is kept for symmetry with subpqr(), where rr = g*e is the
        // source spline.  The doublet spline uses only f.
        (void)e;
        out[s].coefficients = quad.c * f;
        out[s].origin = origin;
        out[s].valid = true;
    }

    return out;
}

Eigen::Matrix<double, 6, 1> doubletTaylorCoefficients(
    const DoubletSubpanelSpline& spline, const PanelDoubletValues& values) {
    return spline.coefficients * values;
}

double evaluateDoublet(const DoubletSubpanelSpline& spline,
                       const Eigen::Vector2d& local,
                       const PanelDoubletValues& values) {
    if (!spline.valid)
        return 0.0;
    const Eigen::Matrix<double, 6, 1> c =
        doubletTaylorCoefficients(spline, values);
    return c(0) + c(1) * local.x() + c(2) * local.y() +
           c(3) * local.x() * local.x() / 2.0 +
           c(4) * local.x() * local.y() +
           c(5) * local.y() * local.y() / 2.0;
}

} // namespace a502
