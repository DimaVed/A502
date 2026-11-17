// NearField.cpp
// -------------
// Subsonic near-field source/doublet subpanel influence coefficients.  This
// is a direct port of subsbi() from panair.f90.

#include "NearField.h"

#include "Compressibility.h"
#include "PanelGeometry.h"
#include "SubpanelSplines.h"

#include <cmath>

namespace a502 {

namespace {
constexpr double kPi = 3.14159265358979323846;
}

SubsonicSubpanelInfluence subsonicSubpanelInfluence(
    const std::array<Eigen::Vector3d, 4>& vertices, int ns, int ics, int its,
    const Eigen::Vector3d& fieldPoint, double aj, int nf) {
    SubsonicSubpanelInfluence out;
    if (ns < 3 || ns > 4 || (nf != 6 && nf != 10) || its < 1 || its > 3)
        return out;

    const double h = fieldPoint.z();
    const double signh = (h > 0.0) ? 1.0 : (h < 0.0 ? -1.0 : 0.0);
    const double hh = h * h;
    const double hm = std::abs(h);

    std::array<double, 13> b{};
    for (int edge = 1; edge <= ns; ++edge) {
        if (edge == ics)
            continue;
        int next = (edge % ns) + 1;
        if (next == ics)
            next = (next % ns) + 1;

        const Eigen::Vector3d& p1 =
            vertices[static_cast<std::size_t>(edge - 1)];
        const Eigen::Vector3d& p2 =
            vertices[static_cast<std::size_t>(next - 1)];
        const double aks1 = p1.x() - fieldPoint.x();
        const double aet1 = p1.y() - fieldPoint.y();
        const double aks2 = p2.x() - fieldPoint.x();
        const double aet2 = p2.y() - fieldPoint.y();

        const double dks = aks2 - aks1;
        const double det = aet2 - aet1;
        const double drm = std::sqrt(dks * dks + det * det);
        if (!(drm > 0.0))
            continue;
        const double drmi = 1.0 / drm;
        const double ank = drmi * det;
        const double ane = -drmi * dks;
        const double a = aks1 * ank + aet1 * ane;
        const double aa = a * a;
        const double gg = aa + hh;
        const double el1 = aet1 * ank - aks1 * ane;
        const double el2 = aet2 * ank - aks2 * ane;
        const double el12 = el1 * el2;
        const double el1s = el1 * el1;
        const double el2s = el2 * el2;
        const double s1 = std::sqrt(el1s + gg);
        const double s2 = std::sqrt(el2s + gg);
        const double s12 = s1 * s2;
        const double s1p2 = s1 + s2;
        if (!(s12 > 0.0) || !(s1p2 > 0.0))
            continue;
        const double s = drm * (el1 + el2) / s1p2;
        const double si = -s / s12;

        double elsp = 0.0;
        double sina = 0.0;
        double cosa = 0.0;
        double ratio = 0.0;
        if (el12 > 0.0) {
            const double elgs = -si * s1p2 / (el2 * s1 + el1 * s2);
            elsp = drm * (el1 + el2) * (gg + el1s + el2s) /
                   (el2 * s2 + el1 * s1);
            sina = a * (drm + hm * s12 * elgs);
            cosa = gg + hm * s1p2 + el12 +
                   hh * (gg + el1s + el2s) / (s12 + el12);
            ratio = (el2 > 0.0) ? (s2 + el2) / (s1 + el1)
                                : (s1 - el1) / (s2 - el2);
        } else {
            const double els = el2 * s1 - el1 * s2;
            elsp = el2 * s2 - el1 * s1;
            const double elgs = els / (gg * s12);
            sina = a * (gg * drm + hm * els);
            cosa = (gg + hm * s1) * (gg + hm * s2) + aa * el12;
            ratio = (s1 - el1) * (s2 + el2) / gg;
        }

        if (!(ratio > 0.0))
            return out;

        const double hh113 = signh * std::atan2(sina, cosa);
        const double f111 = std::log(ratio);
        const double f121 = a * ane * f111 + ank * s;
        const double f211 = a * ank * f111 - ane * s;

        b[0] += hh113;
        b[1] += ank * f111;
        b[2] += ane * f111;
        b[3] += a * f111;
        b[4] += ank * f121;
        b[5] += ane * f121;
        b[6] += a * f211;
        b[7] += a * f121;

        if (nf > 6) {
            const double f221 =
                0.5 * ank * ane * (gg + 2.0 * aa) * f111 +
                a * (ank + ane) * (ank - ane) * s - 0.5 * ank * ane * elsp;
            b[8] += ank * f221;
            b[9] += ane * f221;
            if (its != 2) {
                b[10] += 0.5 * a *
                         (3.0 * a * ank * f211 -
                          (aa + ane * ane * hh) * f111 -
                          ane * (a * ank * s - ane * elsp));
                b[11] += 0.5 * a *
                         (3.0 * a * ane * f121 -
                          (aa + ank * ank * hh) * f111 +
                          ank * (a * ane * s + ank * elsp));
                b[12] += a * f221;
            }
        }
    }

    out.source(0, 0) = -b[3] + h * b[0];
    out.source(0, 1) = -0.5 * (b[6] + hh * b[1]);
    out.source(0, 2) = -0.5 * (b[7] + hh * b[2]);
    if (nf > 6 && its != 2) {
        out.source(0, 3) = -(hh * (h * b[0] - b[5]) + b[10]) / 3.0;
        out.source(0, 4) = -(hh * b[4] + b[12]) / 3.0;
        out.source(0, 5) =
            -(hh * (h * b[0] - b[3] + b[5]) + b[11]) / 3.0;
    }

    if (its != 1) {
        out.doublet(0, 0) = b[0];
        out.doublet(0, 1) = -h * b[1];
        out.doublet(0, 2) = -h * b[2];
        out.doublet(0, 3) = h * (b[5] - h * b[0]);
        out.doublet(0, 4) = -h * b[4];
        out.doublet(0, 5) = -h * (b[5] + out.source(0, 0));
        if (nf > 6) {
            out.doublet(0, 6) = h * (b[9] + hh * b[1]);
            out.doublet(0, 7) = -h * (b[8] + out.source(0, 2));
            out.doublet(0, 8) = -h * (b[9] + out.source(0, 1));
            out.doublet(0, 9) = h * (b[8] + hh * b[2]);
        }
    }

    out.source(1, 0) = b[1];
    out.source(1, 1) = h * b[0] - b[5];
    out.source(1, 2) = b[4];
    out.source(2, 0) = b[2];
    out.source(2, 1) = b[4];
    out.source(2, 2) = b[5] + out.source(0, 0);
    out.source(3, 0) = b[0];
    out.source(3, 1) = -h * b[1];
    out.source(3, 2) = -h * b[2];
    if (nf > 6 && its != 2) {
        out.source(1, 3) = -b[9] - hh * b[1];
        out.source(1, 4) = b[8] + out.source(0, 2);
        out.source(1, 5) = b[9] + out.source(0, 1);
        out.source(2, 3) = out.source(1, 4);
        out.source(2, 4) = out.source(1, 5);
        out.source(2, 5) = -b[8] - hh * b[2];
        out.source(3, 3) = -h * out.source(1, 1);
        out.source(3, 4) = -h * b[4];
        out.source(3, 5) = -h * out.source(2, 2);
    }

    if (its != 1) {
        out.doublet(1, 0) = 0.0;
        out.doublet(1, 1) = b[0];
        out.doublet(1, 2) = 0.0;
        out.doublet(1, 3) = 2.0 * out.source(3, 1);
        out.doublet(1, 4) = out.source(3, 2);
        out.doublet(1, 5) = 0.0;
        out.doublet(2, 0) = 0.0;
        out.doublet(2, 1) = 0.0;
        out.doublet(2, 2) = b[0];
        out.doublet(2, 3) = 0.0;
        out.doublet(2, 4) = out.source(3, 1);
        out.doublet(2, 5) = 2.0 * out.source(3, 2);
        out.doublet(3, 0) = 0.0;
        out.doublet(3, 1) = -b[1];
        out.doublet(3, 2) = -b[2];
        out.doublet(3, 3) = -2.0 * out.source(1, 1);
        out.doublet(3, 4) = -2.0 * b[4];
        out.doublet(3, 5) = -2.0 * out.source(2, 2);
        if (nf > 6) {
            out.doublet(1, 6) = 3.0 * out.source(3, 3);
            out.doublet(1, 7) = 2.0 * out.source(3, 4);
            out.doublet(1, 8) = out.source(3, 5);
            out.doublet(1, 9) = 0.0;
            out.doublet(2, 6) = 0.0;
            out.doublet(2, 7) = out.source(3, 3);
            out.doublet(2, 8) = out.doublet(1, 7);
            out.doublet(2, 9) = 3.0 * out.source(3, 5);
            out.doublet(3, 6) = -3.0 * out.source(1, 3);
            out.doublet(3, 7) = -3.0 * out.source(1, 4);
            out.doublet(3, 8) = -3.0 * out.source(2, 4);
            out.doublet(3, 9) = -3.0 * out.source(2, 5);
        }
    }

    const double pi4i = 1.0 / (4.0 * kPi);
    const double pi8i = 0.5 * pi4i;
    const double pi4aj = pi4i * aj;
    const double pi24i = pi8i / 3.0;
    const double x1 = fieldPoint.x();
    const double y1 = fieldPoint.y();
    const double x2 = 0.5 * x1;
    const double y2 = 0.5 * y1;
    const double x3 = x1 / 3.0;
    const double y3 = y1 / 3.0;

    if (its != 2) {
        const int ng = (nf == 10) ? 6 : 3;
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < ng; ++j)
                out.source(i, j) *= pi4aj;
        if (ng == 6) {
            for (int i = 0; i < 4; ++i) {
                out.source(i, 3) =
                    0.5 * out.source(i, 3) +
                    x1 * (x2 * out.source(i, 0) + out.source(i, 1));
                out.source(i, 4) =
                    out.source(i, 4) +
                    x2 * (x1 * out.source(i, 0) + out.source(i, 1)) +
                    x1 * out.source(i, 2);
                out.source(i, 5) =
                    0.5 * out.source(i, 5) +
                    x2 * (y2 * out.source(i, 0) + out.source(i, 2));
            }
        }
        for (int i = 0; i < 4; ++i) {
            out.source(i, 1) += x1 * out.source(i, 0);
            out.source(i, 2) += y1 * out.source(i, 0);
        }
    }

    if (its != 1) {
        for (int i = 0; i < 4; ++i)
            for (int j = 0; j < 6; ++j)
                out.doublet(i, j) *= pi4i;
        for (int i = 0; i < 4; ++i) {
            out.doublet(i, 1) += x1 * out.doublet(i, 0);
            out.doublet(i, 2) += y1 * out.doublet(i, 0);
            const double dvx = out.doublet(i, 1) - x2 * out.doublet(i, 0);
            out.doublet(i, 3) = 0.5 * out.doublet(i, 3) + x1 * dvx;
            const double dvy = out.doublet(i, 2) - y2 * out.doublet(i, 0);
            out.doublet(i, 5) = 0.5 * out.doublet(i, 5) + y1 * dvy;
            out.doublet(i, 4) += x1 * dvy + y1 * dvx;
        }
        if (nf > 6) {
            for (int i = 0; i < 4; ++i) {
                const double dvy = out.doublet(i, 2) - y3 * out.doublet(i, 0);
                const double dvx = out.doublet(i, 1) - x3 * out.doublet(i, 0);
                const double dvxx = out.doublet(i, 3) - x2 * dvx;
                const double dvxy =
                    out.doublet(i, 4) - x2 * dvy - y2 * dvx;
                const double dvyy = out.doublet(i, 5) - y2 * dvy;
                out.doublet(i, 6) = pi24i * out.doublet(i, 6) + x1 * dvxx;
                out.doublet(i, 7) =
                    pi8i * out.doublet(i, 7) + x1 * dvxy + x2 * dvxx;
                out.doublet(i, 8) =
                    pi8i * out.doublet(i, 8) + x1 * dvyy + x2 * dvxy;
                out.doublet(i, 9) = pi24i * out.doublet(i, 9) + y1 * dvyy;
            }
        }
    }

    out.valid = true;
    return out;
}

PanelNearFieldInfluence subsonicPanelNearFieldInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int nf) {
    PanelNearFieldInfluence out;
    if (panel.degenerate || nf != 6)
        return out;

    const auto sourceSplines = buildSourceSubpanelSplines(panel);
    const auto doubletSplines = buildDoubletSubpanelSplines(panel, flow);

    for (int s = 0; s < 8; ++s) {
        const std::size_t si = static_cast<std::size_t>(s);
        const SubpanelGeometry& sp = panel.subpanels[si];
        if (sp.degenerate)
            continue;
        if (!sp.scaledValid ||
            !sourceSplines[si].valid || !doubletSplines[si].valid)
            return PanelNearFieldInfluence{};

        const Eigen::Matrix3d& frame = sp.scaledToLocal;
        // subpqr origin convention: P9 for inner subpanels, first corner
        // for the outer four.
        const Eigen::Vector3d origin = (s >= 4) ? panel.cp[8] : sp.origin;

        std::array<Eigen::Vector3d, 4> vertices{};
        for (int k = 0; k < 3; ++k)
            vertices[static_cast<std::size_t>(k)] =
                frame * (sp.vertex[static_cast<std::size_t>(k)] - origin);
        // All eight subpanels are triangles.  subsbi() uses ics=4 to skip
        // the fourth edge and close the triangle from vertex 3 to vertex 1.
        vertices[3] = vertices[0];

        const Eigen::Vector3d localField = frame * (fieldPoint - origin);
        const SubsonicSubpanelInfluence sub =
            subsonicSubpanelInfluence(vertices, 4, 4, 3, localField,
                                      sp.jacobian, nf);
        if (!sub.valid)
            return PanelNearFieldInfluence{};

        const Eigen::Matrix<double, 4, 3> localSource =
            sub.source.block<4, 3>(0, 0);
        const Eigen::Matrix<double, 4, 6> localDoublet =
            sub.doublet.block<4, 6>(0, 0);

        Eigen::Matrix<double, 4, 5> sourceContribution =
            localSource * sourceSplines[si].coefficients;
        Eigen::Matrix<double, 4, 9> doubletContribution =
            localDoublet * doubletSplines[si].coefficients;

        // Transform the local velocity rows back to global coordinates
        // before accumulating this subpanel.
        for (int col = 0; col < 5; ++col) {
            const Eigen::Vector3d v =
                sp.scaledToGlobal *
                Eigen::Vector3d(sourceContribution(1, col),
                                sourceContribution(2, col),
                                sourceContribution(3, col));
            sourceContribution(1, col) = v.x();
            sourceContribution(2, col) = v.y();
            sourceContribution(3, col) = v.z();
        }
        for (int col = 0; col < 9; ++col) {
            const Eigen::Vector3d v =
                sp.scaledToGlobal *
                Eigen::Vector3d(doubletContribution(1, col),
                                doubletContribution(2, col),
                                doubletContribution(3, col));
            doubletContribution(1, col) = v.x();
            doubletContribution(2, col) = v.y();
            doubletContribution(3, col) = v.z();
        }

        out.source += sourceContribution;
        out.doublet += doubletContribution;
    }

    out.valid = true;
    return out;
}

PanelIntermediateInfluence subsonicPanelIntermediateInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int nf) {
    PanelIntermediateInfluence out;
    if (panel.degenerate || nf != 6)
        return out;

    const SubpanelGeometry& inner = panel.subpanels[4];
    if (!inner.scaledValid)
        return out;

    const Eigen::Matrix3d& frame = inner.scaledToLocal;
    const Eigen::Vector3d origin = panel.cp[8];
    std::array<Eigen::Vector3d, 4> vertices{};
    for (int k = 0; k < 4; ++k)
        vertices[static_cast<std::size_t>(k)] =
            frame * (panel.cp[static_cast<std::size_t>(k)] - origin);

    const Eigen::Vector3d localField = frame * (fieldPoint - origin);
    const SubsonicSubpanelInfluence infl = subsonicSubpanelInfluence(
        vertices, 4, panel.collapsedEdge, 3, localField, inner.jacobian, nf);
    if (!infl.valid)
        return out;

    out.source = infl.source;
    out.doublet = infl.doublet;

    // Local -> global velocity transformation.
    for (int col = 0; col < 6; ++col) {
        const Eigen::Vector3d v =
            inner.scaledToGlobal *
            Eigen::Vector3d(out.source(1, col), out.source(2, col),
                            out.source(3, col));
        out.source(1, col) = v.x();
        out.source(2, col) = v.y();
        out.source(3, col) = v.z();
    }
    for (int col = 0; col < 10; ++col) {
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

} // namespace a502
