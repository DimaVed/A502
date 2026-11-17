// NearField.cpp
// -------------
// Subsonic near-field source/doublet subpanel influence coefficients.  This
// is a direct port of subsbi() from panair.f90.

#include "NearField.h"

#include "Compressibility.h"
#include "PanelGeometry.h"
#include "SubpanelSplines.h"

#include <cmath>
#include <vector>

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

SupersonicSubpanelInfluence supersonicSubpanelInfluence(
    const std::array<Eigen::Vector3d, 4>& vertices, int ns, int ics, int its,
    const Eigen::Vector3d& fieldPoint, double aj, int nf) {
    SupersonicSubpanelInfluence out;
    if (ns < 3 || ns > 4 || (nf != 6 && nf != 10) || its < 1 || its > 3)
        return out;

    // supsbi(): edge-integral form of the linear supersonic kernel.
    constexpr double pi = 3.14159265358979323846;
    constexpr double dltmch = 1.0e2;
    constexpr double deltbt = 1.0e-2;
    const double h = fieldPoint.z();
    const double hh = h * h;

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
        const double bet = (ank - ane) * (ank + ane);
        const double gg = aa - bet * hh;
        const double sbet = std::sqrt(std::abs(bet));
        double el1 = aks1 * ane + aet1 * ank;
        double el2 = aks2 * ane + aet2 * ank;

        double hh113 = 0.0;
        const double rr1 = aks1 * aks1 - aet1 * aet1 - hh;
        const double rr2 = aks2 * aks2 - aet2 * aet2 - hh;
        double r1 = 0.0;
        double r2 = 0.0;
        if (rr1 > 0.0 && aks1 < 0.0)
            r1 = std::sqrt(rr1);
        if (rr2 > 0.0 && aks2 < 0.0)
            r2 = std::sqrt(rr2);

        double f111 = 0.0;
        double f121 = 0.0;
        double f211 = 0.0;
        double f221 = 0.0;

        if (!(r1 > 0.0 || r2 > 0.0)) {
            if (bet <= 0.0 || el1 * el2 >= 0.0 || gg <= 0.0 ||
                a * ank >= 0.0)
                continue;
            if (h != 0.0)
                hh113 = (h * ank >= 0.0) ? pi : -pi;
            f111 = pi / sbet;
            f121 = -a * ane * f111 / bet;
            f211 = a * ank * f111 / bet;
            f221 = -ank * ane * f111 * (gg + 2.0 * aa) /
                   (2.0 * bet * bet);
        } else {
            const double ggb = std::abs(gg);
            const double g = std::sqrt(ggb);
            if (r1 == 0.0)
                el1 = -g;
            if (r2 == 0.0)
                el2 = g;
            double fact1 = 0.0;
            double fact2 = 0.0;
            if (bet > 0.0) {
                fact1 = (el1 * r2 - el2 * r1) / ggb;
                fact2 = (bet * r1 * r2 + el1 * el2) / ggb;
            } else {
                fact1 = (r2 - r1) * (r2 + r1) /
                        (el1 * r2 + el2 * r1);
                fact2 = (ggb - el1 * el1 - el2 * el2) /
                        (bet * r1 * r2 - el1 * el2);
            }
            if (h != 0.0)
                hh113 = std::atan2(h * a * fact1,
                                   r1 * r2 + hh * fact2);

            if (std::abs(fact2) < dltmch * std::abs(sbet * fact1)) {
                if (bet < 0.0)
                    f111 = -((ane >= 0.0) ? 1.0 : -1.0) *
                           std::log((sbet * r1 + std::abs(el1)) /
                                    (sbet * r2 + std::abs(el2))) /
                           sbet;
                else
                    f111 = -std::atan2(sbet * fact1, fact2) / sbet;
                f121 = -(ank * (r2 - r1) + a * ane * f111) / bet;
                f211 = ank * (a * f111 - 2.0 * ane * f121) -
                       ane * (r2 - r1);
                f221 = (3.0 * a * (ank * f121 - ane * f211) +
                        ane * (r2 * aet2 - r1 * aet1) -
                        ank * (r2 * aks2 - r1 * aks1) +
                        2.0 * hh * ank * ane * f111) /
                       (4.0 * bet);
            } else {
                const double sig = fact1 / fact2;
                const double sigs = sig * sig;
                const double seris =
                    sig * sigs *
                    (1.0 / 3.0 - bet * sigs / 5.0 +
                     (bet * sigs) * (bet * sigs) / 7.0);
                f111 = -sig + bet * seris;
                f121 = (el2 * r1 * aet1 - el1 * r2 * aet2 -
                        ank * r1 * r2 * (r2 - r1)) /
                           (ggb * fact2) -
                       a * ane * seris;
                f211 = ank * (a * f111 - 2.0 * ane * f121) -
                       ane * (r2 - r1);
                const double zet1 = ank * aet1 - ane * aks1;
                const double zet2 = ank * aet2 - ane * aks2;
                if (sbet > deltbt) {
                    f221 = (3.0 * a * (ank * f121 - ane * f211) +
                            ane * (r2 * aet2 - r1 * aet1) -
                            ank * (r2 * aks2 - r1 * aks1) +
                            2.0 * hh * ank * ane * f111) /
                           (4.0 * bet);
                } else {
                    f221 = -0.5 * a * (r2 - r1) +
                           a * (hh * hh * (r2 - r1) -
                                2.0 * hh * ank * ane * a *
                                    (zet2 * r2 - zet1 * r1) +
                                1.5 * aa * (zet2 * zet2 * r2 -
                                            zet1 * zet1 * r1)) /
                               (15.0 * ggb * ggb);
                }
            }
        }

        b[0] += hh113;
        b[1] += ank * f111;
        b[2] += ane * f111;
        b[3] += a * f111;
        b[4] += ank * f121;
        b[5] += ane * f121;
        b[6] += a * f211;
        b[7] += a * f121;
        if (nf > 6) {
            b[8] += ank * f221;
            b[9] += ane * f221;
            if (its != 2) {
                const double fact1 =
                    (aa + 0.5 * std::abs(gg)) * f111 -
                    4.0 * ank * ane * f221 -
                    0.5 * (el2 * r2 - el1 * r1);
                const double fact2 = a * (ank * f211 - ane * f121);
                b[10] += a * (ane * ane * fact1 + fact2);
                b[11] += a * (ank * ank * fact1 - fact2);
                b[12] += a * f221;
            }
        }
    }

    out.source(0, 0) = -b[3] - h * b[0];
    out.source(0, 1) = -0.5 * (b[6] - hh * b[1]);
    out.source(0, 2) = -0.5 * (b[7] + hh * b[2]);
    if (nf > 6 && its != 2) {
        out.source(0, 3) = -(hh * (b[5] + h * b[0]) + b[10]) / 3.0;
        out.source(0, 4) = (hh * b[4] - b[12]) / 3.0;
        out.source(0, 5) =
            -(hh * (b[5] - h * b[0] - b[3]) + b[11]) / 3.0;
    }

    if (its != 1) {
        out.doublet(0, 0) = -b[0];
        out.doublet(0, 1) = h * b[1];
        out.doublet(0, 2) = -h * b[2];
        out.doublet(0, 3) = -h * (b[5] + h * b[0]);
        out.doublet(0, 4) = h * b[4];
        out.doublet(0, 5) = -h * (b[5] + out.source(0, 0));
        if (nf > 6) {
            out.doublet(0, 6) = h * (hh * b[1] - b[9]);
            out.doublet(0, 7) = h * (b[8] + out.source(0, 2));
            out.doublet(0, 8) = -h * (b[9] + out.source(0, 1));
            out.doublet(0, 9) = h * (b[8] + hh * b[2]);
        }
    }

    out.source(1, 0) = b[1];
    out.source(1, 1) = -b[5] - h * b[0];
    out.source(1, 2) = b[4];
    out.source(2, 0) = b[2];
    out.source(2, 1) = -b[4];
    out.source(2, 2) = b[5] + out.source(0, 0);
    out.source(3, 0) = -b[0];
    out.source(3, 1) = h * b[1];
    out.source(3, 2) = -h * b[2];
    if (nf > 6) {
        out.source(1, 3) = hh * b[1] - b[9];
        out.source(1, 4) = b[8] + out.source(0, 2);
        out.source(1, 5) = -b[9] - out.source(0, 1);
        out.source(2, 3) = -b[8] - out.source(0, 2);
        out.source(2, 4) = -out.source(1, 5);
        out.source(2, 5) = -b[8] - hh * b[2];
        out.source(3, 3) = h * out.source(1, 1);
        out.source(3, 4) = h * b[4];
        out.source(3, 5) = -h * out.source(2, 2);
    }

    if (its != 1) {
        out.doublet(1, 0) = 0.0;
        out.doublet(1, 1) = -b[0];
        out.doublet(1, 2) = 0.0;
        out.doublet(1, 3) = 2.0 * out.source(3, 1);
        out.doublet(1, 4) = out.source(3, 2);
        out.doublet(1, 5) = 0.0;
        out.doublet(2, 0) = 0.0;
        out.doublet(2, 1) = 0.0;
        out.doublet(2, 2) = -b[0];
        out.doublet(2, 3) = 0.0;
        out.doublet(2, 4) = out.source(3, 1);
        out.doublet(2, 5) = 2.0 * out.source(3, 2);
        out.doublet(3, 0) = 0.0;
        out.doublet(3, 1) = b[1];
        out.doublet(3, 2) = -b[2];
        out.doublet(3, 3) = 2.0 * out.source(1, 1);
        out.doublet(3, 4) = 2.0 * b[4];
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

    const double pi4i = 1.0 / (4.0 * pi);
    const double pi2i = 2.0 * pi4i;
    const double pi2aj = pi2i * aj;
    const double pi4aj = pi4i * aj;
    const double pi12i = pi4i / 3.0;
    const double x1 = fieldPoint.x();
    const double y1 = fieldPoint.y();
    const double x2 = 0.5 * x1;
    const double y2 = 0.5 * y1;
    const double x3 = x1 / 3.0;
    const double y3 = y1 / 3.0;

    for (int i = 0; i < 4; ++i) {
        if (its != 2) {
            const int ng = (nf > 6) ? 6 : 3;
            for (int j = 0; j < ng; ++j)
                out.source(i, j) *= pi2aj;
            out.source(i, 1) += x1 * out.source(i, 0);
            out.source(i, 2) += y1 * out.source(i, 0);
            if (ng == 6) {
                const double dux = out.source(i, 1) - x2 * out.source(i, 0);
                const double duy = out.source(i, 2) - y2 * out.source(i, 0);
                out.source(i, 3) =
                    pi4aj * out.source(i, 3) + x1 * dux;
                out.source(i, 4) =
                    pi2aj * out.source(i, 4) + x1 * duy + y1 * dux;
                out.source(i, 5) =
                    pi4aj * out.source(i, 5) + y1 * duy;
            }
        }
        if (its != 1) {
            for (int j = 0; j < 6; ++j)
                out.doublet(i, j) *= pi2i;
            out.doublet(i, 1) += x1 * out.doublet(i, 0);
            out.doublet(i, 2) += y1 * out.doublet(i, 0);
            const double dvx = out.doublet(i, 1) - x2 * out.doublet(i, 0);
            const double dvy = out.doublet(i, 2) - y2 * out.doublet(i, 0);
            out.doublet(i, 3) =
                pi4i * out.doublet(i, 3) + x1 * dvx;
            out.doublet(i, 4) =
                pi2i * out.doublet(i, 4) + x1 * dvy + y1 * dvx;
            out.doublet(i, 5) =
                pi4i * out.doublet(i, 5) + y1 * dvy;
            if (nf > 6) {
                const double dvx3 =
                    out.doublet(i, 1) - x3 * out.doublet(i, 0);
                const double dvy3 =
                    out.doublet(i, 2) - y3 * out.doublet(i, 0);
                const double dvxx = out.doublet(i, 3) - x2 * dvx3;
                const double dvxy =
                    out.doublet(i, 4) - x2 * dvy3 - y2 * dvx3;
                const double dvyy = out.doublet(i, 5) - y2 * dvy3;
                out.doublet(i, 6) =
                    pi12i * out.doublet(i, 6) + x1 * dvxx;
                out.doublet(i, 7) =
                    pi4i * out.doublet(i, 7) + x1 * dvxy + y1 * dvxx;
                out.doublet(i, 8) =
                    pi4i * out.doublet(i, 8) + x1 * dvyy + y1 * dvxy;
                out.doublet(i, 9) =
                    pi12i * out.doublet(i, 9) + y1 * dvyy;
            }
        }
    }

    out.valid = true;
    return out;
}

namespace {

// ---------------------------------------------------------------------------
// aicsup()/shftic(): closed-form supersonic superinclined subpanel AICs.
// ---------------------------------------------------------------------------

namespace {

struct AicSupResult {
    Eigen::Matrix<double, 4, 6> source{
        Eigen::Matrix<double, 4, 6>::Zero()};
    Eigen::Matrix<double, 4, 10> doublet{
        Eigen::Matrix<double, 4, 10>::Zero()};
    bool valid{true};
};

// shftic(): translate polynomial coefficients from the Mach-cone-centred
// basis to the panel origin.  Columns are ordered as
// [1, x, y, x^2/2, xy, y^2/2, x^3/6, x^2*y/2, x*y^2/2, y^3/6].
template <int Cols>
void shiftAic(Eigen::Matrix<double, 4, Cols>& z, int degree,
              double x, double y) {
    for (int k = 1; k <= degree; ++k) {
        for (int lbk = 1; lbk <= k; ++lbk) {
            const int l = k + 1 - lbk;
            const int lbs = l * (l + 1) / 2;
            const int lbsml = lbs - l;
            const double f = static_cast<double>(degree + 1 - l);
            const double xf = x / f;
            const double yf = y / f;
            for (int i = 1; i <= l; ++i) {
                const int iplbs = i + lbs;
                const int ilbsml = i + lbsml;
                for (int iz = 0; iz < 4; ++iz) {
                    z(iz, iplbs - 1) += xf * z(iz, ilbsml - 1);
                    z(iz, iplbs) += yf * z(iz, ilbsml - 1);
                }
            }
        }
    }
}

AicSupResult computeAicSup(const std::vector<std::array<double, 2>>& q,
                           int n, const std::array<double, 2>& p, double x,
                           int nf, bool doublt, double btinv) {
    AicSupResult out;
    constexpr double pi = 3.14159265358979323846;
    const double pi2 = 2.0 * pi;
    const int ng = 3 * (1 + nf / 10);
    if (x < 0.0 || n < 3)
        return out;
    if (n > 20) {
        out.valid = false;
        return out;
    }

    const bool xeqzro = x == 0.0;
    const double xsq = x * x;

    // Signed panel area in the y-z (cone) plane.
    double area = 0.0;
    int km1 = n;
    for (int k = 1; k <= n; ++k) {
        area += q[static_cast<std::size_t>(km1)][0] *
                    q[static_cast<std::size_t>(k)][1] -
                q[static_cast<std::size_t>(km1)][1] *
                    q[static_cast<std::size_t>(k)][0];
        km1 = km1 % n + 1;
    }
    if (!(area > 0.0))
        return out;

    std::vector<std::array<double, 2>> rho(
        static_cast<std::size_t>(n + 1));
    std::vector<std::array<double, 2>> tg(
        static_cast<std::size_t>(n + 1));
    std::vector<std::array<double, 2>> nm(
        static_cast<std::size_t>(n + 1));
    std::vector<double> rhosq(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<double> d(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<double> vm(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<double> vp(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<double> r(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<double> phi(static_cast<std::size_t>(n + 1), 0.0);
    std::vector<char> edge(static_cast<std::size_t>(n + 1), 0);
    std::vector<char> corner(static_cast<std::size_t>(n + 1), 0);

    km1 = n;
    for (int k = 1; k <= n; ++k) {
        const std::size_t ks = static_cast<std::size_t>(k);
        const std::size_t kms = static_cast<std::size_t>(km1);
        rho[ks][0] = q[ks][0] - p[0];
        rho[ks][1] = q[ks][1] - p[1];
        rhosq[ks] = rho[ks][0] * rho[ks][0] + rho[ks][1] * rho[ks][1];
        corner[ks] = (!xeqzro && rhosq[ks] < xsq) ? 1 : 0;

        const double dy = q[ks][0] - q[kms][0];
        const double dz = q[ks][1] - q[kms][1];
        const double dnm = std::sqrt(dy * dy + dz * dz);
        if (!(dnm > 0.0)) {
            out.valid = false;
            return out;
        }
        tg[ks][0] = dy / dnm;
        tg[ks][1] = dz / dnm;
        nm[ks][0] = tg[ks][1];
        nm[ks][1] = -tg[ks][0];
        d[ks] = rho[ks][0] * nm[ks][0] + rho[ks][1] * nm[ks][1];
        km1 = km1 % n + 1;
    }

    bool intsct = false;
    bool within = true;
    bool convex = true;
    bool subset = true;
    area = 0.0;
    km1 = n;
    int ic2 = 2;
    for (int k = 1; k <= n; ++k) {
        const std::size_t ks = static_cast<std::size_t>(k);
        const std::size_t kms = static_cast<std::size_t>(km1);
        subset = subset && (corner[ks] != 0);
        vp[ks] = tg[ks][0] * rho[ks][0] + tg[ks][1] * rho[ks][1];
        vm[ks] = tg[ks][0] * rho[kms][0] + tg[ks][1] * rho[kms][1];
        edge[ks] = (corner[kms] || corner[ks]) ? 1 : 0;
        const double vpvm = vp[ks] * vm[ks];
        if (std::abs(d[ks]) < x && x != 0.0 && vpvm <= 0.0)
            edge[ks] = 1;
        intsct = intsct || edge[ks] != 0;

        const double areak = rho[kms][0] * rho[ks][1] -
                             rho[kms][1] * rho[ks][0];
        area += areak;

        const double snalfa = tg[kms][0] * nm[ks][0] +
                              tg[kms][1] * nm[ks][1];
        if (snalfa < 0.0)
            convex = false;
        const int isgn = (snalfa < 0.0) ? -1 : 1;
        const int ik = (d[ks] < 0.0) ? -1 : 1;
        const int ikm1 = (d[kms] < 0.0) ? -1 : 1;
        if (ik * ikm1 < 0)
            ic2 -= isgn;
        km1 = km1 % n + 1;
    }
    within = ic2 == 2;
    if (!convex && !subset) {
        // A502 aborts on this case; report it as an unusable influence.
        out.valid = false;
        return out;
    }
    if (area <= 0.0)
        return out;
    if (!intsct && !within)
        return out;

    // Edge functions r(k) and panel function psi.
    int kp1 = 2;
    double psi = pi2;
    for (int k = 1; k <= n; ++k) {
        const std::size_t ks = static_cast<std::size_t>(k);
        const std::size_t kp1s = static_cast<std::size_t>(kp1);
        r[ks] = 0.0;
        if (corner[ks])
            psi += pi;
        if (edge[ks])
            psi -= pi;
        if (corner[ks]) {
            r[ks] = std::sqrt(xsq - rhosq[ks]);
            const double tktkp1 = tg[ks][0] * tg[kp1s][0] +
                                  tg[ks][1] * tg[kp1s][1];
            const double tcross = tg[ks][0] * tg[kp1s][1] -
                                  tg[ks][1] * tg[kp1s][0];
            const double xval = d[ks] * d[kp1s] - xsq * tktkp1;
            const double yval = x * r[ks] * tcross;
            psi -= std::atan2(yval, -xval);
        }
        kp1 = kp1 % n + 1;
    }

    // Edge function phi(k).
    km1 = n;
    for (int k = 1; k <= n; ++k) {
        const std::size_t ks = static_cast<std::size_t>(k);
        const std::size_t kms = static_cast<std::size_t>(km1);
        phi[ks] = 0.0;
        if (edge[ks]) {
            phi[ks] = -pi;
            if (corner[ks] || corner[kms]) {
                double vpk = vp[ks];
                if (!corner[ks])
                    vpk = 1.0;
                double vmk = vm[ks];
                if (!corner[kms])
                    vmk = -1.0;
                const double xval = vpk * vmk + r[ks] * r[kms];
                const double yval = r[ks] * vmk - r[kms] * vpk;
                phi[ks] = std::atan2(yval, xval);
            }
        }
        km1 = km1 % n + 1;
    }

    Eigen::Matrix<double, 4, 6>& saic = out.source;
    Eigen::Matrix<double, 4, 10>& daic = out.doublet;

    const double psix = psi * x;
    const double psixx2 = 0.5 * x * psix;
    saic(0, 0) = psi;
    saic(1, 1) = psix;
    saic(2, 2) = psix;
    saic(3, 0) = psix;
    if (doublt) {
        daic(0, 3) = psix;
        daic(0, 5) = psix;
        daic(1, 1) = psi;
        daic(2, 2) = psi;
        daic(3, 0) = psi;
        daic(3, 3) = psixx2;
        daic(3, 5) = psixx2;
    }

    double h[2][2][2] = {};
    double bv[2] = {0.0, 0.0};
    double btw[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
    double aipsrr[2][2] = {{0.0, 0.0}, {0.0, 0.0}};

    if (intsct) {
        km1 = n;
        for (int k = 1; k <= n; ++k) {
            const std::size_t ks = static_cast<std::size_t>(k);
            const std::size_t kms = static_cast<std::size_t>(km1);
            if (edge[ks]) {
                const double phdmns = -phi[ks] * d[ks];
                const double phny = phi[ks] * nm[ks][0];
                const double phnz = phi[ks] * nm[ks][1];
                const double phvsqh =
                    0.5 * phi[ks] * (xsq - d[ks] * d[ks]);
                const double delr = r[ks] - r[kms];
                const double drdmn2 = -delr * d[ks] * 0.5;

                saic(0, 1) += phny * x;
                saic(0, 2) += phnz * x;
                saic(1, 0) += phny;
                saic(1, 1) += phdmns * tg[ks][0] * tg[ks][0];
                saic(1, 2) += phdmns * tg[ks][0] * tg[ks][1];
                saic(2, 0) += phnz;
                saic(2, 1) += phdmns * tg[ks][1] * tg[ks][0];
                saic(2, 2) += phdmns * tg[ks][1] * tg[ks][1];
                saic(3, 0) += phdmns;
                saic(3, 1) += phvsqh * nm[ks][0];
                saic(3, 2) += phvsqh * nm[ks][1];
                saic(1, 1) += delr * tg[ks][0] * nm[ks][0];
                saic(1, 2) += delr * tg[ks][0] * nm[ks][1];
                saic(2, 1) += delr * tg[ks][1] * nm[ks][0];
                saic(2, 2) += delr * tg[ks][1] * nm[ks][1];
                saic(3, 1) += drdmn2 * tg[ks][0];
                saic(3, 2) += drdmn2 * tg[ks][1];

                const double phdxm2 = phdmns * x * 0.5;
                const double delrx2 = delr * x * 0.5;
                const double tyy = tg[ks][0] * tg[ks][0];
                const double tyz = tg[ks][0] * tg[ks][1];
                const double tyz2 = 2.0 * tyz;
                const double tzz = tg[ks][1] * tg[ks][1];
                const double tskw = tzz - tyy;
                if (doublt) {
                    daic(0, 1) += phny;
                    daic(0, 2) += phnz;
                    daic(0, 3) += phdmns * tyy;
                    daic(0, 4) += phdmns * tyz2;
                    daic(0, 5) += phdmns * tzz;
                    daic(1, 3) += phny * x;
                    daic(1, 4) += phnz * x;
                    daic(2, 4) += phny * x;
                    daic(2, 5) += phnz * x;
                    daic(3, 1) += phny * x;
                    daic(3, 2) += phnz * x;
                    daic(3, 3) += phdxm2 * tyy;
                    daic(3, 4) += phdxm2 * tyz2;
                    daic(3, 5) += phdxm2 * tzz;
                    daic(0, 3) += delr * tyz;
                    daic(0, 4) += delr * tskw;
                    daic(0, 5) += -delr * tyz;
                    daic(3, 3) += delrx2 * tyz;
                    daic(3, 4) += delrx2 * tskw;
                    daic(3, 5) += -delrx2 * tyz;
                }

                if (nf >= 10) {
                    const double cnn = -phi[ks] * d[ks] * d[ks];
                    const double cnt = -d[ks] * delr;
                    const double ctt =
                        -phvsqh -
                        0.5 * (r[ks] * vp[ks] - r[kms] * vm[ks]);
                    btw[0][0] = cnn * tzz + cnt * tyz2 + ctt * tyy;
                    btw[0][1] = -cnn * tyz + cnt * tskw + ctt * tyz;
                    btw[1][1] = cnn * tyy - cnt * tyz2 + ctt * tzz;
                    btw[1][0] = btw[0][1];
                    for (int l = 0; l < 2; ++l)
                        for (int i = 0; i < 2; ++i)
                            for (int j = 0; j < 2; ++j)
                                h[i][j][l] += btw[i][j] * nm[ks][l];
                    aipsrr[0][0] += d[ks] * btw[0][0];
                    aipsrr[0][1] += d[ks] * btw[0][1];
                    aipsrr[1][1] += d[ks] * btw[1][1];
                }
            }
            km1 = km1 % n + 1;
        }
    }

    if (nf >= 10) {
        btw[0][0] = -saic(1, 1);
        btw[0][1] = -saic(1, 2);
        btw[1][0] = -saic(2, 1);
        btw[1][1] = -saic(2, 2);
        bv[0] = saic(3, 1);
        bv[1] = saic(3, 2);
        for (int ik = 0; ik < 2; ++ik)
            for (int j = 0; j < 2; ++j)
                h[ik][j][ik] -= bv[j];
        for (int jk = 0; jk < 2; ++jk)
            for (int i = 0; i < 2; ++i)
                h[i][jk][jk] -= bv[i];

        if (doublt) {
            daic(0, 6) = -0.5 * h[0][0][0];
            daic(0, 7) =
                -0.5 * (h[0][0][1] + h[0][1][0] + h[1][0][0]);
            daic(0, 8) =
                -0.5 * (h[1][1][0] + h[1][0][1] + h[0][1][1]);
            daic(0, 9) = -0.5 * h[1][1][1];
            daic(1, 6) = -0.5 * x * btw[0][0];
            daic(1, 7) = -0.5 * x * (btw[0][1] + btw[1][0]);
            daic(1, 8) = -0.5 * x * btw[1][1];
            for (int j = 7; j <= 9; ++j)
                daic(2, j) = daic(1, j - 1);
            for (int j = 6; j <= 9; ++j)
                daic(3, j) = (x / 3.0) * daic(0, j);
        }

        saic(0, 3) = 0.5 * x * saic(1, 1);
        saic(0, 4) = x * saic(1, 2);
        saic(0, 5) = 0.5 * x * saic(2, 2);
        saic(1, 3) = -0.5 * h[0][0][0];
        saic(1, 4) = -h[0][0][1];
        saic(1, 5) = -0.5 * h[0][1][1];
        saic(2, 3) = -0.5 * h[0][0][1];
        saic(2, 4) = -h[0][1][1];
        saic(2, 5) = -0.5 * h[1][1][1];
        saic(3, 3) =
            (1.0 / 6.0) * (aipsrr[0][0] - xsq * btw[0][0]);
        saic(3, 4) =
            (1.0 / 3.0) * (aipsrr[0][1] - xsq * btw[0][1]);
        saic(3, 5) =
            (1.0 / 6.0) * (aipsrr[1][1] - xsq * btw[1][1]);
    }

    double fn = 1.0 / pi2;
    if (xeqzro && within)
        fn *= 0.5;
    const double sn = -fn * btinv;
    saic *= sn;
    daic *= fn;

    if (nf > 6) {
        shiftAic(saic, 2, p[0], p[1]);
        if (doublt)
            shiftAic(daic, 3, p[0], p[1]);
    } else {
        shiftAic(saic, 1, p[0], p[1]);
        if (doublt)
            shiftAic(daic, 2, p[0], p[1]);
    }

    (void)ng;
    return out;
}

} // namespace

using SubpanelKernel = SubsonicSubpanelInfluence (*)(
    const std::array<Eigen::Vector3d, 4>&, int, int, int,
    const Eigen::Vector3d&, double, int);

PanelNearFieldInfluence panelNearFieldInfluenceInternal(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int nf, bool selectByMach) {
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

        // nftpic() dispatch: subsonic flows always use subsbi(); supersonic
        // flows use supsbi() for subinclined and supspi() for superinclined
        // subpanels.
        SubpanelKernel kernel = &subsonicSubpanelInfluence;
        bool superinclined = false;
        if (selectByMach && flow.mach > 1.0) {
            if (sp.inclination > 0) {
                kernel = &supersonicSubpanelInfluence;
            } else if (sp.inclination < 0) {
                superinclined = true;
            } else {
                return PanelNearFieldInfluence{};
            }
        }

        const Eigen::Vector3d localField = frame * (fieldPoint - origin);
        SubsonicSubpanelInfluence sub;
        if (superinclined) {
            const double sgnx =
                (flow.compd.dot(sp.normal) >= 0.0) ? 1.0 : -1.0;
            sub = supersonicSuperinclinedSubpanelInfluence(
                vertices, 4, 4, 3, localField, sp.jacobian, nf, sgnx);
        } else {
            sub = kernel(vertices, 4, 4, 3, localField, sp.jacobian, nf);
        }
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

} // namespace

SupersonicSubpanelInfluence supersonicSuperinclinedSubpanelInfluence(
    const std::array<Eigen::Vector3d, 4>& vertices, int ns, int ics, int its,
    const Eigen::Vector3d& fieldPoint, double aj, int nf, double sgnx) {
    SupersonicSubpanelInfluence out;
    if (ns < 3 || ns > 4 || (nf != 6 && nf != 10) || its < 1 || its > 3)
        return out;

    std::vector<std::array<double, 2>> q(
        static_cast<std::size_t>(ns + 1));
    int n = 0;
    for (int i = 1; i <= ns; ++i) {
        if (i == ics)
            continue;
        ++n;
        q[static_cast<std::size_t>(n)] = {vertices[static_cast<std::size_t>(
                                               i - 1)]
                                              .x(),
                                          vertices[static_cast<std::size_t>(
                                               i - 1)]
                                              .y()};
    }
    const std::array<double, 2> p = {fieldPoint.x(), fieldPoint.y()};
    const double x = sgnx * fieldPoint.z();
    const bool doublt = its >= 2;
    const AicSupResult aic = computeAicSup(q, n, p, x, nf, doublt, aj);
    if (!aic.valid)
        return out;

    // supspi(): arrange the aicsup rows into nftpic's (ph,vx,vy,vz) layout.
    const int ng = 3 * (1 + nf / 10);
    for (int i = 0; i < ng; ++i) {
        out.source(0, i) = aic.source(3, i);
        out.source(1, i) = aic.source(1, i);
        out.source(2, i) = aic.source(2, i);
        out.source(3, i) = aic.source(0, i) * sgnx;
    }
    const int nd = doublt ? nf : 0;
    for (int i = 0; i < nd; ++i) {
        out.doublet(0, i) = aic.doublet(3, i) * sgnx;
        out.doublet(1, i) = aic.doublet(1, i) * sgnx;
        out.doublet(2, i) = aic.doublet(2, i) * sgnx;
        out.doublet(3, i) = aic.doublet(0, i);
    }
    out.valid = true;
    return out;
}

PanelNearFieldInfluence subsonicPanelNearFieldInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int nf) {
    return panelNearFieldInfluenceInternal(panel, flow, fieldPoint, nf,
                                           false);
}

PanelNearFieldInfluence panelNearFieldInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int nf) {
    return panelNearFieldInfluenceInternal(panel, flow, fieldPoint, nf, true);
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
