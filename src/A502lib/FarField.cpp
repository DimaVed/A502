// FarField.cpp
// ------------
// Far-field AIC assembly (ffpic/ffdqgv) and canonical-value mapping.

#include "FarField.h"

#include "Compressibility.h"
#include "PanelGeometry.h"
#include "PanelParameters.h"
#include "SubpanelSplines.h"

#include <algorithm>
#include <cmath>

namespace a502 {

namespace {

constexpr double kPi = 3.14159265358979323846;

constexpr int kMi[10] = {0, 1, 0, 2, 1, 0, 3, 2, 1, 0};
constexpr int kNj[10] = {0, 0, 1, 0, 1, 2, 0, 1, 2, 3};
constexpr double kFcc[10] = {1.0, 1.0, 1.0, 2.0, 1.0,
                             2.0, 6.0, 2.0, 2.0, 6.0};

} // namespace

FarFieldPanelData buildFarFieldPanelData(
    const PanelGeometry& panel, const Compressibility& flow, int its) {
    FarFieldPanelData out;
    if (panel.degenerate)
        return out;
    const SubpanelGeometry& inner = panel.subpanels[4];
    if (!inner.scaledValid)
        return out;

    out.center = panel.cp[8];
    out.toLocal = inner.scaledToLocal;
    out.toGlobal = inner.scaledToGlobal;
    out.ajf = inner.jacobian;
    out.sf = flow.sbetam;
    out.rf = static_cast<double>(inner.inclination);
    out.its = its;

    std::array<Eigen::Vector3d, 4> corners;
    for (int k = 0; k < 4; ++k)
        corners[static_cast<std::size_t>(k)] =
            out.toLocal *
            (panel.cp[static_cast<std::size_t>(k)] - out.center);

    out.cf = panelAreaMoments(corners, 4, panel.collapsedEdge, 5);
    if (!(out.cf(0, 0) > 0.0))
        return out;

    const double pif = 1.0 / (kPi * (3.0 + out.sf));
    const double fh = out.sf * pif;
    for (int ka = 0; ka < 10; ++ka) {
        for (int kb = 0; kb < 6; ++kb) {
            const int i = kMi[ka] + kMi[kb];
            const int j = kNj[ka] + kNj[kb];
            out.hm(ka, kb) =
                fh * out.cf(i, j) / (kFcc[ka] * kFcc[kb]);
        }
    }

    for (int ka = 0; ka < 6; ++ka) {
        out.hbm[static_cast<std::size_t>(ka)](0, 0) = out.hm(ka, 1);
        out.hbm[static_cast<std::size_t>(ka)](0, 1) =
            2.0 * out.hm(ka, 3);
        out.hbm[static_cast<std::size_t>(ka)](0, 2) = out.hm(ka, 4);
        out.hbm[static_cast<std::size_t>(ka)](1, 0) = out.hm(ka, 2);
        out.hbm[static_cast<std::size_t>(ka)](1, 1) = out.hm(ka, 4);
        out.hbm[static_cast<std::size_t>(ka)](1, 2) =
            2.0 * out.hm(ka, 5);
    }

    out.valid = true;
    return out;
}

FarFieldAic evaluateFarFieldAic(
    const FarFieldPanelData& data, const Eigen::Vector3d& fieldPoint,
    int iflu, int nf) {
    FarFieldAic out;
    if (!data.valid || iflu < 1 || iflu > 3 || (nf != 6 && nf != 10))
        return out;

    const int itsx = data.its;
    const int itsf = data.its;
    const int ng = 3 * (1 + nf / 10);
    constexpr int nt[4] = {0, 1, 3, 6};

    const Eigen::Vector3d pc = fieldPoint - data.center;
    const Eigen::Vector3d local = data.toLocal * pc;
    const double x = local.x();
    const double y = local.y();
    const double z = local.z();

    const double denom =
        std::sqrt(data.rf * x * x + data.sf * y * y +
                  data.rf * data.sf * z * z);
    if (!(denom > 0.0))
        return out;
    const double p1 = 1.0 / denom;
    const double p2 = p1 * p1;
    const double p3 = p2 * p1;
    const double fs = -data.ajf * data.sf;
    const double yp = y * p1 * data.sf;
    const double xp = x * p1 * data.rf;

    double gk[7] = {};
    double hk[7] = {};
    if (iflu >= 3) {
        const double p5t3 = 3.0 * p2 * p3;
        hk[3] = p5t3 * (5.0 * xp * xp - data.rf);
        hk[4] = p5t3 * 5.0 * xp * yp;
        hk[5] = p5t3 * (5.0 * yp * yp - data.sf);
        if (itsx != 2) {
            gk[3] = fs * p3 * (3.0 * xp * xp - data.rf);
            gk[4] = fs * p3 * 3.0 * xp * yp;
            gk[5] = fs * p3 * (3.0 * yp * yp - data.sf);
        }
    }
    if (iflu >= 2) {
        const double p4t3 = 3.0 * p2 * p2;
        hk[1] = p4t3 * xp;
        hk[2] = p4t3 * yp;
        if (itsx != 2) {
            gk[1] = fs * p2 * xp;
            gk[2] = fs * p2 * yp;
        }
    }
    hk[0] = p3;
    gk[0] = fs * p1;

    const int nmom = nt[iflu];
    double phsdsn[6] = {};
    if (itsx != 2) {
        for (int j = 0; j < ng; ++j) {
            double sum = 0.0;
            for (int k = 0; k < nmom; ++k)
                sum += data.hm(j, k) * gk[k];
            phsdsn[j] = sum;
        }
    }

    const int nh = (itsx == 1) ? ng : std::max(nf, ng);
    double h[10] = {};
    for (int j = 0; j < nh; ++j) {
        double sum = 0.0;
        for (int k = 0; k < nmom; ++k)
            sum += data.hm(j, k) * hk[k];
        h[j] = sum;
    }

    double hb[2][6] = {};
    if (iflu > 1) {
        const int nmomHb = nt[iflu - 1];
        for (int j = 0; j < ng; ++j) {
            double s0 = 0.0;
            double s1 = 0.0;
            const Eigen::Matrix<double, 2, 3>& hbmj =
                data.hbm[static_cast<std::size_t>(j)];
            for (int k = 0; k < nmomHb; ++k) {
                s0 += hbmj(0, k) * hk[k];
                s1 += hbmj(1, k) * hk[k];
            }
            hb[0][j] = s0;
            hb[1][j] = s1;
        }
    }
    for (int j = 0; j < ng; ++j) {
        hb[0][j] -= x * h[j];
        hb[1][j] -= y * h[j];
    }

    if (itsx != 2) {
        const double ajrs = data.ajf * data.rf * data.sf;
        const double ajrz = data.ajf * data.rf * z;
        for (int j = 0; j < ng; ++j) {
            const Eigen::Vector3d vs(
                -ajrs * hb[0][j], -data.ajf * hb[1][j],
                ajrz * h[j]);
            const Eigen::Vector3d global = data.toGlobal * vs;
            out.source(0, j) = phsdsn[j];
            out.source(1, j) = global.x();
            out.source(2, j) = global.y();
            out.source(3, j) = global.z();
        }
        if (itsf == 1) {
            out.valid = true;
            return out;
        }
    }

    if (itsx == 1) {
        out.valid = true;
        return out;
    }

    double zh[10] = {};
    for (int j = 0; j < nf; ++j)
        zh[j] = z * h[j];

    for (int j = 0; j < 6; ++j)
        out.doublet(0, j) = zh[j];
    const double row1[6] = {0.0, zh[0], 0.0, zh[1], zh[2], 0.0};
    const double row2[6] = {0.0, 0.0, zh[0], 0.0, zh[1], zh[2]};
    const double row3[6] = {0.0,
                            hb[0][0],
                            hb[1][0],
                            hb[0][1],
                            hb[0][2] + hb[1][1],
                            hb[1][2]};
    for (int j = 0; j < 6; ++j) {
        out.doublet(1, j) = row1[j];
        out.doublet(2, j) = row2[j];
        out.doublet(3, j) = row3[j];
    }

    if (nf > 6) {
        for (int j = 6; j < 10; ++j)
            out.doublet(0, j) = zh[j];
        out.doublet(1, 6) = zh[3];
        out.doublet(1, 7) = zh[4];
        out.doublet(1, 8) = zh[5];
        out.doublet(1, 9) = 0.0;
        out.doublet(2, 6) = 0.0;
        out.doublet(2, 7) = zh[3];
        out.doublet(2, 8) = zh[4];
        out.doublet(2, 9) = zh[5];
        out.doublet(3, 6) = hb[0][3];
        out.doublet(3, 7) = hb[0][4] + hb[1][3];
        out.doublet(3, 8) = hb[0][5] + hb[1][4];
        out.doublet(3, 9) = hb[1][5];
    }

    for (int j = 0; j < nf; ++j) {
        const Eigen::Vector3d v = data.toGlobal *
                                  Eigen::Vector3d(out.doublet(1, j),
                                                  out.doublet(2, j),
                                                  out.doublet(3, j));
        out.doublet(1, j) = v.x();
        out.doublet(2, j) = v.y();
        out.doublet(3, j) = v.z();
    }

    out.valid = true;
    return out;
}

PanelFarFieldInfluence farFieldPanelInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int iflu, int nf, int its) {
    PanelFarFieldInfluence out;
    if (nf != 6)
        return out;

    const FarFieldPanelData data =
        buildFarFieldPanelData(panel, flow, its);
    const FarFieldAic aic =
        evaluateFarFieldAic(data, fieldPoint, iflu, nf);
    if (!aic.valid)
        return out;

    const auto sourceSplines = buildSourceSubpanelSplines(panel);
    if (!sourceSplines[4].valid)
        return out;

    const NearPlaneQuadratic qa =
        buildNearPlaneQuadratic(panel.cp, data.toLocal);
    if (!qa.valid)
        return out;

    // ffpic gives AICs for the far-field Taylor basis.  The source basis is
    // the linear distribution fitted by the outer spline, which the inner
    // subpanel spline reproduces exactly.  The doublet basis is the
    // near-plane quadratic qa of qcof(), as in the A502 ffpic caller.
    out.source = aic.source.block<4, 3>(0, 0) *
                 sourceSplines[4].coefficients;
    out.doublet = aic.doublet.block<4, 6>(0, 0) * qa.coefficients;
    out.valid = true;
    return out;
}

} // namespace a502
