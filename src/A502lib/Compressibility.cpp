// Compressibility.cpp
// -------------------
// Port of fsvcal() metric setup and refloc() local-frame construction.

#include "Compressibility.h"
#include "PanelGeometry.h"

#include <cmath>
#include <stdexcept>

namespace a502 {

namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kTol = 1.0e-14;

Eigen::Matrix3d rotateMatrix(double alpha, double beta) {
    // rotate() in panair.f90, written out in row-major order.
    const double ca = std::cos(alpha);
    const double sa = std::sin(alpha);
    const double cb = std::cos(beta);
    const double sb = std::sin(beta);

    Eigen::Matrix3d r;
    r <<  ca * cb, -ca * sb,  sa,
          sb,       cb,       0.0,
         -sa * cb,  sa * sb,  ca;
    return r;
}

} // namespace

Compressibility Compressibility::fromAngles(double mach, double alphaDeg,
                                            double betaDeg) {
    if (!std::isfinite(mach) || !std::isfinite(alphaDeg) ||
        !std::isfinite(betaDeg))
        throw std::invalid_argument("Compressibility::fromAngles: non-finite input");

    Compressibility c;
    c.mach = mach;
    c.betams = 1.0 - mach * mach;
    c.abetms = std::abs(c.betams);
    if (c.abetms < kTol)
        throw std::invalid_argument(
            "Compressibility::fromAngles: transonic flow (M = 1) is not supported");

    c.betam = std::sqrt(c.abetms);
    c.sbetam = (c.betams >= 0.0) ? 1.0 : -1.0;
    c.btsqi = 1.0 / c.abetms;

    const double aarg = alphaDeg * kPi / 180.0;
    const double barg = betaDeg * kPi / 180.0;
    c.compd = Eigen::Vector3d(std::cos(aarg) * std::cos(barg),
                              -std::cos(aarg) * std::sin(barg),
                              std::sin(aarg));

    c.arotc = rotateMatrix(aarg, barg);
    c.arotci = c.arotc.transpose(); // orthogonal

    // ggcp = reference -> P-G scaled coordinates; ggcpit is its companion.
    // Row i is scaled by gfac = {1, betam, betam}.
    const double gfac[3] = {1.0, c.betam, c.betam};
    for (int i = 0; i < 3; ++i) {
        c.ggcp.row(i) = c.arotc.row(i) * gfac[i];
        c.ggcpit.row(i) = c.arotc.row(i) / gfac[i];
    }

    return c;
}

ScaledFrame buildReflocFrame(const Compressibility& flow,
                             const Eigen::Vector3d& normal,
                             const Eigen::Vector3d& reference) {
    ScaledFrame out;

    const Eigen::Vector3d en = normal.normalized();
    const Eigen::Vector3d enb = flow.ggcpit * en;

    const double afsq = enb.y() * enb.y() + enb.z() * enb.z();
    const double gmsq = flow.sbetam * enb.x() * enb.x() + afsq;
    if (!(std::abs(gmsq) > kTol) || !std::isfinite(gmsq)) {
        // Mach-inclined panel: refloc aborts in the original code.
        return out;
    }
    const int rpntyp = (gmsq > 0.0) ? 1 : -1;

    const double gm = std::sqrt(std::abs(gmsq));
    const double af = std::sqrt(afsq);

    double cosa = 0.0;
    double sina = 1.0;
    if (afsq > kTol) {
        cosa = enb.z() / af;
        sina = enb.y() / af;
    }

    Eigen::Matrix3d aa = Eigen::Matrix3d::Zero();
    aa(0, 0) = af / gm;
    aa(0, 1) = -flow.sbetam * enb.x() / gm * sina;
    aa(0, 2) = -flow.sbetam * enb.x() / gm * cosa;
    aa(1, 1) = rpntyp * cosa;
    aa(1, 2) = -rpntyp * sina;
    aa(2, 0) = enb.x() / gm;
    aa(2, 1) = af / gm * sina;
    aa(2, 2) = af / gm * cosa;

    Eigen::Matrix3d a = aa * flow.ggcp;

    // Inverse transformation, as assembled in refloc().
    const double rs = rpntyp * flow.sbetam;
    Eigen::Matrix3d aaa = aa;
    aaa(0, 0) = rpntyp * aaa(0, 0);
    aaa(2, 1) = rpntyp * aaa(2, 1);
    aaa(2, 2) = rpntyp * aaa(2, 2);
    aaa(0, 1) = rs * aaa(0, 1);
    aaa(0, 2) = rs * aaa(0, 2);
    aaa(2, 0) = rs * aaa(2, 0);

    Eigen::Matrix3d ai = (aaa * flow.ggcpit).transpose();

    if (flow.mach <= 1.0 && reference.norm() > kTol) {
        // Rotate the local system so that 'reference' is parallel to x'.
        const Eigen::Vector3d gp = a * reference;
        const double gpnorm = std::hypot(gp.x(), gp.y());
        if (gpnorm > kTol) {
            const double cosr = gp.x() / gpnorm;
            const double sinr = gp.y() / gpnorm;
            for (int j = 0; j < 3; ++j) {
                const double t = cosr * a(0, j) + sinr * a(1, j);
                a(1, j) = -sinr * a(0, j) + cosr * a(1, j);
                a(0, j) = t;

                const double w = cosr * ai(j, 0) + sinr * ai(j, 1);
                ai(j, 1) = -sinr * ai(j, 0) + cosr * ai(j, 1);
                ai(j, 0) = w;
            }
        }
    }

    if (!a.allFinite() || !ai.allFinite())
        return out;

    out.toLocal = a;
    out.toGlobal = ai;
    out.jacobian = flow.btsqi / gm;
    out.inclination = rpntyp;
    out.valid = true;
    return out;
}

void applyCompressibleFrames(PanelGeometry& panel,
                             const Compressibility& flow) {
    if (panel.degenerate)
        return;

    // Reference directions used by geomc(): genref 1,2,3.
    const Eigen::Vector3d gen1 = panel.cp[3] - panel.cp[0]; // P4-P1
    const Eigen::Vector3d gen2 = panel.cp[6] - panel.cp[4]; // P7-P5
    const Eigen::Vector3d gen3 = panel.cp[2] - panel.cp[1]; // P3-P2
    static const int kgen[5] = {1, 3, 3, 1, 2}; // 1-based, as in geomc

    for (int s = 0; s < 8; ++s) {
        SubpanelGeometry& sp = panel.subpanels[s];
        const int kg = (s < 4) ? kgen[s] : kgen[4];
        const Eigen::Vector3d& gen =
            (kg == 1) ? gen1 : (kg == 2 ? gen2 : gen3);

        const ScaledFrame f = buildReflocFrame(flow, sp.normal, gen);
        sp.scaledValid = f.valid;
        if (f.valid) {
            sp.scaledToLocal = f.toLocal;
            sp.scaledToGlobal = f.toGlobal;
            sp.jacobian = f.jacobian;
            sp.inclination = f.inclination;
        }
    }
}

} // namespace a502
