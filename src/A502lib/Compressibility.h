#pragma once

// Compressibility.h
// -----------------
// Freestream / Prandtl-Glauert metric used by PAN AIR:
//
//   betams = 1 - M^2,  betam = sqrt(|betams|),  sbetam = sign(betams)
//   compd  = compressibility direction (unit vector)
//   ggcp   = scaled reference -> compressibility-axis transformation
//   ggcpit = its companion used by refloc()
//
// The formulas follow fsvcal() and refloc() in panair.f90.

#define EIGEN_NO_DEBUG

#include <Eigen/Dense>

namespace a502 {

struct Compressibility {
    double mach{0.0};
    double betams{1.0};   // 1 - M^2
    double abetms{1.0};   // |1 - M^2|
    double betam{1.0};    // sqrt(abetms)
    double sbetam{1.0};   // sign(1 - M^2)
    double btsqi{1.0};    // 1 / abetms

    Eigen::Vector3d compd{Eigen::Vector3d::UnitX()};
    Eigen::Matrix3d arotc{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d arotci{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d ggcp{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d ggcpit{Eigen::Matrix3d::Identity()};

    // alphaDeg = alpc, betaDeg = betc (degrees).  Throws std::invalid_argument
    // for transonic M = 1, where the Prandtl-Glauert transformation is
    // singular.
    static Compressibility fromAngles(double mach, double alphaDeg,
                                      double betaDeg);
};

// Local frame produced by refloc(): local = toLocal * (global - origin),
// global = origin + toGlobal * local.  The frame is built for the
// compressibility-scaled metric, so its z axis is the subpanel normal and
// in-plane points have local z = 0.  'inclination' is rpntyp (+1 subinclined,
// -1 superinclined).
struct ScaledFrame {
    Eigen::Matrix3d toLocal{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d toGlobal{Eigen::Matrix3d::Identity()};
    double          jacobian{1.0};
    int             inclination{1};
    bool            valid{false};
};

ScaledFrame buildReflocFrame(const Compressibility& flow,
                             const Eigen::Vector3d& normal,
                             const Eigen::Vector3d& reference);

struct PanelGeometry;
void applyCompressibleFrames(PanelGeometry& panel, const Compressibility& flow);

} // namespace a502
