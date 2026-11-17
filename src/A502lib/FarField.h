#pragma once

// FarField.h
// ----------
// Far-field panel influence coefficients (ffpic / ffdqgv / ffdqg in
// panair.f90).  The far-field AICs are expressed in the central subpanel
// frame and then mapped back to the panel canonical values (sigma at
// P1..P4,P9 and mu at P1..P9) so they can be used with the same outer
// splines as the near-field PIC.
//
// The panel moment matrix cf is built with panelAreaMoments() (panmom), the
// h/hb moment matrices follow ffdqgv(), and the doublet Taylor coefficients
// are mapped back to the canonical P1..P9 values with qcof().  Only the
// quadratic doublet basis (nf = 6) is available; the cubic terms are not
// ported yet.

#define EIGEN_NO_DEBUG

#include <array>

#include <Eigen/Dense>

namespace a502 {

struct PanelGeometry;
struct Compressibility;

struct FarFieldPanelData {
    Eigen::Vector3d center{Eigen::Vector3d::Zero()};
    Eigen::Matrix3d toLocal{Eigen::Matrix3d::Identity()};
    Eigen::Matrix3d toGlobal{Eigen::Matrix3d::Identity()};
    double ajf{1.0};
    double sf{1.0}; // sbetam
    double rf{1.0}; // central subpanel inclination
    Eigen::Matrix<double, 6, 6> cf{
        Eigen::Matrix<double, 6, 6>::Zero()};
    Eigen::Matrix<double, 10, 6> hm{
        Eigen::Matrix<double, 10, 6>::Zero()};
    // hbm(ka,kb) for ka = 0..5 (basis function) and kb = 0..2 (moment).
    std::array<Eigen::Matrix<double, 2, 3>, 6> hbm{};
    int its{3}; // 1 source only, 2 doublet only, 3 both
    bool valid{false};
};

// its follows the panel singularity type.  The panel must have valid
// compressibility-scaled subpanel frames.
FarFieldPanelData buildFarFieldPanelData(
    const PanelGeometry& panel, const Compressibility& flow, int its);

struct FarFieldAic {
    // rows: potential, local velocity vx,vy,vz (transformed to global for
    // velocity rows as in ffpic); source columns: 3 or 6 far-field source
    // coefficients; doublet columns: nf coefficients.
    Eigen::Matrix<double, 4, 6> source{
        Eigen::Matrix<double, 4, 6>::Zero()};
    Eigen::Matrix<double, 4, 10> doublet{
        Eigen::Matrix<double, 4, 10>::Zero()};
    bool valid{false};
};

// iflu: 1 monopole, 2 dipole, 3 quadrupole.  nf is 6 or 10.
FarFieldAic evaluateFarFieldAic(
    const FarFieldPanelData& data, const Eigen::Vector3d& fieldPoint,
    int iflu, int nf);

struct PanelFarFieldInfluence {
    // Same canonical-value layout as PanelNearFieldInfluence.
    Eigen::Matrix<double, 4, 5> source{
        Eigen::Matrix<double, 4, 5>::Zero()};
    Eigen::Matrix<double, 4, 9> doublet{
        Eigen::Matrix<double, 4, 9>::Zero()};
    bool valid{false};
};

// Far-field panel influence mapped to the canonical panel singularity
// values.  Only nf = 6 is supported.
PanelFarFieldInfluence farFieldPanelInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int iflu, int nf, int its);

} // namespace a502
