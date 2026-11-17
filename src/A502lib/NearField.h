#pragma once

// NearField.h
// -----------
// Subsonic near-field subpanel influence coefficients (subsbi in panair.f90).
//
// The local subpanel coordinates have z along the subpanel normal and the
// field point is x.  The coefficients map the source and doublet Taylor
// series in the local subpanel frame to the potential and velocity at x:
//
//   source:  (a0 + ax*x' + ay*y' + ...)  -> du(0..3, 0..5)
//   doublet: (b0 + bx*x' + by*y' + ...)  -> dv(0..3, 0..9)
//
// Row 0 is the potential coefficient; rows 1..3 are the local x,y,z
// velocity coefficients.  The polynomial basis is
// [1, x, y, x^2/2, x*y, y^2/2] for nf=6 and adds
// [x^3/6, x^2*y/2, x*y^2/2, y^3/6] for nf=10.  The factor aj is the ratio of
// global to local area elements.  The kernel normalization is 1/(4*pi).

#define EIGEN_NO_DEBUG

#include <array>

#include <Eigen/Dense>

namespace a502 {

struct PanelGeometry;
struct Compressibility;

struct SubsonicSubpanelInfluence {
    Eigen::Matrix<double, 4, 6> source{Eigen::Matrix<double, 4, 6>::Zero()};
    Eigen::Matrix<double, 4, 10> doublet{Eigen::Matrix<double, 4, 10>::Zero()};
    bool valid{false};
};

// vertices: local coordinates of the subpanel corners (ns = 3 or 4).
// ics: 0 for a quadrilateral, 1..4 when edge ics is collapsed (triangle).
// its: 1 source only, 2 doublet only, 3 source and doublet.
// nf: 6 (quadratic) or 10 (cubic) doublet terms.
SubsonicSubpanelInfluence subsonicSubpanelInfluence(
    const std::array<Eigen::Vector3d, 4>& vertices, int ns, int ics, int its,
    const Eigen::Vector3d& fieldPoint, double aj, int nf = 6);

// Supersonic near-field subpanel influence (supsbi in panair.f90) for
// subinclined subpanels.  Same local-frame conventions and output layout as
// the subsonic kernel.
using SupersonicSubpanelInfluence = SubsonicSubpanelInfluence;
SupersonicSubpanelInfluence supersonicSubpanelInfluence(
    const std::array<Eigen::Vector3d, 4>& vertices, int ns, int ics, int its,
    const Eigen::Vector3d& fieldPoint, double aj, int nf = 6);

// Supersonic superinclined subpanel influence (supspi + aicsup).  sgnx is
// the sign of the freestream component along the subpanel normal (A502
// /norx/ sgnx).  The returned row layout matches nftpic's expectation:
// row 0 is the potential and rows 1..3 are the local velocity components
// after supspi's component rearrangement.
SupersonicSubpanelInfluence supersonicSuperinclinedSubpanelInfluence(
    const std::array<Eigen::Vector3d, 4>& vertices, int ns, int ics, int its,
    const Eigen::Vector3d& fieldPoint, double aj, int nf = 6,
    double sgnx = 1.0);

// Near-field influence of one panel on a field point, summed over its eight
// triangular subpanels.  Rows are potential and global x,y,z velocity; source
// columns are sigma1..sigma4,sigma9 and doublet columns are mu1..mu9.
struct PanelNearFieldInfluence {
    Eigen::Matrix<double, 4, 5> source{Eigen::Matrix<double, 4, 5>::Zero()};
    Eigen::Matrix<double, 4, 9> doublet{Eigen::Matrix<double, 4, 9>::Zero()};
    bool valid{false};
};

// Only the quadratic (nf=6) panel splines are currently available.
// The panel must already have compressibility-scaled subpanel frames
// (applyCompressibleFrames), as built by AeroModel::solve().
PanelNearFieldInfluence subsonicPanelNearFieldInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int nf = 6);

// nftpic() kernel dispatch: subsonic flows use subsbi(), supersonic flows
// use supsbi() for subinclined subpanels.  Supersonic superinclined
// subpanels need supspi(), which is not ported yet, so they return invalid.
PanelNearFieldInfluence panelNearFieldInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int nf = 6);

// Quasi-far (intermediate) influence of qffcal(): one quadrilateral/panel
// frame evaluation instead of the eight subpanel calls.  The result maps
// the panel Taylor coefficients to the field point.
struct PanelIntermediateInfluence {
    Eigen::Matrix<double, 4, 6> source{Eigen::Matrix<double, 4, 6>::Zero()};
    Eigen::Matrix<double, 4, 10> doublet{Eigen::Matrix<double, 4, 10>::Zero()};
    bool valid{false};
};

PanelIntermediateInfluence subsonicPanelIntermediateInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int nf = 6);

} // namespace a502
