#pragma once

// Pic.h
// -----
// Influence-level dispatch of the panel influence coefficients, following
// the iflun branch of the AIC generator in panair.f90:
//
//   iflun = 1,2,3  ffpic()  far field (monopole/dipole/quadrupole)
//   iflun = 4      qffcal() one sub-panel intermediate field
//   iflun = 5      qnfcal() two sub-panel intermediate field
//   iflun = 6      qnfcal()/pifcal() eight sub-panel near field
//
// The two-subpanel qnfcal() branch is not ported yet; level 5 currently
// falls back to the eight-subpanel near field, which is the same continuous
// integral at a higher cost.

#include "NearField.h"

namespace a502 {

struct PanelGeometry;
struct Compressibility;

// qffcal() equivalent: one evaluation over the four panel corners in the
// inner-panel frame, mapped back to canonical source/doublet values.
PanelNearFieldInfluence panelIntermediateInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int its, int nf = 6);

// Selects the PIC method from influ().  level < 0 computes influenceLevel().
// its is the panel singularity type: 1 source, 2 doublet, 3 both.
PanelNearFieldInfluence panelPicInfluence(
    const PanelGeometry& panel, const Compressibility& flow,
    const Eigen::Vector3d& fieldPoint, int its, int nf = 6,
    int level = -1);

} // namespace a502
