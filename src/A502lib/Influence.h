#pragma once

// Influence.h
// -----------
// Panel/field-point influence classification (influ in panair.f90).
//
// PAN AIR chooses one of six PIC evaluation methods for every
// (field point, panel) pair:
//
//   0  no influence
//   1  monopole far field
//   2  dipole far field
//   3  quadrupole far field
//   4  one-subpanel intermediate field
//   5  two-subpanel intermediate field
//   6  eight-subpanel near field
//
// The classification is based on the compressibility-scaled panel diameter
// and the field-point distance.  The steady-flow case has zero phase, so the
// frequency upgrade in influ() is not applied here.

#define EIGEN_NO_DEBUG

#include <Eigen/Dense>

#include "Compressibility.h"

namespace a502 {

struct PanelGeometry;

// Compressible panel diameter from geomc():
//   diam = 2 * max_i sqrt((cp_i - cp_9, H (cp_i - cp_9)))
// where H scales the compressibility direction by |1-M^2|.
double compressedPanelDiameter(const PanelGeometry& panel,
                               const Compressibility& flow);

// Influence indicator of influ().  Returns 6 for a point at the panel centre
// (pcm <= 0 in the Fortran routine).
int influenceLevel(const PanelGeometry& panel, const Compressibility& flow,
                   const Eigen::Vector3d& fieldPoint);

} // namespace a502
