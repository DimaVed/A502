#pragma once

// SubpanelSplines.h
// -----------------
// Subpanel interpolation matrices for the PAN AIR higher-order singularity
// distributions.
//
// The source strength is linear on each triangular subpanel:
//
//   sigma(x',y') = a0 + ax*x' + ay*y'
//
// where (x',y') are the compressibility-scaled local coordinates used by
// refloc() and the origin is the first subpanel vertex (subpqr in
// panair.f90).  The five panel source values are sigma at P1..P4 and P9;
// values at the edge midpoints P5..P8 are the averages of the adjacent
// corner values (equation 5.5.5 of the Theory Document).
//
// This module builds the 3x5 matrix that maps [sigma1,sigma2,sigma3,sigma4,
// sigma9] to [a0,ax,ay] for each of the eight subpanels.  The doublet
// quadratic splines require the edge-lambda construction of subpqr/alam and
// will be added next.

#define EIGEN_NO_DEBUG

#include <array>

#include <Eigen/Dense>

#include "PanelGeometry.h"

namespace a502 {

using PanelSourceValues = Eigen::Matrix<double, 5, 1>;

struct SourceSubpanelSpline {
    // rows: a0, ax, ay; columns: sigma1..sigma4, sigma9
    Eigen::Matrix<double, 3, 5> coefficients{Eigen::Matrix<double, 3, 5>::Zero()};
    bool valid{false};
};

std::array<SourceSubpanelSpline, 8>
buildSourceSubpanelSplines(const PanelGeometry& panel);

// Returns [a0, ax, ay] for the given panel source values.
Eigen::Vector3d sourceTaylorCoefficients(
    const SourceSubpanelSpline& spline, const PanelSourceValues& values);

// Evaluates sigma at local coordinates (x',y').
double evaluateSource(const SourceSubpanelSpline& spline,
                      const Eigen::Vector2d& local,
                      const PanelSourceValues& values);

} // namespace a502
