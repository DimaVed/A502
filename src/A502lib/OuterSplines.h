#pragma once

// OuterSplines.h
// --------------
// Outer spline matrices that express the panel singularity values in terms
// of the global singularity parameters.
//
// First stage: the source analysis spline for interior panels.  The panel
// source distribution is linear in the inner refloc frame, and the three
// linear coefficients are obtained by a weighted least-squares fit to the
// 3x3 block of neighboring panel-center source parameters (sing/lsqsf in
// panair.f90).  The panel's own parameter has the large weight 1e4.
//
// The result is a 5x9 matrix whose rows are the source values at P1..P4,P9
// and whose columns correspond to the nine neighboring panels (di,dj from
// -1..1).  Edge/corner panels need the special neighbor selection of
// Appendix I.1.1.2 and are marked invalid for now.

#define EIGEN_NO_DEBUG

#include <array>

#include <Eigen/Dense>

namespace a502 {

struct Network; // defined in AeroModel.h

struct SourceOuterSpline {
    // Network-local panel indices i + (nm-1)*j, ordered di=-1..1,
    // dj=-1..1.  -1 means that the neighbor does not exist.
    std::array<int, 9> neighborPanel{{-1, -1, -1, -1, -1,
                                      -1, -1, -1, -1}};
    // rows: P1..P4,P9; columns: the nine neighbor source parameters
    Eigen::Matrix<double, 5, 9> coefficients{
        Eigen::Matrix<double, 5, 9>::Zero()};
    bool valid{false};
};

SourceOuterSpline buildSourceOuterSpline(const Network& network,
                                         int panelI, int panelJ);

} // namespace a502
