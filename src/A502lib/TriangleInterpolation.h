#pragma once

// TriangleInterpolation.h
// -----------------------
// Port of bqbfun() and the linear/quadratic parts of tcof() from panair.f90.
//
// bqbfun() evaluates the nine biquadratic basis functions of the panel
// isoparametric coordinates (s,t) in [-1,1]:
//
//   P1(1,1)   P5(0,1)   P2(-1,1)
//   P8(1,0)   P9(0,0)   P6(-1,0)
//   P4(1,-1)  P7(0,-1)  P3(-1,-1)
//
// The quadratic triangle map uses the same coefficient ordering as the
// PAN AIR subpanel splines:
//
//   f(x,y) = c0 + cx*x + cy*y + cxx*x^2/2 + cxy*x*y + cyy*y^2/2.
//
// The edge lambdas are  lambda = 2*f(midpoint) - (f(end1)+f(end2))/2.

#define EIGEN_NO_DEBUG

#include <Eigen/Dense>

namespace a502 {

struct BiquadraticBasis {
    Eigen::Matrix<double, 9, 1> phi{Eigen::Matrix<double, 9, 1>::Zero()};
    Eigen::Matrix<double, 9, 1> phis{Eigen::Matrix<double, 9, 1>::Zero()};
    Eigen::Matrix<double, 9, 1> phit{Eigen::Matrix<double, 9, 1>::Zero()};
};

BiquadraticBasis bqbfun(double s, double t);

struct TriangleLinearMap {
    // rows: coefficients of [1, x, y] for each vertex basis function
    Eigen::Matrix3d b{Eigen::Matrix3d::Zero()};
    bool valid{false};
};

struct TriangleQuadraticMap {
    // rows: [c0, cx, cy, cxx, cxy, cyy]
    // columns: vertex values v1,v2,v3 then edge lambdas l1,l2,l3
    // (l1 = edge v2-v3, l2 = edge v3-v1, l3 = edge v1-v2)
    Eigen::Matrix<double, 6, 6> c{Eigen::Matrix<double, 6, 6>::Zero()};
    bool valid{false};
};

// p.col(k) is the local (x,y) position of triangle vertex k, k=0,1,2.
TriangleLinearMap buildLinearTriangleMap(const Eigen::Matrix3d& p);
TriangleQuadraticMap buildQuadraticTriangleMap(const Eigen::Matrix3d& p);

} // namespace a502
