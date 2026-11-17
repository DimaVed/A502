#pragma once

// EdgeSplines.h
// -------------
// Port of edgls(): value of a quadratic at an edge grid point from the four
// neighbouring edge parameters (two on each side), with the two adjacent
// values fitted exactly and the outer two in a weighted least-squares sense.
//
// The four parameters are at the midpoints of intervals whose lengths are
// d[0..3]; the evaluation point is the common endpoint of intervals 2 and 3:
//
//   x1 = -d1 - d0/2,  x2 = -d1/2,  x3 = d2/2,  x4 = d2 + d3/2.
//
// f(0) = a[0] f(x1) + a[1] f(x2) + a[2] f(x3) + a[3] f(x4).

#define EIGEN_NO_DEBUG

#include <array>

namespace a502 {

struct EdgeLeastSquares {
    std::array<double, 4> a{{0.0, 0.5, 0.5, 0.0}};
    bool valid{true};
};

EdgeLeastSquares edgls(const std::array<double, 4>& d,
                       const std::array<double, 4>& w);

} // namespace a502
