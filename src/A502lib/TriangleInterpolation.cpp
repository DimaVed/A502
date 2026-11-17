// TriangleInterpolation.cpp
// -------------------------
// bqbfun() and tcof() (linear/quadratic parts) from panair.f90.

#include "TriangleInterpolation.h"

#include <cmath>

namespace a502 {

namespace {
constexpr double kEps = 1.0e-14;
}

BiquadraticBasis bqbfun(double s, double t) {
    BiquadraticBasis out;

    const double bs[3] = {0.5 * s * (s + 1.0), 1.0 - s * s,
                          0.5 * s * (s - 1.0)};
    const double bt[3] = {0.5 * t * (t + 1.0), 1.0 - t * t,
                          0.5 * t * (t - 1.0)};
    const double ds[3] = {s + 0.5, -2.0 * s, s - 0.5};
    const double dt[3] = {t + 0.5, -2.0 * t, t - 0.5};

    // phi(1..9) correspond to P1..P9; the Fortran ordering is retained:
    // phi(1)=bs1*bt1, phi(5)=bs2*bt1, phi(2)=bs3*bt1,
    // phi(8)=bs1*bt2, phi(9)=bs2*bt2, phi(6)=bs3*bt2,
    // phi(4)=bs1*bt3, phi(7)=bs2*bt3, phi(3)=bs3*bt3.
    out.phi(0) = bs[0] * bt[0]; // P1
    out.phi(4) = bs[1] * bt[0]; // P5
    out.phi(1) = bs[2] * bt[0]; // P2
    out.phi(7) = bs[0] * bt[1]; // P8
    out.phi(8) = bs[1] * bt[1]; // P9
    out.phi(5) = bs[2] * bt[1]; // P6
    out.phi(3) = bs[0] * bt[2]; // P4
    out.phi(6) = bs[1] * bt[2]; // P7
    out.phi(2) = bs[2] * bt[2]; // P3

    out.phis(0) = ds[0] * bt[0];
    out.phis(4) = ds[1] * bt[0];
    out.phis(1) = ds[2] * bt[0];
    out.phis(7) = ds[0] * bt[1];
    out.phis(8) = ds[1] * bt[1];
    out.phis(5) = ds[2] * bt[1];
    out.phis(3) = ds[0] * bt[2];
    out.phis(6) = ds[1] * bt[2];
    out.phis(2) = ds[2] * bt[2];

    out.phit(0) = bs[0] * dt[0];
    out.phit(4) = bs[1] * dt[0];
    out.phit(1) = bs[2] * dt[0];
    out.phit(7) = bs[0] * dt[1];
    out.phit(8) = bs[1] * dt[1];
    out.phit(5) = bs[2] * dt[1];
    out.phit(3) = bs[0] * dt[2];
    out.phit(6) = bs[1] * dt[2];
    out.phit(2) = bs[2] * dt[2];

    return out;
}

TriangleLinearMap buildLinearTriangleMap(const Eigen::Matrix3d& p) {
    TriangleLinearMap out;

    const double det = (p(0, 1) - p(0, 0)) * (p(1, 2) - p(1, 0)) -
                       (p(0, 2) - p(0, 0)) * (p(1, 1) - p(1, 0));
    if (std::abs(det) <= kEps)
        return out;

    for (int i = 0; i < 3; ++i) {
        const int ip1 = (i + 1) % 3;
        const int ip2 = (ip1 + 1) % 3;
        out.b(0, i) =
            (p(0, ip1) * p(1, ip2) - p(0, ip2) * p(1, ip1)) / det;
        out.b(1, i) = (p(1, ip1) - p(1, ip2)) / det;
        out.b(2, i) = (p(0, ip2) - p(0, ip1)) / det;
    }
    out.valid = true;
    return out;
}

TriangleQuadraticMap buildQuadraticTriangleMap(const Eigen::Matrix3d& p) {
    TriangleQuadraticMap out;
    const TriangleLinearMap lin = buildLinearTriangleMap(p);
    if (!lin.valid)
        return out;

    const Eigen::Matrix3d& b = lin.b;
    for (int i = 0; i < 3; ++i) {
        const int ip1 = (i + 1) % 3;
        const int ip2 = (ip1 + 1) % 3;

        out.c(0, i) = b(0, i) * b(0, i);
        out.c(1, i) = 2.0 * b(0, i) * b(1, i);
        out.c(2, i) = 2.0 * b(0, i) * b(2, i);
        out.c(3, i) = 2.0 * b(1, i) * b(1, i);
        out.c(4, i) = 2.0 * b(1, i) * b(2, i);
        out.c(5, i) = 2.0 * b(2, i) * b(2, i);

        out.c(0, i + 3) = 2.0 * b(0, ip1) * b(0, ip2);
        out.c(1, i + 3) =
            2.0 * (b(0, ip1) * b(1, ip2) + b(1, ip1) * b(0, ip2));
        out.c(2, i + 3) =
            2.0 * (b(0, ip1) * b(2, ip2) + b(2, ip1) * b(0, ip2));
        out.c(3, i + 3) = 4.0 * b(1, ip1) * b(1, ip2);
        out.c(4, i + 3) =
            2.0 * (b(1, ip1) * b(2, ip2) + b(2, ip1) * b(1, ip2));
        out.c(5, i + 3) = 4.0 * b(2, ip1) * b(2, ip2);
    }
    out.valid = true;
    return out;
}

} // namespace a502
