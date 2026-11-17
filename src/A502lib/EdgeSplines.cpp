// EdgeSplines.cpp
// ---------------
// Port of edgls() from panair.f90.

#include "EdgeSplines.h"

#include <cmath>

namespace a502 {

EdgeLeastSquares edgls(const std::array<double, 4>& d,
                       const std::array<double, 4>& w) {
    EdgeLeastSquares out;

    if (d[1] == 0.0 && d[2] == 0.0) {
        out.a = {0.0, 0.5, 0.5, 0.0};
        return out;
    }

    const double x1 = -d[1] - 0.5 * d[0];
    const double x2 = -0.5 * d[1];
    const double x3 = 0.5 * d[2];
    const double x4 = d[2] + 0.5 * d[3];

    const double w1 = w[0];
    const double w4 = w[3];

    const double f =
        w1 * std::pow((x1 - x3) * (x1 - x2), 2) +
        w4 * std::pow((x4 - x3) * (x4 - x2), 2);

    if (!(std::abs(f) > 0.0)) {
        out.a = {0.0, 0.0, 0.0, 0.0};
        out.valid = false;
        return out;
    }

    out.a[0] = w1 * x2 * x3 * (x1 - x3) * (x1 - x2) / f;
    out.a[1] = x3 * (w1 * x1 * (x1 - x2) * std::pow(x1 - x3, 2) +
                     w4 * x4 * (x4 - x2) * std::pow(x4 - x3, 2)) /
               (f * (x3 - x2));
    out.a[2] = -x2 * (w1 * x1 * (x1 - x3) * std::pow(x1 - x2, 2) +
                      w4 * x4 * (x4 - x3) * std::pow(x4 - x2, 2)) /
               (f * (x3 - x2));
    out.a[3] = w4 * x2 * x3 * (x4 - x3) * (x4 - x2) / f;
    return out;
}

} // namespace a502
