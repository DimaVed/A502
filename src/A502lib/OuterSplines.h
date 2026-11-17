#pragma once

// OuterSplines.h
// --------------
// Outer spline matrices that express the panel singularity values in terms
// of the global singularity parameters.
//
// Source outer spline: the panel source distribution is linear in the inner
// refloc frame, and the three linear coefficients are obtained by a weighted
// least-squares fit to the 3x3 block of neighboring panel-center source
// parameters (sing/lsqsf in panair.f90).  The panel's own parameter has the
// large weight 1e4.  The result is a 5x9 matrix whose rows are the source
// values at P1..P4,P9 and whose columns correspond to the nine neighboring
// panels (di,dj from -1..1).
//
// Doublet outer spline: the doublet parameters live on the staggered
// (nm+1)x(nn+1) grid of gcpcal()/grdind() - panel centres in the interior,
// panel-edge midpoints on network edges, and network corners.  The sparse
// NetworkDoubletSpline covers all panels, including edge and corner panels,
// using blccal()'s one-dimensional edge fits and clamped quadratic stencils.
// The older 9x25 DoubletOuterSpline is retained for deep interior panels.

#define EIGEN_NO_DEBUG

#include <array>
#include <vector>

#include <Eigen/Dense>

#include "Compressibility.h"

namespace a502 {

struct Network; // defined in AeroModel.h
struct DoubletRow {
    std::vector<int> index;
    std::vector<double> coefficient;

    // Adds value to the entry for parameter id, combining duplicate ids.
    void add(int id, double value);
};

// Doublet singularity parameters of one network, in the staggered
// (nm+1) x (nn+1) grid used by gcpcal()/grdind().  Interior grid points are
// panel centres, edge interior points are panel-edge midpoints, and the four
// corner points are network corners.  The ids are network-local and unique.
struct NetworkDoubletParameters {
    int nPanelM{0};
    int nPanelN{0};
    int nCenters{0};
    int count{0};
    // corners: (m=1,n=1), (1,nn+1), (nm+1,nn+1), (nm+1,1)
    std::array<int, 4> cornerId{{-1, -1, -1, -1}};
    // edge midpoints: top row, right column, bottom row, left column
    std::array<int, 4> edgeMidBase{{-1, -1, -1, -1}};
    std::vector<Eigen::Vector3d> location;

    int centerId(int panelI, int panelJ) const {
        return panelI + nPanelM * panelJ;
    }
    // (m,n) are 1-based singularity-parameter grid indices.
    int parameterId(int m, int n) const;
};

NetworkDoubletParameters buildNetworkDoubletParameters(const Network& network);

// Sparse form of the 9xN doublet outer-spline matrix for one panel.  Each
// canonical panel point depends only on the nearby panel-centre and network-
// edge parameters, so the rows are stored as index/coefficient lists.
struct DoubletPanelSpline {
    std::array<DoubletRow, 9> rows;
    bool valid{false};
};

DoubletPanelSpline buildDoubletPanelSpline(
    const Network& network, const Compressibility& flow, int panelI, int panelJ,
    const NetworkDoubletParameters& parameters);

// Network-level assembly of the sparse B_D matrix.  Rows are stored panel by
// panel in A502 column-major panel order (panelI fastest), nine canonical
// points per panel.
struct NetworkDoubletSpline {
    int nPanelM{0};
    int nPanelN{0};
    int nParameters{0};
    std::vector<DoubletRow> rows;
    bool valid{false};
};

NetworkDoubletSpline buildNetworkDoubletSpline(
    const Network& network, const Compressibility& flow,
    const NetworkDoubletParameters& parameters);

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

struct DoubletOuterSpline {
    // Network-local panel-center indices i + (nm-1)*j for the 5x5 block.
    std::array<int, 25> neighborPanel{};
    // rows: P1..P9; columns: the 25 panel-center doublet parameters
    Eigen::Matrix<double, 9, 25> coefficients{
        Eigen::Matrix<double, 9, 25>::Zero()};
    bool valid{false};
};

// Dense 9x25 doublet outer spline for deep interior panels, where all four
// panel corners and all four edge midpoints have the full blccal stencils.
DoubletOuterSpline buildDoubletOuterSpline(const Network& network,
                                           const Compressibility& flow,
                                           int panelI, int panelJ);

} // namespace a502
