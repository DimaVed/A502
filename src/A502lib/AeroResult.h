#pragma once

// AeroResult.h
// -----------
// Parser / ParaView export for a PanAir (A502) text solution file (*.out).
// Extracts per-panel flow quantities, off-body samples and configuration
// force/moment summaries from the marked report sections.

#define EIGEN_NO_DEBUG

#include <cstddef>
#include <string>
#include <vector>

#include <Eigen/Dense>

#include "AeroModel.h"

namespace a502 {

struct PanelResult {
    int             jc{0};
    int             ip{0};       // panel index within the network (1-based)
    int             network{0};  // 1-based network index from the report
    std::string     netname;
    Eigen::Vector3d xyz{Eigen::Vector3d::Zero()};
    Eigen::Vector3d w{Eigen::Vector3d::Zero()}; // velocity / mass-flux
    double          cp2nd{0.0};
    double          cpisn{0.0};
    double          mach{0.0};
    double          source{0.0};
    double          doublet{0.0};
};

struct OffBodyPoint {
    int             soln{0};
    int             pt{0};
    Eigen::Vector3d xyz{Eigen::Vector3d::Zero()};
    Eigen::Vector3d w{Eigen::Vector3d::Zero()};
    double          ppot{0.0};
    double          cp{0.0};
    double          mach{0.0};
};

struct ForceMomentSummary {
    int    soln{0};
    double alpha{0.0};
    double beta{0.0};
    double cl{0.0};
    double cdi{0.0};
    double cy{0.0};
    double fx{0.0};
    double fy{0.0};
    double fz{0.0};
    double mx{0.0};
    double my{0.0};
    double mz{0.0};
    double area{0.0};
    bool   fullConfig{false}; // true = full configuration (with symmetry)
};

class AeroResult {
public:
    void clear();
    void parse(const std::string& outPath);
    void printSummary() const;

    // Point cloud of panel CPs + off-body points in one VTU (region flag).
    void PrintParaview(const std::string& vtkPath) const;

    // Surface quads from mesh with Cp mapped by (netname, ip), j-major order.
    void PrintParaview(const std::string& vtkPath,
                       const AeroModel& mesh) const;

    std::vector<PanelResult>        panels;
    std::vector<OffBodyPoint>       offBody;
    std::vector<ForceMomentSummary> forceMoments;
    int                             nSolutions{0};
};

} // namespace a502
