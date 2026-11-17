#pragma once

// FlowSolution.h
// --------------
// Surface velocity and pressure coefficient from a solved singularity
// vector.  The Cp rules follow cpcal() in panair.f90: linear, slender-body,
// second-order, and isentropic.  Material ratios default to air (mat 0).

#define EIGEN_NO_DEBUG

#include <vector>

#include <Eigen/Dense>

#include "GlobalBoundarySystem.h"

namespace a502 {

struct SurfaceFlowResult {
    std::vector<Eigen::Vector3d> velocity;
    std::vector<double> cp;
    bool valid{false};
};

enum class CpFormula {
    Linear,      // cp(1) in cpcal
    SlenderBody, // cp(2)
    SecondOrder, // cp(3)
    Isentropic   // cp(4)
};

// Material ratios from A502 matprp/matvfs.  Air (material 0) has all ratios
// equal to one.
struct CpMaterialRatios {
    double velocity{1.0};    // vfmat
    double temperature{1.0}; // tratio
    double pressure{1.0};    // pratio
    double cpFactor{1.0};    // cpfmat = vfmat*wfmat
};

struct MaterialProperties; // defined in AeroModel.h

// Port of matvfs(): material ratios for the current freestream.  Below
// Mach 0.01 the ratios of air are returned regardless of the material.
CpMaterialRatios materialRatios(const MaterialProperties& material,
                                const Compressibility& flow);

// cpcal(): perturbationVelocity is the perturbation velocity (total minus
// freestream), freestream is the freestream velocity vector, mach is the
// freestream Mach number.
double pressureCoefficient(
    const Eigen::Vector3d& perturbationVelocity,
    const Eigen::Vector3d& freestream, double mach, CpFormula formula,
    const CpMaterialRatios& material = {});

// Linear pressure coefficient for the total velocity:
//   cp = -2 (V - Vinf) . Vinf / |Vinf|^2.
double linearPressureCoefficient(const Eigen::Vector3d& velocity,
                                 const Eigen::Vector3d& freestream);

// Evaluates potential/velocity rows at all control points of all networks
// and applies the selected Cp rule.
SurfaceFlowResult evaluateSurfaceFlow(
    const std::vector<Network>& networks, const Compressibility& flow,
    const GlobalSingularityLayout& layout,
    const GlobalParameterValues& values, bool includeEdgeAndCorner = true,
    CpFormula formula = CpFormula::Linear,
    const CpMaterialRatios& material = {});

} // namespace a502
