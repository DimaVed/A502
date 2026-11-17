// BoundaryConditions.cpp
// ----------------------
// Port of the nct=1/nct=4 core of aical().

#include "BoundaryConditions.h"

#include "AeroModel.h"
#include "OuterSplines.h"
#include "PanelGeometry.h"
#include "SubpanelSplines.h"

#include <cmath>

namespace a502 {

namespace {

Eigen::Vector3d compressibilityScale(double bs,
                                     const Eigen::Vector3d& direction,
                                     const Eigen::Vector3d& v) {
    const double c2 = direction.squaredNorm();
    if (c2 <= 1.0e-30)
        return v;
    const double f = (bs - 1.0) * v.dot(direction) / c2;
    return v + f * direction;
}

} // namespace

LocalSingularityTerms buildPanelPointLocalTerms(
    const Network& network, const Compressibility& flow, int panelI,
    int panelJ, const Eigen::Vector3d& fieldPoint,
    const NetworkDoubletParameters& doubletParameters,
    bool hasSource, bool hasDoublet) {
    LocalSingularityTerms out;
    const int nPanelM = network.nm - 1;
    const int nPanelN = network.nn - 1;
    if (nPanelM <= 0 || nPanelN <= 0 || panelI < 0 || panelJ < 0 ||
        panelI >= nPanelM || panelJ >= nPanelN)
        return out;
    if (!hasSource && !hasDoublet) {
        out.valid = true;
        return out;
    }

    PanelGeometry panel = buildPanelGeometry(
        network.coordinates, network.nm, network.nn, panelI, panelJ);
    if (panel.degenerate)
        return out;
    applyCompressibleFrames(panel, flow);

    const auto sourceSubpanels = buildSourceSubpanelSplines(panel);
    const auto doubletSubpanels = buildDoubletSubpanelSplines(panel, flow);

    auto insideTriangle = [](const Eigen::Vector2d& p,
                             const Eigen::Vector2d& a,
                             const Eigen::Vector2d& b,
                             const Eigen::Vector2d& c,
                             double& minWeight) {
        const double den =
            (b.x() - a.x()) * (c.y() - a.y()) -
            (c.x() - a.x()) * (b.y() - a.y());
        if (std::abs(den) <= 1.0e-30)
            return false;
        const double w1 =
            ((p.x() - a.x()) * (c.y() - a.y()) -
             (c.x() - a.x()) * (p.y() - a.y())) /
            den;
        const double w2 =
            ((b.x() - a.x()) * (p.y() - a.y()) -
             (p.x() - a.x()) * (b.y() - a.y())) /
            den;
        const double w0 = 1.0 - w1 - w2;
        minWeight = std::min({w0, w1, w2});
        return minWeight >= -1.0e-8;
    };

    double scale = 1.0;
    for (int k = 0; k < 9; ++k)
        scale = std::max(scale, (panel.cp[k] - panel.cp[8]).norm());

    int sub = -1;
    double bestInterior = -1.0e30;
    double bestZ = 1.0e30;
    for (int s = 0; s < 8; ++s) {
        const std::size_t si = static_cast<std::size_t>(s);
        if (panel.subpanels[si].degenerate)
            continue;
        if (hasSource && !sourceSubpanels[si].valid)
            continue;
        if (hasDoublet && !doubletSubpanels[si].valid)
            continue;
        const SubpanelGeometry& sp = panel.subpanels[si];
        const Eigen::Matrix3d& frame =
            sp.scaledValid ? sp.scaledToLocal : sp.toLocal;
        const Eigen::Vector3d p0 =
            frame * (sp.vertex[0] - sp.origin);
        const Eigen::Vector3d p1 =
            frame * (sp.vertex[1] - sp.origin);
        const Eigen::Vector3d p2 =
            frame * (sp.vertex[2] - sp.origin);
        const Eigen::Vector3d localField =
            frame * (fieldPoint - sp.origin);
        double minWeight = -1.0;
        if (!insideTriangle(Eigen::Vector2d(localField.x(), localField.y()),
                            Eigen::Vector2d(p0.x(), p0.y()),
                            Eigen::Vector2d(p1.x(), p1.y()),
                            Eigen::Vector2d(p2.x(), p2.y()),
                            minWeight))
            continue;
        const double zt = std::abs(localField.z());
        if (zt > 1.0e-6 * scale)
            continue;
        if (minWeight > bestInterior ||
            (minWeight == bestInterior && zt < bestZ)) {
            bestInterior = minWeight;
            bestZ = zt;
            sub = s;
        }
    }
    if (sub < 0)
        return out;

    const std::size_t si = static_cast<std::size_t>(sub);
    const SubpanelGeometry& sp = panel.subpanels[si];
    const Eigen::Matrix3d& frame =
        sp.scaledValid ? sp.scaledToLocal : sp.toLocal;

    if (hasSource) {
        // Source strength at the field point in terms of canonical values.
        const SourceSubpanelSpline& sourceSpline = sourceSubpanels[si];
        const Eigen::Vector3d sourceLocal =
            frame * (fieldPoint - sourceSpline.origin);
        // The source spline rows are [a0, ax, ay]; columns are sigma1..4,sigma9.
        Eigen::Matrix<double, 1, 5> valueBasis =
            (Eigen::Matrix<double, 1, 3>() << 1.0, sourceLocal.x(),
             sourceLocal.y()).finished() *
            sourceSpline.coefficients;

        const SourceOuterSpline sourceOuter =
            buildSourceOuterSpline(network, panelI, panelJ);
        if (!sourceOuter.valid)
            return LocalSingularityTerms{};
        const Eigen::Matrix<double, 1, 9> sourceRow =
            valueBasis * sourceOuter.coefficients;
        for (int k = 0; k < 9; ++k) {
            const int id =
                sourceOuter.neighborPanel[static_cast<std::size_t>(k)];
            if (id < 0)
                continue;
            out.sourceId.push_back(id);
            out.sourceValue.conservativeResize(out.sourceValue.size() + 1);
            out.sourceValue(out.sourceValue.size() - 1) = sourceRow(k);
        }
    }

    if (hasDoublet) {
        const DoubletPanelSpline doubletOuter = buildDoubletPanelSpline(
            network, flow, panelI, panelJ, doubletParameters);
        if (!doubletOuter.valid)
            return LocalSingularityTerms{};

        const DoubletSubpanelSpline& doubletSpline = doubletSubpanels[si];
        const Eigen::Vector3d doubletLocal =
            frame * (fieldPoint - doubletSpline.origin);
        const double x = doubletLocal.x();
        const double y = doubletLocal.y();
        Eigen::Matrix<double, 1, 6> valueBasis6;
        valueBasis6 << 1.0, x, y, 0.5 * x * x, x * y, 0.5 * y * y;
        Eigen::Matrix<double, 1, 9> valueCanonical =
            valueBasis6 * doubletSpline.coefficients;

        Eigen::Matrix<double, 1, 6> dx;
        dx << 0.0, 1.0, 0.0, x, y, 0.0;
        Eigen::Matrix<double, 1, 6> dy;
        dy << 0.0, 0.0, 1.0, 0.0, x, y;
        Eigen::Matrix<double, 2, 9> gradCanonicalLocal;
        gradCanonicalLocal.row(0) = dx * doubletSpline.coefficients;
        gradCanonicalLocal.row(1) = dy * doubletSpline.coefficients;
        const Eigen::Matrix3d toGlobal =
            sp.scaledValid ? sp.scaledToGlobal : sp.toGlobal;
        Eigen::Matrix<double, 3, 9> gradCanonical;
        for (int j = 0; j < 9; ++j) {
            gradCanonical.col(j) =
                toGlobal * Eigen::Vector3d(gradCanonicalLocal(0, j),
                                           gradCanonicalLocal(1, j), 0.0);
        }

        auto findOrAdd = [&](int id) -> int {
            for (std::size_t k = 0; k < out.doubletId.size(); ++k) {
                if (out.doubletId[k] == id)
                    return static_cast<int>(k);
            }
            out.doubletId.push_back(id);
            return static_cast<int>(out.doubletId.size() - 1);
        };
        for (int row = 0; row < 9; ++row) {
            for (int id :
                 doubletOuter.rows[static_cast<std::size_t>(row)].index)
                findOrAdd(id);
        }

        out.doubletValue = Eigen::Matrix<double, 4, Eigen::Dynamic>::Zero(
            4, static_cast<Eigen::Index>(out.doubletId.size()));
        for (int row = 0; row < 9; ++row) {
            const DoubletRow& splineRow =
                doubletOuter.rows[static_cast<std::size_t>(row)];
            for (std::size_t k = 0; k < splineRow.index.size(); ++k) {
                const int col = findOrAdd(splineRow.index[k]);
                const double c = splineRow.coefficient[k];
                out.doubletValue(0, col) +=
                    c * valueCanonical(row);
                out.doubletValue.col(col).segment<3>(1) +=
                    c * gradCanonical.col(row);
            }
        }
    }

    out.valid = true;
    return out;
}

LocalSingularityTerms buildPanelCenterLocalTerms(
    const Network& network, const Compressibility& flow, int panelI,
    int panelJ, const NetworkDoubletParameters& doubletParameters,
    bool hasSource, bool hasDoublet) {
    PanelGeometry panel = buildPanelGeometry(
        network.coordinates, network.nm, network.nn, panelI, panelJ);
    if (panel.degenerate)
        return LocalSingularityTerms{};
    return buildPanelPointLocalTerms(
        network, flow, panelI, panelJ, panel.cp[8], doubletParameters,
        hasSource, hasDoublet);
}

BoundaryRow buildBoundaryRow(
    const GlobalInfluence& influence, const Eigen::Vector3d& normal,
    const Compressibility& flow,
    const BoundaryConditionCoefficients& coefficients,
    const LocalSingularityTerms& local) {
    BoundaryRow out;
    if (!influence.valid || !local.valid)
        return out;
    if (coefficients.nct != 1 && coefficients.nct != 2 &&
        coefficients.nct != 4)
        return out;

    const double nnorm = normal.norm();
    if (!(nnorm > 0.0))
        return out;
    const Eigen::Vector3d n = normal / nnorm;

    const Eigen::Index nSource = influence.source.cols();
    const Eigen::Index nDoublet = influence.doublet.cols();
    out.source = Eigen::VectorXd::Zero(nSource);
    out.doublet = Eigen::VectorXd::Zero(nDoublet);
    out.rhs = coefficients.bet;

    if (coefficients.nct == 1) {
        const double ca = coefficients.cu + coefficients.cl;
        const double cd = 0.5 * (coefficients.cu - coefficients.cl);
        const Eigen::Vector3d tn =
            compressibilityScale(flow.betams, flow.compd, n);

        for (Eigen::Index j = 0; j < nSource; ++j) {
            const double v = tn.x() * influence.source(1, j) +
                             tn.y() * influence.source(2, j) +
                             tn.z() * influence.source(3, j);
            out.source(j) = ca * v;
        }
        for (Eigen::Index j = 0; j < nDoublet; ++j) {
            const double v = tn.x() * influence.doublet(1, j) +
                             tn.y() * influence.doublet(2, j) +
                             tn.z() * influence.doublet(3, j);
            out.doublet(j) = ca * v;
        }

        for (std::size_t k = 0; k < local.sourceId.size(); ++k) {
            const int id = local.sourceId[k];
            if (id < 0 || id >= nSource)
                return BoundaryRow{};
            out.source(id) +=
                cd * local.sourceValue(static_cast<Eigen::Index>(k));
        }
    } else if (coefficients.nct == 4) {
        const double da = coefficients.du + coefficients.dl;
        const double dd = 0.5 * (coefficients.du - coefficients.dl);

        out.source = da * influence.source.row(0).transpose();
        out.doublet = da * influence.doublet.row(0).transpose();

        for (std::size_t k = 0; k < local.doubletId.size(); ++k) {
            const int id = local.doubletId[k];
            if (id < 0 || id >= nDoublet)
                return BoundaryRow{};
            out.doublet(id) +=
                dd * local.doubletValue(0, static_cast<Eigen::Index>(k));
        }
    } else {
        // General (nct=2) coefficients: aical's ta/td projection.
        const double ca = coefficients.cu + coefficients.cl;
        const double cd = 0.5 * (coefficients.cu - coefficients.cl);
        const double da = coefficients.du + coefficients.dl;
        const double dd = 0.5 * (coefficients.du - coefficients.dl);
        const Eigen::Vector3d tn =
            compressibilityScale(flow.betams, flow.compd, n);

        auto projectTangential =
            [&](const Eigen::Vector3d& v) -> Eigen::Vector3d {
            const double v2 = v.squaredNorm();
            const double vn = v.dot(n);
            if (v2 > vn * vn + 1.0e-30)
                return std::sqrt(1.0 / (1.0 - vn * vn / v2)) *
                       (v - vn * n);
            return Eigen::Vector3d::Zero();
        };

        const Eigen::Vector3d ta =
            projectTangential(coefficients.tu + coefficients.tl) +
            ca * tn;
        const Eigen::Vector3d td =
            projectTangential(0.5 * (coefficients.tu - coefficients.tl)) +
            cd * tn;

        for (Eigen::Index j = 0; j < nSource; ++j) {
            out.source(j) =
                da * influence.source(0, j) +
                ta.x() * influence.source(1, j) +
                ta.y() * influence.source(2, j) +
                ta.z() * influence.source(3, j);
        }
        for (Eigen::Index j = 0; j < nDoublet; ++j) {
            out.doublet(j) =
                da * influence.doublet(0, j) +
                ta.x() * influence.doublet(1, j) +
                ta.y() * influence.doublet(2, j) +
                ta.z() * influence.doublet(3, j);
        }

        const double factor = tn.dot(n);
        if (std::abs(factor) <= 1.0e-14)
            return BoundaryRow{};
        const double sourceFactor = td.dot(n) / factor;
        for (std::size_t k = 0; k < local.sourceId.size(); ++k) {
            const int id = local.sourceId[k];
            if (id < 0 || id >= nSource)
                return BoundaryRow{};
            out.source(id) +=
                sourceFactor * local.sourceValue(static_cast<Eigen::Index>(k));
        }

        for (std::size_t k = 0; k < local.doubletId.size(); ++k) {
            const int id = local.doubletId[k];
            if (id < 0 || id >= nDoublet)
                return BoundaryRow{};
            const Eigen::Index kk = static_cast<Eigen::Index>(k);
            out.doublet(id) +=
                dd * local.doubletValue(0, kk) +
                td.dot(local.doubletValue.col(kk).segment<3>(1));
        }
    }

    out.valid = true;
    return out;
}

BoundaryConditionCoefficients standardNormalMassFlux(
    double cu, double cl, const Eigen::Vector3d& normal,
    const Eigen::Vector3d& freestream) {
    BoundaryConditionCoefficients out;
    out.nct = 1;
    out.cu = cu;
    out.cl = cl;
    const double nnorm = normal.norm();
    out.bet = (nnorm > 0.0)
                  ? -normal.dot(freestream) / nnorm
                  : 0.0;
    return out;
}

BoundaryRightHandSide buildBoundaryRightHandSide(
    int nropt, const Compressibility& flow,
    const Eigen::Vector3d& normal,
    const Eigen::Vector3d& controlPoint,
    const Eigen::Vector3d& reference) {
    BoundaryRightHandSide out;
    const int option = std::abs(nropt);
    switch (option) {
    case 0:
        out.supported = true;
        out.hasValue = false;
        return out;
    case 1: // right-hand side explicitly zero
        out.value = 0.0;
        break;
    case 2: // minus freestream normal component
        out.value = -normal.dot(flow.compd);
        break;
    case 3: { // scaled co-position vector dotted into freestream
        const Eigen::Vector3d offset = controlPoint - reference;
        const Eigen::Vector3d scaled = compressibilityScale(
            1.0 / flow.betams, flow.compd, offset);
        out.value = -scaled.dot(flow.compd);
        break;
    }
    case 5:
    case 6:
    case 7: // freestream normal component
        out.value = normal.dot(flow.compd);
        break;
    case 8: // co-position vector dotted into freestream
        out.value = -(controlPoint - reference).dot(flow.compd);
        break;
    default:
        // nropt=4 uses the Fortran-coded cbet(); nropt=9 keeps the input
        // value.  Neither is ported.
        return out;
    }
    out.supported = true;
    out.hasValue = true;
    return out;
}

LeftHandBoundaryOption buildLeftHandBoundaryOption(
    int nlopt, const Compressibility& flow,
    const Eigen::Vector3d& normal,
    double upperCpFactor, double lowerCpFactor) {
    LeftHandBoundaryOption out;
    if (nlopt == 0) {
        out.supported = true;
        out.hasCondition = false;
        return out;
    }
    if (nlopt < 0) {
        // Negative option values select upper/lower surfaces in the BC map;
        // the map is not parsed yet.
        return out;
    }
    out.hasCondition = true;

    BoundaryConditionCoefficients bc;
    switch (nlopt) {
    // ---- normal perturbation mass flux (nct = 1) ------------------------
    case 2: // upper surface
        bc.nct = 1;
        bc.cu = 1.0;
        break;
    case 3: // lower surface
        bc.nct = 1;
        bc.cl = 1.0;
        break;
    case 4: // average
        bc.nct = 1;
        bc.cu = 0.5;
        bc.cl = 0.5;
        break;
    case 5: // difference
        bc.nct = 1;
        bc.cu = 1.0;
        bc.cl = -1.0;
        break;

    // ---- perturbation potential (nct = 4) -------------------------------
    case 6: // upper surface
        bc.nct = 4;
        bc.du = 1.0;
        break;
    case 7: // lower surface
        bc.nct = 4;
        bc.dl = 1.0;
        break;
    case 8: // average
        bc.nct = 4;
        bc.du = 0.5;
        bc.dl = 0.5;
        break;
    case 9: // difference
        bc.nct = 4;
        bc.du = 1.0;
        bc.dl = -1.0;
        break;
    case 15:
    case 16:
    case 17: // delta(cp...)-style conditions use the average potential
        bc.nct = 4;
        bc.du = 0.5;
        bc.dl = 0.5;
        break;

    // ---- material-based pressure conditions (nct = 2) -------------------
    case 18:
    case 19:
    case 20: {
        const double n2 = normal.squaredNorm();
        const double znznu =
            n2 + (flow.betams - 1.0) * normal.dot(flow.compd) *
                     normal.dot(flow.compd);
        if (std::abs(znznu) <= 1.0e-30)
            return out;
        const Eigen::Vector3d zncs =
            compressibilityScale(flow.betams, flow.compd, normal);

        bc.nct = 2;
        for (int side = 0; side < 2; ++side) {
            const double sgn = (side == 0) ? 1.0 : -1.0;
            const double cpFactor =
                (side == 0) ? upperCpFactor : lowerCpFactor;
            // fsvm(1) = 1: the freestream vector is the unit direction.
            const Eigen::Vector3d tg =
                (-2.0 * sgn * cpFactor) * flow.compd;
            const double fac = normal.dot(tg) / znznu;
            const Eigen::Vector3d tangential = tg - fac * zncs;
            if (side == 0) {
                bc.cu = fac;
                bc.tu = tangential;
            } else {
                bc.cl = fac;
                bc.tl = tangential;
            }
        }
        break;
    }

    // ---- normal perturbation velocity (nct = 2) -------------------------
    case 11:
    case 12:
    case 13:
    case 14: {
        const double nnorm = normal.norm();
        if (!(nnorm > 0.0))
            return out;
        const Eigen::Vector3d zn = normal / nnorm;
        const Eigen::Vector3d zncs =
            compressibilityScale(flow.betams, flow.compd, zn);

        // bcopt() constants al1, al2, enx1, enx2 default to 0,0,1,1; the
        // $compressibility axis-angle input options are not parsed yet.
        const double al1 = 0.0;
        const double al2 = 0.0;
        const double enx1 = 1.0;
        const double enx2 = 1.0;
        const double enx = flow.compd.dot(zn);
        const double enxm = flow.mach * enx;
        double alam = al1;
        if (enxm >= enx2) {
            alam = al2;
        } else if (enxm > enx1 && enx2 > enx1) {
            alam = (al2 * (enxm - enx1) + al1 * (enx2 - enxm)) /
                   (enx2 - enx1);
        }
        alam *= flow.mach * flow.mach * enx;

        const Eigen::Vector3d h = zn - alam * flow.compd;
        const double hznc = h.dot(zn);
        const double zncsn = zncs.dot(zn);
        if (std::abs(zncsn) <= 1.0e-30)
            return out;
        const double factor = hznc / zncsn;
        const Eigen::Vector3d tangential = h - factor * zncs;

        bc.nct = 2;
        if (nlopt == 11) { // upper surface
            bc.cu = factor;
            bc.tu = tangential;
        } else if (nlopt == 12) { // lower surface
            bc.cl = factor;
            bc.tl = tangential;
        } else if (nlopt == 13) { // average
            bc.cu = 0.5 * factor;
            bc.cl = 0.5 * factor;
            bc.tu = 0.5 * tangential;
            bc.tl = 0.5 * tangential;
        } else { // 14: difference
            bc.cu = factor;
            bc.cl = -factor;
            bc.tu = tangential;
            bc.tl = -tangential;
        }
        break;
    }

    default:
        // nlopt=1 (input coefficients), 10 (ccof), and 21..25 are not
        // ported yet.
        return out;
    }

    out.coefficients = bc;
    out.supported = true;
    return out;
}

BoundarySystemRows buildNetworkNormalFluxRows(
    const Network& network, const Compressibility& flow,
    const NetworkSingularityLayout& layout,
    const std::vector<ControlPoint>& controlPoints, double cu, double cl) {
    BoundarySystemRows out;

    std::vector<BoundaryRow> rows;
    std::vector<int> indices;
    rows.reserve(controlPoints.size());
    indices.reserve(controlPoints.size());

    for (std::size_t k = 0; k < controlPoints.size(); ++k) {
        const ControlPoint& cp = controlPoints[k];
        if (cp.kind != ControlPointKind::PanelCenter)
            continue;

        const NetworkInfluence raw = evaluateNetworkNearFieldInfluence(
            network, flow, cp.position, layout);
        if (!raw.valid)
            return BoundarySystemRows{};

        GlobalInfluence global;
        global.source = raw.source;
        global.doublet = raw.doublet;
        global.valid = true;

        const LocalSingularityTerms local = buildPanelPointLocalTerms(
            network, flow, cp.panelI, cp.panelJ, cp.position,
            layout.doublet, layout.hasSource, layout.hasDoublet);
        if (!local.valid)
            return BoundarySystemRows{};

        const BoundaryConditionCoefficients bc =
            standardNormalMassFlux(cu, cl, cp.normal, flow.compd);
        BoundaryRow row = buildBoundaryRow(
            global, cp.normal, flow, bc, local);
        if (!row.valid)
            return BoundarySystemRows{};

        rows.push_back(std::move(row));
        indices.push_back(static_cast<int>(k));
    }

    if (rows.empty())
        return out;

    out.nRows = static_cast<int>(rows.size());
    out.nSource = layout.source.count;
    out.nDoublet = layout.doublet.count;
    out.source = Eigen::MatrixXd::Zero(out.nRows, out.nSource);
    out.doublet = Eigen::MatrixXd::Zero(out.nRows, out.nDoublet);
    out.rhs = Eigen::VectorXd::Zero(out.nRows);
    out.controlPointIndex = std::move(indices);

    for (int r = 0; r < out.nRows; ++r) {
        const BoundaryRow& row = rows[static_cast<std::size_t>(r)];
        if (row.source.size() != out.nSource ||
            row.doublet.size() != out.nDoublet)
            return BoundarySystemRows{};
        out.source.row(r) = row.source.transpose();
        out.doublet.row(r) = row.doublet.transpose();
        out.rhs(r) = row.rhs;
    }

    out.valid = true;
    return out;
}

} // namespace a502
