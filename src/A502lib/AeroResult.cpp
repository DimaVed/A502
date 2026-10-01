// AeroResult.cpp
// -------------
// Implementation of the A502 / PanAir *.out solution parser and VTK export.

#include "AeroResult.h"
#include "Logger.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>

#ifdef WITH_VTK
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellTypes.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkNew.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkStringArray.h>
#include <vtkUnstructuredGrid.h>
#include <vtkXMLUnstructuredGridWriter.h>
#endif

namespace a502 {

namespace {

std::string trim(const std::string& s) {
    const char* ws = " \t\r";
    std::size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos)
        return std::string();
    std::size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

bool startsWith(const std::string& s, const char* prefix) {
    const std::size_t n = std::char_traits<char>::length(prefix);
    return s.size() >= n && s.compare(0, n, prefix) == 0;
}

bool containsCI(const std::string& hay, const char* needle) {
    auto lower = [](unsigned char c) {
        return static_cast<char>(std::tolower(c));
    };
    std::string h;
    h.reserve(hay.size());
    for (unsigned char c : hay)
        h.push_back(lower(c));
    std::string n;
    for (const char* p = needle; *p; ++p)
        n.push_back(lower(static_cast<unsigned char>(*p)));
    return h.find(n) != std::string::npos;
}

std::vector<std::string> tokenize(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok)
        out.push_back(tok);
    return out;
}

bool parseDouble(const std::string& tok, double& v) {
    if (tok.empty())
        return false;
    std::string t = tok;
    char last = t.back();
    if (last == 'D' || last == 'd')
        t.back() = 'E';
    try {
        std::size_t pos = 0;
        v = std::stod(t, &pos);
        return pos == t.size();
    } catch (...) {
        return false;
    }
}

bool parseInt(const std::string& tok, int& v) {
    double d = 0.0;
    if (!parseDouble(tok, d))
        return false;
    v = static_cast<int>(std::lround(d));
    return true;
}

// "network id:winga       index:   1    ..."
bool parseNetworkHeader(const std::string& line, std::string& name, int& index) {
    auto p = line.find("network id:");
    if (p == std::string::npos)
        return false;
    p += 11;
    while (p < line.size() && std::isspace(static_cast<unsigned char>(line[p])))
        ++p;
    std::size_t q = p;
    while (q < line.size() && !std::isspace(static_cast<unsigned char>(line[q])))
        ++q;
    name = trim(line.substr(p, q - p));
    auto ip = line.find("index:");
    if (ip == std::string::npos)
        return !name.empty();
    ip += 6;
    try {
        index = std::stoi(trim(line.substr(ip)));
    } catch (...) {
        index = 0;
    }
    return !name.empty();
}

} // anonymous namespace

void AeroResult::clear() {
    panels.clear();
    offBody.clear();
    forceMoments.clear();
    nSolutions = 0;
}

void AeroResult::parse(const std::string& outPath) {
    clear();
    BLOG(info) << "Parsing A502 result file: " << outPath;

    std::ifstream in(outPath);
    if (!in.is_open()) {
        BLOG(error) << "Cannot open result file \"" << outPath << "\"";
        throw std::runtime_error("AeroResult::parse: cannot open \"" + outPath +
                                 "\"");
    }

    enum class Mode { Scan, PanelTable, OffBody, ForceFull };
    Mode mode = Mode::Scan;
    std::string curNet;
    int curNetIndex = 0;
    bool awaitPanelHeader = false;
    bool extendedPanelFormat = false; // nacelle-style multi-line block
    int extPhase = 0;                 // 0=geom line, 1=upper Cp line, 2..3=skip
    PanelResult pendingPanel;
    bool forceAwaitData = false;
    bool forceIsFull = false;

    std::string line;
    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.empty())
            continue;

        if (containsCI(t, "simultaneous solution number")) {
            ++nSolutions;
            continue;
        }

        if (startsWith(t, "0*b*off-body") || startsWith(t, "*b*off-body")) {
            mode = Mode::OffBody;
            continue;
        }
        if (startsWith(t, "0*e*off-body") || startsWith(t, "*e*off-body")) {
            mode = Mode::Scan;
            continue;
        }
        if (startsWith(t, "0*b*streaml") || startsWith(t, "*b*streaml") ||
            startsWith(t, "0*b*section") || startsWith(t, "*b*section") ||
            startsWith(t, "0*b*for-mom") || startsWith(t, "*b*for-mom")) {
            if (mode == Mode::PanelTable)
                mode = Mode::Scan;
            awaitPanelHeader = false;
            extPhase = 0;
            continue;
        }

        if (containsCI(t, "full configuration forces and moments summary")) {
            mode = Mode::ForceFull;
            forceIsFull = true;
            forceAwaitData = false;
            continue;
        }
        if (containsCI(t, "input configuration forces and moments summary")) {
            mode = Mode::ForceFull;
            forceIsFull = false;
            forceAwaitData = false;
            continue;
        }

        std::string netName;
        int netIdx = 0;
        if (parseNetworkHeader(t, netName, netIdx)) {
            curNet = netName;
            curNetIndex = netIdx;
            mode = Mode::PanelTable;
            awaitPanelHeader = true;
            extendedPanelFormat = false;
            extPhase = 0;
            continue;
        }

        if (mode == Mode::PanelTable) {
            if (awaitPanelHeader) {
                // Short: "... wx wy wz cp2ndu ..."
                // Extended (nac6): separate rows with wxu / cplinu / cp2ndu
                if (containsCI(t, "wxu") || containsCI(t, "cplinu")) {
                    extendedPanelFormat = true;
                    awaitPanelHeader = false;
                } else if (containsCI(t, "cp2ndu") || containsCI(t, "cpisnu")) {
                    extendedPanelFormat = false;
                    awaitPanelHeader = false;
                }
                continue;
            }

            auto toks = tokenize(t);
            if (extendedPanelFormat) {
                if (extPhase == 0) {
                    // jc ip x y z ...
                    if (toks.size() < 5)
                        continue;
                    pendingPanel = PanelResult{};
                    if (!parseInt(toks[0], pendingPanel.jc) ||
                        !parseInt(toks[1], pendingPanel.ip))
                        continue;
                    double x = 0, y = 0, z = 0;
                    if (!parseDouble(toks[2], x) || !parseDouble(toks[3], y) ||
                        !parseDouble(toks[4], z))
                        continue;
                    pendingPanel.network = curNetIndex;
                    pendingPanel.netname = curNet;
                    pendingPanel.xyz = Eigen::Vector3d(x, y, z);
                    extPhase = 1;
                } else if (extPhase == 1) {
                    // lmachu wxu wyu wzu pheu vxu vyu vzu cplinu cpslnu cp2ndu cpisnu
                    if (toks.size() < 12) {
                        extPhase = 0;
                        continue;
                    }
                    double vals[12];
                    bool ok = true;
                    for (int i = 0; i < 12; ++i) {
                        if (!parseDouble(toks[static_cast<std::size_t>(i)],
                                         vals[i])) {
                            ok = false;
                            break;
                        }
                    }
                    if (!ok) {
                        extPhase = 0;
                        continue;
                    }
                    pendingPanel.mach = vals[0];
                    pendingPanel.w = Eigen::Vector3d(vals[1], vals[2], vals[3]);
                    pendingPanel.cp2nd = vals[10];
                    pendingPanel.cpisn = vals[11];
                    panels.push_back(pendingPanel);
                    extPhase = 2; // skip lower + delta rows
                } else {
                    ++extPhase;
                    if (extPhase > 3)
                        extPhase = 0;
                }
                continue;
            }

            // Short data rows:
            // jc ip x y z wx wy wz cp2ndu cpisnu lmachu source doublet
            if (toks.size() < 13)
                continue;
            PanelResult p;
            if (!parseInt(toks[0], p.jc) || !parseInt(toks[1], p.ip))
                continue;
            double vals[11];
            bool ok = true;
            for (int i = 0; i < 11; ++i) {
                if (!parseDouble(toks[static_cast<std::size_t>(i + 2)], vals[i])) {
                    ok = false;
                    break;
                }
            }
            if (!ok)
                continue;
            p.network = curNetIndex;
            p.netname = curNet;
            p.xyz = Eigen::Vector3d(vals[0], vals[1], vals[2]);
            p.w = Eigen::Vector3d(vals[3], vals[4], vals[5]);
            p.cp2nd = vals[6];
            p.cpisn = vals[7];
            p.mach = vals[8];
            p.source = vals[9];
            p.doublet = vals[10];
            panels.push_back(std::move(p));
            continue;
        }

        if (mode == Mode::OffBody) {
            if (containsCI(t, "soln") && containsCI(t, "pt") &&
                containsCI(t, "ppot"))
                continue;
            auto toks = tokenize(t);
            // soln pt x y z wx wy wz ppot cp mach  => 11 tokens
            if (toks.size() < 11)
                continue;
            OffBodyPoint o;
            if (!parseInt(toks[0], o.soln) || !parseInt(toks[1], o.pt))
                continue;
            double vals[9];
            bool ok = true;
            for (int i = 0; i < 9; ++i) {
                if (!parseDouble(toks[static_cast<std::size_t>(i + 2)], vals[i])) {
                    ok = false;
                    break;
                }
            }
            if (!ok)
                continue;
            o.xyz = Eigen::Vector3d(vals[0], vals[1], vals[2]);
            o.w = Eigen::Vector3d(vals[3], vals[4], vals[5]);
            o.ppot = vals[6];
            o.cp = vals[7];
            o.mach = vals[8];
            offBody.push_back(std::move(o));
            continue;
        }

        if (mode == Mode::ForceFull) {
            if (containsCI(t, "sol-no") && containsCI(t, "alpha")) {
                forceAwaitData = true;
                continue;
            }
            if (startsWith(t, "------"))
                continue;
            if (!forceAwaitData)
                continue;
            auto toks = tokenize(t);
            // First data line: sol alpha beta cl cdi cy fx fy fz
            // Second: mx my mz area
            if (toks.size() >= 9) {
                ForceMomentSummary fm;
                fm.fullConfig = forceIsFull;
                if (!parseInt(toks[0], fm.soln)) {
                    mode = Mode::Scan;
                    continue;
                }
                parseDouble(toks[1], fm.alpha);
                parseDouble(toks[2], fm.beta);
                parseDouble(toks[3], fm.cl);
                parseDouble(toks[4], fm.cdi);
                parseDouble(toks[5], fm.cy);
                parseDouble(toks[6], fm.fx);
                parseDouble(toks[7], fm.fy);
                parseDouble(toks[8], fm.fz);

                // Look ahead for the moment/area line
                std::string next;
                std::streampos pos = in.tellg();
                while (std::getline(in, next)) {
                    std::string nt = trim(next);
                    if (nt.empty() || startsWith(nt, "------") ||
                        startsWith(nt, "***"))
                        continue;
                    auto ntoks = tokenize(nt);
                    if (ntoks.size() >= 4) {
                        parseDouble(ntoks[0], fm.mx);
                        parseDouble(ntoks[1], fm.my);
                        parseDouble(ntoks[2], fm.mz);
                        parseDouble(ntoks[3], fm.area);
                    }
                    break;
                }
                (void)pos;
                forceMoments.push_back(std::move(fm));
                forceAwaitData = false;
                mode = Mode::Scan;
            }
            continue;
        }
    }

    BLOG(info) << "Parsed solutions=" << nSolutions
               << " panels=" << panels.size()
               << " off-body=" << offBody.size()
               << " F&M summaries=" << forceMoments.size();
}

void AeroResult::printSummary() const {
    BLOG(info) << "solutions: " << nSolutions;
    BLOG(info) << "panel results: " << panels.size();
    BLOG(info) << "off-body points: " << offBody.size();

    std::unordered_map<std::string, int> byNet;
    for (const auto& p : panels)
        ++byNet[p.netname.empty() ? std::string("(unnamed)") : p.netname];
    for (const auto& kv : byNet)
        BLOG(info) << "  network '" << kv.first << "' panels=" << kv.second;

    for (const auto& fm : forceMoments) {
        BLOG(info) << (fm.fullConfig ? "full" : "input")
                   << " config F&M sol=" << fm.soln
                   << " alpha=" << fm.alpha << " beta=" << fm.beta
                   << " CL=" << fm.cl << " CDi=" << fm.cdi << " CY=" << fm.cy
                   << " area=" << fm.area;
    }
}

#ifdef WITH_VTK
namespace {

void addScalar(vtkCellData* cd, const char* name, const std::vector<double>& v) {
    vtkNew<vtkDoubleArray> a;
    a->SetName(name);
    a->SetNumberOfComponents(1);
    a->SetNumberOfTuples(static_cast<vtkIdType>(v.size()));
    for (std::size_t i = 0; i < v.size(); ++i)
        a->SetValue(static_cast<vtkIdType>(i), v[i]);
    cd->AddArray(a);
}

void addInt(vtkCellData* cd, const char* name, const std::vector<int>& v) {
    vtkNew<vtkIntArray> a;
    a->SetName(name);
    a->SetNumberOfComponents(1);
    a->SetNumberOfTuples(static_cast<vtkIdType>(v.size()));
    for (std::size_t i = 0; i < v.size(); ++i)
        a->SetValue(static_cast<vtkIdType>(i), v[i]);
    cd->AddArray(a);
}

} // namespace
#endif

void AeroResult::PrintParaview(const std::string& vtkPath) const {
#ifdef WITH_VTK
    BLOG(info) << "Writing result VTU (point cloud): " << vtkPath;
    if (panels.empty() && offBody.empty()) {
        BLOG(error) << "PrintParaview: no panel or off-body data";
        throw std::runtime_error(
            "AeroResult::PrintParaview: no panel or off-body data");
    }

    vtkNew<vtkPoints> points;
    vtkNew<vtkUnstructuredGrid> grid;

    const std::size_t n = panels.size() + offBody.size();
    std::vector<double> cp2nd(n), cpisn(n), mach(n), source(n), doublet(n),
        ppot(n), region(n);
    std::vector<int> networkId(n), panelIp(n);
    std::vector<double> wx(n), wy(n), wz(n);
    vtkNew<vtkStringArray> netname;
    netname->SetName("netname");

    vtkIdType idx = 0;
    for (const auto& p : panels) {
        points->InsertNextPoint(p.xyz.x(), p.xyz.y(), p.xyz.z());
        vtkIdType pid = idx;
        grid->InsertNextCell(VTK_VERTEX, 1, &pid);
        const std::size_t i = static_cast<std::size_t>(idx);
        cp2nd[i] = p.cp2nd;
        cpisn[i] = p.cpisn;
        mach[i] = p.mach;
        source[i] = p.source;
        doublet[i] = p.doublet;
        ppot[i] = 0.0;
        region[i] = 0.0;
        networkId[i] = p.network;
        panelIp[i] = p.ip;
        wx[i] = p.w.x();
        wy[i] = p.w.y();
        wz[i] = p.w.z();
        netname->InsertNextValue(p.netname.c_str());
        ++idx;
    }
    for (const auto& o : offBody) {
        points->InsertNextPoint(o.xyz.x(), o.xyz.y(), o.xyz.z());
        vtkIdType pid = idx;
        grid->InsertNextCell(VTK_VERTEX, 1, &pid);
        const std::size_t i = static_cast<std::size_t>(idx);
        cp2nd[i] = o.cp;
        cpisn[i] = o.cp;
        mach[i] = o.mach;
        source[i] = 0.0;
        doublet[i] = 0.0;
        ppot[i] = o.ppot;
        region[i] = 1.0;
        networkId[i] = 0;
        panelIp[i] = o.pt;
        wx[i] = o.w.x();
        wy[i] = o.w.y();
        wz[i] = o.w.z();
        netname->InsertNextValue("offbody");
        ++idx;
    }

    grid->SetPoints(points);
    vtkCellData* cd = grid->GetCellData();
    addScalar(cd, "Cp_2nd", cp2nd);
    addScalar(cd, "Cp_isen", cpisn);
    addScalar(cd, "Mach", mach);
    addScalar(cd, "source", source);
    addScalar(cd, "doublet", doublet);
    addScalar(cd, "ppot", ppot);
    addScalar(cd, "region", region);
    addInt(cd, "network_id", networkId);
    addInt(cd, "panel_ip", panelIp);
    addScalar(cd, "Wx", wx);
    addScalar(cd, "Wy", wy);
    addScalar(cd, "Wz", wz);
    cd->AddArray(netname);

    vtkNew<vtkDoubleArray> wvec;
    wvec->SetName("W");
    wvec->SetNumberOfComponents(3);
    wvec->SetNumberOfTuples(static_cast<vtkIdType>(n));
    for (std::size_t i = 0; i < n; ++i)
        wvec->SetTuple3(static_cast<vtkIdType>(i), wx[i], wy[i], wz[i]);
    cd->AddArray(wvec);

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(vtkPath.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToAscii();
    if (writer->Write() == 0) {
        BLOG(error) << "VTK writer failed for \"" << vtkPath << "\"";
        throw std::runtime_error(
            "AeroResult::PrintParaview: VTK writer failed for \"" + vtkPath +
            "\"");
    }
    BLOG(info) << "Result VTU written: " << vtkPath << " (" << n << " cells)";
#else
    (void)vtkPath;
    BLOG(error) << "PrintParaview: built without VTK support (WITH_VTK)";
    throw std::runtime_error(
        "AeroResult::PrintParaview: built without VTK support (WITH_VTK)");
#endif
}

void AeroResult::PrintParaview(const std::string& vtkPath,
                               const AeroModel& mesh) const {
#ifdef WITH_VTK
    BLOG(info) << "Writing result VTU (mesh + Cp): " << vtkPath;
    if (mesh.networks.empty()) {
        BLOG(warning) << "Mesh has no networks; falling back to point cloud";
        PrintParaview(vtkPath);
        return;
    }

    // Panels for each network, sorted by global ip (report order ≈ j-major).
    std::unordered_map<std::string, std::vector<const PanelResult*>> byNet;
    for (const auto& p : panels)
        byNet[p.netname].push_back(&p);
    for (auto& kv : byNet) {
        std::sort(kv.second.begin(), kv.second.end(),
                  [](const PanelResult* a, const PanelResult* b) {
                      return a->ip < b->ip;
                  });
    }

    vtkNew<vtkPoints> points;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkStringArray> cellName;
    vtkNew<vtkIntArray> cellNet;
    vtkNew<vtkDoubleArray> cp2nd;
    vtkNew<vtkDoubleArray> cpisn;
    vtkNew<vtkDoubleArray> mach;
    vtkNew<vtkDoubleArray> source;
    vtkNew<vtkDoubleArray> doublet;
    vtkNew<vtkDoubleArray> wvec;
    cellName->SetName("netname");
    cellNet->SetName("network_id");
    cp2nd->SetName("Cp_2nd");
    cpisn->SetName("Cp_isen");
    mach->SetName("Mach");
    source->SetName("source");
    doublet->SetName("doublet");
    wvec->SetName("W");
    wvec->SetNumberOfComponents(3);
    grid->SetPoints(points);

    const double nan = std::numeric_limits<double>::quiet_NaN();
    int baseIndex = 0;
    std::size_t mapped = 0;
    std::size_t missing = 0;

    for (std::size_t netIdx = 0; netIdx < mesh.networks.size(); ++netIdx) {
        const Network& n = mesh.networks[netIdx];
        const int nm = n.nm > 0 ? n.nm : 0;
        const int nn = n.nn > 0 ? n.nn : 0;
        if (nm < 2 || nn < 2)
            continue;

        std::size_t nAvail = n.coordinates.size() / 3;
        std::size_t nWanted = static_cast<std::size_t>(nm) * nn;
        std::size_t nGridPts = std::min(nWanted, nAvail);
        if (nGridPts < 1)
            continue;

        std::vector<vtkIdType> localIds(static_cast<std::size_t>(nm) * nn, -1);
        for (std::size_t g = 0; g < nGridPts; ++g) {
            points->InsertNextPoint(n.coordinates[g * 3 + 0],
                                    n.coordinates[g * 3 + 1],
                                    n.coordinates[g * 3 + 2]);
            localIds[g] = baseIndex + static_cast<vtkIdType>(g);
        }

        const auto& netPanels = byNet[n.netname];
        std::size_t localPanel = 0;
        // A502 panel order: j-major (panel columns), i-fast (panel rows).
        for (int j = 0; j + 1 < nn; ++j) {
            for (int i = 0; i + 1 < nm; ++i) {
                vtkIdType p00 =
                    localIds[static_cast<std::size_t>(j) * nm + i];
                vtkIdType p10 =
                    localIds[static_cast<std::size_t>(j) * nm + (i + 1)];
                vtkIdType p11 =
                    localIds[static_cast<std::size_t>(j + 1) * nm + (i + 1)];
                vtkIdType p01 =
                    localIds[static_cast<std::size_t>(j + 1) * nm + i];
                if (p00 < 0 || p10 < 0 || p11 < 0 || p01 < 0)
                    continue;

                vtkIdType quad[4] = {p00, p10, p11, p01};
                grid->InsertNextCell(VTK_QUAD, 4, quad);
                cellName->InsertNextValue(n.netname.c_str());
                cellNet->InsertNextValue(static_cast<int>(netIdx));

                const PanelResult* pr = nullptr;
                if (localPanel < netPanels.size()) {
                    pr = netPanels[localPanel];
                    ++mapped;
                } else {
                    ++missing;
                }
                ++localPanel;

                if (pr) {
                    cp2nd->InsertNextValue(pr->cp2nd);
                    cpisn->InsertNextValue(pr->cpisn);
                    mach->InsertNextValue(pr->mach);
                    source->InsertNextValue(pr->source);
                    doublet->InsertNextValue(pr->doublet);
                    wvec->InsertNextTuple3(pr->w.x(), pr->w.y(), pr->w.z());
                } else {
                    cp2nd->InsertNextValue(nan);
                    cpisn->InsertNextValue(nan);
                    mach->InsertNextValue(nan);
                    source->InsertNextValue(nan);
                    doublet->InsertNextValue(nan);
                    wvec->InsertNextTuple3(nan, nan, nan);
                }
            }
        }
        baseIndex += static_cast<int>(nGridPts);
    }

    if (missing > 0)
        BLOG(warning) << "Mesh overlay: mapped=" << mapped
                      << " missing Cp for " << missing << " quads";

    grid->GetCellData()->AddArray(cellName);
    grid->GetCellData()->AddArray(cellNet);
    grid->GetCellData()->AddArray(cp2nd);
    grid->GetCellData()->AddArray(cpisn);
    grid->GetCellData()->AddArray(mach);
    grid->GetCellData()->AddArray(source);
    grid->GetCellData()->AddArray(doublet);
    grid->GetCellData()->AddArray(wvec);
    grid->GetCellData()->SetScalars(cp2nd);

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(vtkPath.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToAscii();
    if (writer->Write() == 0) {
        BLOG(error) << "VTK writer failed for \"" << vtkPath << "\"";
        throw std::runtime_error(
            "AeroResult::PrintParaview(mesh): VTK writer failed for \"" +
            vtkPath + "\"");
    }
    BLOG(info) << "Result mesh VTU written: " << vtkPath << " ("
               << grid->GetNumberOfCells() << " cells)";
#else
    (void)vtkPath;
    (void)mesh;
    BLOG(error) << "PrintParaview: built without VTK support (WITH_VTK)";
    throw std::runtime_error(
        "AeroResult::PrintParaview: built without VTK support (WITH_VTK)");
#endif
}

} // namespace a502
