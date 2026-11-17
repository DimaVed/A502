// AeroModel.cpp
// -------------
// Implementation of the A502 input deck parser / model (see AeroModel.h).

#include "AeroModel.h"
#include "FlowSolution.h"
#include "GlobalBoundarySystem.h"
#include "Logger.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#ifdef WITH_VTK
#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkCellTypes.h>
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

// ---------------------------------------------------------------------------
// small string / value helpers
// ---------------------------------------------------------------------------

std::string trim(const std::string& s) {
    const char* ws = " \t\r";
    std::size_t b = s.find_first_not_of(ws);
    if (b == std::string::npos)
        return std::string();
    std::size_t e = s.find_last_not_of(ws);
    return s.substr(b, e - b + 1);
}

bool isNumeric(const std::string& tok) {
    if (tok.empty())
        return false;
    // PanAir numbers may carry a trailing 'D'/'E' exponent marker.
    std::string t = tok;
    char last = t.back();
    if (last == 'D' || last == 'd' || last == 'E' || last == 'e') {
        t.pop_back();
        if (t.empty())
            return false;
    }
    std::size_t pos = 0;
    try {
        (void)std::stod(t, &pos);
    } catch (...) {
        return false;
    }
    return pos == t.size();
}

double toDouble(const std::string& tok) {
    std::string t = tok;
    char last = t.back();
    if (last == 'D' || last == 'd')
        t.pop_back(); // replace 'D' exponent marker with 'E'
    try {
        return std::stod(t);
    } catch (...) {
        return 0.0;
    }
}

Value toValue(const std::string& tok) {
    Value v;
    v.text = tok;
    v.numeric = isNumeric(tok);
    if (v.numeric)
        v.number = toDouble(tok);
    return v;
}

std::vector<std::string> tokenize(const std::string& s) {
    std::vector<std::string> out;
    std::istringstream iss(s);
    std::string tok;
    while (iss >> tok)
        out.push_back(tok);
    return out;
}
// Lower-case and strip every non-alphanumeric character store them in str.
// Used to make "angles of attack", "angles-of-attack" and "angles_of_attack"
// compare equal.
std::string normalizeKeyword(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s) {
        if (std::isalnum(c))
            out.push_back(static_cast<char>(std::tolower(c)));
    }
    return out;
}

// Strip the array subscripts, e.g. "alpha(1)" -> "alpha".
std::string baseName(const std::string& s) {
    std::size_t p = s.find('(');
    return p == std::string::npos ? s : s.substr(0, p);
}

// Lower-case copy used for A502's case-insensitive name comparisons.
std::string lowerCopy(const std::string& s) {
    std::string out;
    out.reserve(s.size());
    for (unsigned char c : s)
        out.push_back(static_cast<char>(std::tolower(c)));
    return out;
}

// Returns the first numeric Value of an Entry, or nullptr.
const Value* firstNumber(const Entry& e) {
    for (const auto& v : e.values)
        if (v.numeric)
            return &v;
    return nullptr;
}

// Collect every numeric Value of an Entry in order.
std::vector<double> entryNumbers(const Entry& e) {
    std::vector<double> out;
    for (const auto& v : e.values)
        if (v.numeric)
            out.push_back(v.number);
    return out;
}

// Append numeric (and string) tokens into a destination vector of Values.
void appendTokens(std::vector<Value>& dst, const std::vector<std::string>& toks) {
    for (const auto& t : toks)
        dst.push_back(toValue(t));
}

// Append raw numeric tokens only (skip non-numeric words).
void appendNumbers(std::vector<double>& dst, const std::vector<Value>& src) {
    for (const auto& v : src)
        if (v.numeric)
            dst.push_back(v.number);
}

} // anonymous namespace
// ===========================================================================
// public API
// ===========================================================================

void AeroModel::clear() {
    titleLines.clear();
    ndtchk = 0;
    misym = mjsym = 0;
    amach = betc = 0.0;
    nacase = 1;
    alpc = 0.0;
    alpha.clear();
    beta.clear();
    printoutOptions.clear();
    refs = ReferenceData{};
    eat = 0.0;
    ivcorr = 0;
    nflowv = tpoff = 0.0;
    offBodyPoints.clear();
    sections.clear();
    networks.clear();
    materials.assign(1, MaterialProperties{});
    surfaceMaterials.clear();
    flow = Compressibility{};
    panels.clear();
    controlPoints.clear();
    totalSurfaceArea = 0.0;
    totalPanels = 0;
    solved_ = false;
    flowSolved_ = false;
}

void AeroModel::parse(const std::string& inpPath) {
    clear();

    BLOG(info) << "Parsing A502 input deck: " << inpPath;
    std::ifstream in(inpPath);
    if (!in.is_open()) {
        BLOG(error) << "Cannot open input file \"" << inpPath << "\"";
        throw std::runtime_error("AeroModel::parse: cannot open input file \"" +
                                 inpPath + "\"");
    }

    Section     current;
    bool        haveSection = false;
    Entry*      curEntry = nullptr;

    std::string line;
    while (std::getline(in, line)) {
        std::string t = trim(line);
        if (t.empty())
            continue;

        char c = t[0];
        if (c == '$') {
            // ---- start of a new section ----
            if (haveSection)
                sections.push_back(std::move(current));
            current = Section{};
            haveSection = true;

            std::string rest = trim(t.substr(1)); // everything after '$'
            std::vector<std::string> toks = tokenize(rest);
            std::string kw = toks.empty() ? std::string() : toks[0];
            current.keyword = kw;

            std::size_t space = rest.find_first_of(" \t");
            current.description = (space == std::string::npos)
                                      ? std::string()
                                      : trim(rest.substr(space + 1));
            current.keywordNorm = normalizeKeyword(kw);
            curEntry = nullptr;
        } else if (c == '*') {
            // ---- "*" marker card, kept as a Star entry ----
            if (haveSection) {
                Entry e;
                e.kind = EntryKind::Star;
                e.text = t;
                current.entries.push_back(std::move(e));
                curEntry = &current.entries.back();
            }
        } else if (c == '=') {
// ---- "=name ..." declaration or a value-on-the-equals card ----
            if (!haveSection)
                continue;
            std::string rest = trim(t.substr(1));
            std::vector<std::string> toks = tokenize(rest);

            // A declaration card contains at least one alphabetic name.
            bool anyName = false;
            for (const auto& tok : toks)
                if (!isNumeric(tok)) {
                    anyName = true;
                    break;
                }

            if (anyName && !toks.empty()) {
                Entry e;
                e.kind = EntryKind::Variable;
                e.text = rest;
                for (const auto& tok : toks)
                    e.names.push_back(baseName(tok));
                current.entries.push_back(std::move(e));
                curEntry = &current.entries.back();
            } else {
                // Old format (e.g. ellip.inp) puts the value right after '='.
                // Treat it as plain data appended to the current group.
                std::vector<Value> vals;
                appendTokens(vals, toks);
                if (curEntry != nullptr &&
                    (curEntry->kind == EntryKind::Variable ||
                     curEntry->kind == EntryKind::Star)) {
                    curEntry->values.insert(curEntry->values.end(), vals.begin(),
                                            vals.end());
                } else {
                    Entry e;
                    e.kind = EntryKind::FreeData;
                    e.text = t;
                    e.values = std::move(vals);
                    current.entries.push_back(std::move(e));
                    curEntry = &current.entries.back();
                }
            }
        } else if (c == '#') {
            // ---- A502 "#" comment / index marker: keep as data ----
            if (!haveSection)
                continue;
            std::vector<Value> vals;
            appendTokens(vals, tokenize(t.substr(1)));
            if (curEntry != nullptr &&
                (curEntry->kind == EntryKind::Variable ||
                 curEntry->kind == EntryKind::Star)) {
                curEntry->values.insert(curEntry->values.end(), vals.begin(),
                                        vals.end());
            } else {
                Entry e;
                e.kind = EntryKind::FreeData;
                e.text = t;
                e.values = std::move(vals);
                current.entries.push_back(std::move(e));
                curEntry = &current.entries.back();
            }
        } else {
            // ---- ordinary data card ----
            if (!haveSection) {
                // Data before the first '$' section.
                continue;
            }
            Entry e;
            e.kind = EntryKind::FreeData;
            e.text = t;
            appendTokens(e.values, tokenize(t));
            current.entries.push_back(std::move(e));
            curEntry = &current.entries.back();
        }
    }
    if (haveSection)
        sections.push_back(std::move(current));

    extractParameters();
    extractNetworks();
    extractMaterials();
    BLOG(info) << "Parsed " << sections.size() << " sections, "
               << networks.size() << " networks"
               << " (mach=" << amach << ", nacase=" << nacase << ")";
}
void AeroModel::parse(const std::string& inpPath, const std::string& outPath) {
    parse(inpPath);
    writeOut(outPath);
}

void AeroModel::writeOut(const std::string& outPath) const {
    BLOG(info) << "Writing model report: " << outPath;
    std::ofstream out(outPath);
    if (!out.is_open()) {
        BLOG(error) << "Cannot create report file \"" << outPath << "\"";
        throw std::runtime_error("AeroModel::writeOut: cannot create file \"" +
                                 outPath + "\"");
    }

    out << std::fixed << std::setprecision(5);
    writeHeader(out);

    // ---- scalar parameters ----
    out << "\n::: SCALAR PARAMETERS :::\n";
    out << "  ndtchk       = " << ndtchk << "\n";
    out << "  symmetry     = misym=" << misym << "  mjsym=" << mjsym << "\n";
    out << "  amach        = " << amach << "\n";
    out << "  betc         = " << betc << "\n";
    out << "  nacase       = " << nacase << "\n";
    out << "  alpc         = " << alpc << "\n";
    out << "  alpha[" << alpha.size() << "] =";
    for (double a : alpha) out << " " << a;
    out << "\n  beta[" << beta.size() << "] =";
    for (double b : beta) out << " " << b;
    out << "\n  eat          = " << eat << "\n";
    out << "  ivcorr       = " << ivcorr << "\n";
    out << "  nflowv,tpoff = " << nflowv << " " << tpoff << "\n";
    out << "  reference    = xref=" << refs.xref
        << " yref=" << refs.yref << " zref=" << refs.zref
        << " nref=" << refs.nref << " nprcof=" << refs.nprcof << "\n";
    out << "  reference    = sref=" << refs.sref
        << " bref=" << refs.bref << " cref=" << refs.cref
        << " dref=" << refs.dref << "\n";
    if (!printoutOptions.empty()) {
        out << "  printout[" << printoutOptions.size() << "] =";
        for (double p : printoutOptions) out << " " << p;
        out << "\n";
    }

    // ---- networks ----
    out << "\n::: NETWORKS (" << networks.size() << ") :::\n";
    if (networks.empty())
        out << "  (none)\n";
    for (const auto& n : networks) {
        out << "  net '" << n.netname << "'  kt=" << n.kt
            << "  nm=" << n.nm << "  nn=" << n.nn
            << "  cpnorm=" << n.cpnorm
            << "  coords=" << n.coordinates.size() << "\n";
        if (!n.coordinates.empty()) {
            double xmin = 1e30, xmax = -1e30, ymin = 1e30, ymax = -1e30,
                   zmin = 1e30, zmax = -1e30;
            for (std::size_t i = 0; i + 2 < n.coordinates.size(); i += 3) {
                xmin = std::min(xmin, n.coordinates[i]);
                xmax = std::max(xmax, n.coordinates[i]);
                ymin = std::min(ymin, n.coordinates[i + 1]);
                ymax = std::max(ymax, n.coordinates[i + 1]);
                zmin = std::min(zmin, n.coordinates[i + 2]);
                zmax = std::max(zmax, n.coordinates[i + 2]);
            }
            out << "      bbox x[" << xmin << ", " << xmax << "]"
                << " y[" << ymin << ", " << ymax << "]"
                << " z[" << zmin << ", " << zmax << "]\n";
        }
    }

    // ---- materials ----
    out << "\n::: MATERIALS (" << (materials.size() - 1) << ") :::\n";
    for (std::size_t i = 0; i < materials.size(); ++i) {
        out << "  [" << i << "] '" << materials[i].name
            << "'  tratio=" << materials[i].temperatureRatio
            << "  pratio=" << materials[i].pressureRatio << "\n";
    }
    for (const auto& sm : surfaceMaterials) {
        out << "  surface material: net '" << sm.network << "' "
            << (sm.surface > 0 ? "upper" : sm.surface < 0 ? "lower"
                                                            : "unknown")
            << " -> material " << sm.material << "\n";
    }

    // ---- full section listing ----
    out << "\n::: SECTIONS (" << sections.size() << ") :::\n";
    for (const auto& sec : sections) {
        out << "\n  $" << sec.keyword;
        if (!sec.description.empty())
            out << "  -- " << sec.description;
        out << "\n";
        for (const auto& e : sec.entries) {
            switch (e.kind) {
            case EntryKind::Variable:
                out << "    = ";
                for (std::size_t k = 0; k < e.names.size(); ++k)
                    out << (k ? " " : "") << e.names[k];
                break;
            case EntryKind::Star:
                out << "    *  " << e.text;
                break;
            case EntryKind::FreeData:
                out << "    .  " << e.text;
                break;
            }
            if (!e.values.empty()) {
                out << "   [[ ";
                for (const auto& v : e.values)
                    out << (v.numeric ? std::to_string(v.number) : v.text) << " ";
                out << "]]";
            }
            out << "\n";
        }
    }
    out << "\n... end of model report ...\n";
    BLOG(info) << "Model report written: " << outPath;
}

void AeroModel::printSummary() const {
    if (titleLines.empty()) {
        BLOG(info) << "title: (none)";
    } else {
        BLOG(info) << "title:";
        for (const auto& line : titleLines)
            BLOG(info) << "  " << line;
    }

    BLOG(info) << "mach=" << amach
               << "  nacase=" << nacase
               << "  alpc=" << alpc
               << "  misym=" << misym
               << "  mjsym=" << mjsym;

    BLOG(info) << "networks: " << networks.size();
    for (std::size_t i = 0; i < networks.size(); ++i) {
        const Network& n = networks[i];
        const std::size_t nPts = n.coordinates.size() / 3;
        const int nCells =
            (n.nm > 1 && n.nn > 1) ? (n.nm - 1) * (n.nn - 1) : 0;
        BLOG(info) << "  [" << i << "] '" << n.netname << "'"
                   << "  kt=" << n.kt
                   << "  nm=" << n.nm << " nn=" << n.nn
                   << "  points=" << nPts
                   << "  quads~=" << nCells
                   << "  (section kn=" << n.kn << ")";
    }
}

void AeroModel::PrintParaview(const std::string& vtkPath) const {
#ifdef WITH_VTK
    BLOG(info) << "Writing ParaView VTU: " << vtkPath;
    if (networks.empty()) {
        BLOG(error) << "PrintParaview: model contains no grid networks";
        throw std::runtime_error(
            "AeroModel::PrintParaview: model contains no grid networks");
    }

    vtkNew<vtkPoints>           points;
    vtkNew<vtkUnstructuredGrid> grid;
    vtkNew<vtkIntArray>         pointNet;      // network index, per point
    vtkNew<vtkStringArray>      cellName;      // netname, per cell
    vtkNew<vtkIntArray>         cellNet;

    pointNet->SetName("network_id");
    cellName->SetName("netname");
    cellNet->SetName("network_id");
    grid->SetPoints(points);

    // Merge every grid network into one unstructured grid.  Each network's
    // point array is stored column-major: point(i,j) = coordinates[(j*nm+i)*3].
    int baseIndex = 0; // first global point index of the current network
    for (std::size_t netIdx = 0; netIdx < networks.size(); ++netIdx) {
        const Network& n = networks[netIdx];
        int nm = n.nm > 0 ? n.nm : 0;
        int nn = n.nn > 0 ? n.nn : 0;
        if (nm < 1 || nn < 1)
            continue;

        std::size_t nAvail = n.coordinates.size() / 3;
        std::size_t nWanted = static_cast<std::size_t>(nm) * nn;
        std::size_t nGridPts = std::min(nWanted, nAvail);
        if (nGridPts < 1)
            continue;

        // ---- insert the points ----
        std::vector<vtkIdType> localIds(static_cast<std::size_t>(nm) * nn,
                                        -1);
        for (std::size_t g = 0; g < nGridPts; ++g) {
            double x = n.coordinates[g * 3 + 0];
            double y = n.coordinates[g * 3 + 1];
            double z = n.coordinates[g * 3 + 2];
            points->InsertNextPoint(x, y, z);
            pointNet->InsertNextValue(static_cast<int>(netIdx));
            localIds[g] = baseIndex + static_cast<vtkIdType>(g);
        }

        // ---- build quads (and degenerate lines / vertices) ----
        for (int i = 0; i < nm; ++i)
            for (int j = 0; j < nn; ++j) {
                vtkIdType p = localIds[static_cast<std::size_t>(j) * nm + i];
                if (p < 0)
                    continue;
                bool hasRight = (i + 1 < nm) &&
                                localIds[static_cast<std::size_t>(j) * nm +
                                         (i + 1)] >= 0;
                bool hasUp = (j + 1 < nn) &&
                             localIds[static_cast<std::size_t>(j + 1) * nm +
                                      i] >= 0;
                if (hasRight && hasUp) {
                    vtkIdType quad[4] = {
                        p,
                        localIds[static_cast<std::size_t>(j) * nm + (i + 1)],
                        localIds[static_cast<std::size_t>(j + 1) * nm + (i + 1)],
                        localIds[static_cast<std::size_t>(j + 1) * nm + i]};
                    grid->InsertNextCell(VTK_QUAD, 4, quad);
                    cellName->InsertNextValue(n.netname.c_str());
                    cellNet->InsertNextValue(static_cast<int>(netIdx));
                } else if (hasRight) {
                    vtkIdType line[2] = {
                        p, localIds[static_cast<std::size_t>(j) * nm + (i + 1)]};
                    grid->InsertNextCell(VTK_LINE, 2, line);
                    cellName->InsertNextValue(n.netname.c_str());
                    cellNet->InsertNextValue(static_cast<int>(netIdx));
                } else if (hasUp) {
                    vtkIdType line[2] = {
                        p, localIds[static_cast<std::size_t>(j + 1) * nm + i]};
                    grid->InsertNextCell(VTK_LINE, 2, line);
                    cellName->InsertNextValue(n.netname.c_str());
                    cellNet->InsertNextValue(static_cast<int>(netIdx));
                } else if (nGridPts == 1) {
                    vtkIdType vtx[1] = {p};
                    grid->InsertNextCell(VTK_VERTEX, 1, vtx);
                    cellName->InsertNextValue(n.netname.c_str());
                    cellNet->InsertNextValue(static_cast<int>(netIdx));
                }
            }

        baseIndex += static_cast<int>(nGridPts);
    }

    grid->GetPointData()->AddArray(pointNet);
    grid->GetCellData()->AddArray(cellName);
    grid->GetCellData()->AddArray(cellNet);

    vtkNew<vtkXMLUnstructuredGridWriter> writer;
    writer->SetFileName(vtkPath.c_str());
    writer->SetInputData(grid);
    writer->SetDataModeToAscii();
    if (writer->Write() == 0) {
        BLOG(error) << "VTK writer failed for \"" << vtkPath << "\"";
        throw std::runtime_error(
            "AeroModel::PrintParaview: VTK writer failed for \"" + vtkPath +
            "\"");
    }
    BLOG(info) << "ParaView VTU written: " << vtkPath
               << " (" << grid->GetNumberOfPoints() << " points, "
               << grid->GetNumberOfCells() << " cells)";
#else
    (void)vtkPath;
    BLOG(error) << "PrintParaview: built without VTK support (WITH_VTK)";
    throw std::runtime_error(
        "AeroModel::PrintParaview: built without VTK support (WITH_VTK)");
#endif
}
const Section* AeroModel::findSection(const std::string& kwFragment) const {
    std::string frag = normalizeKeyword(kwFragment);
    for (const auto& s : sections)
        if (s.keywordNorm.find(frag) != std::string::npos)
            return &s;
    return nullptr;
}

void AeroModel::writeHeader(std::ostream& out) const {
    out << "==================================================================\n";
    out << "  A502 / PANAIR input deck  --  model report (AeroModel)\n";
    out << "==================================================================\n";
    if (!titleLines.empty()) {
        out << "  title:";
        for (const auto& tl : titleLines)
            out << "\n    " << tl;
        out << "\n";
    }
    out << "==================================================================\n";
}

// ===========================================================================
// typed parameter extraction
// ===========================================================================

void AeroModel::extractParameters() {
    // ----- title -----
    if (const Section* t = findSection("title")) {
        for (const auto& e : t->entries)
            if (e.kind == EntryKind::FreeData && !e.text.empty())
                titleLines.push_back(e.text);
    }

    // ----- datacheck -----
    if (const Section* s = findSection("datacheck")) {
        for (const auto& e : s->entries)
            if (const Value* v = firstNumber(e)) {
                ndtchk = static_cast<int>(std::lround(v->number));
                break;
            }
    }

    // ----- symmetry -----
    if (const Section* s = findSection("sym")) { // matches $sym, $symmetric ...
        int got = 0;
        for (const auto& e : s->entries) {
            for (const auto& v : e.values) {
                if (!v.numeric)
                    continue;
                if (got == 0)
                    misym = static_cast<int>(std::lround(v.number));
                else if (got == 1)
                    mjsym = static_cast<int>(std::lround(v.number));
                ++got;
                if (got >= 2)
                    break;
            }
            if (got >= 2)
                break;
        }
    }

    // ----- mach number -----
    if (const Section* s = findSection("mach")) {
        for (const auto& e : s->entries)
            if (const Value* v = firstNumber(e)) {
                amach = v->number;
                break;
            }
    }

    // ----- cases -----
    if (const Section* s = findSection("cases")) {
        for (const auto& e : s->entries)
            if (const Value* v = firstNumber(e)) {
                nacase = static_cast<int>(std::lround(v->number));
                break;
            }
    }

    // ----- angles of attack (alpc + alpha array) -----
    if (const Section* s = findSection("angles")) {
        std::vector<double> vals;
        for (const auto& e : s->entries)
            appendNumbers(vals, e.values);
        if (!vals.empty()) {
            alpc = vals[0];
            alpha.assign(vals.begin() + 1, vals.end());
        }
    }

    // ----- yaw / sideslip angles (betc + beta array) -----
    if (const Section* s = findSection("yaw")) {
        std::vector<double> vals;
        for (const auto& e : s->entries)
            appendNumbers(vals, e.values);
        if (!vals.empty()) {
            betc = vals[0];
            beta.assign(vals.begin() + 1, vals.end());
        }
    }

    // ----- references -----
    if (const Section* s = findSection("refer")) {
        for (const auto& e : s->entries) {
            if (e.kind != EntryKind::Variable || e.names.empty())
                continue;
            std::vector<double> vals = entryNumbers(e);
            std::string first = e.names.front();
            if (first == "xref" || first == "xoff") {
                refs.xref = vals.size() > 0 ? vals[0] : 0.0;
                refs.yref = vals.size() > 1 ? vals[1] : 0.0;
                refs.zref = vals.size() > 2 ? vals[2] : 0.0;
                if (vals.size() > 3)
                    refs.nref = static_cast<int>(std::lround(vals[3]));
            } else if (first == "sref") {
                refs.sref = vals.size() > 0 ? vals[0] : 0.0;
                refs.bref = vals.size() > 1 ? vals[1] : 0.0;
                refs.cref = vals.size() > 2 ? vals[2] : 0.0;
                refs.dref = vals.size() > 3 ? vals[3] : 0.0;
            } else if (first == "nref") {
                if (!vals.empty())
                    refs.nref = static_cast<int>(std::lround(vals[0]));
            }
        }
    }

    // ----- printout options -----
    if (const Section* s = findSection("printout")) {
        for (const auto& e : s->entries)
            appendNumbers(printoutOptions, e.values);
    }

    // ----- eat -----
    if (const Section* s = findSection("eat"))
        for (const auto& e : s->entries)
            if (const Value* v = firstNumber(e)) {
                eat = v->number;
                break;
            }

    // ----- velocity correction -----
    if (const Section* s = findSection("velocity"))
        for (const auto& e : s->entries)
            if (const Value* v = firstNumber(e)) {
                ivcorr = static_cast<int>(std::lround(v->number));
                break;
            }

    // ----- flowfield properties -----
    if (const Section* s = findSection("flow")) {
        int got = 0;
        for (const auto& e : s->entries) {
            for (const auto& v : e.values) {
                if (!v.numeric)
                    continue;
                if (got == 0)
                    nflowv = v.number;
                else if (got == 1)
                    tpoff = v.number;
                ++got;
                if (got >= 2)
                    break;
            }
            if (got >= 2)
                break;
        }
    }
}
// ===========================================================================
// grid network extraction ($points sections)
// ===========================================================================

namespace {

bool entryHasName(const Entry& e, const char* name) {
    for (const auto& n : e.names)
        if (n == name)
            return true;
    return false;
}

bool entryAllNumeric(const Entry& e) {
    if (e.values.empty())
        return false;
    for (const auto& v : e.values)
        if (!v.numeric)
            return false;
    return true;
}

// Card A of a $points network: "nm nn ... netname" (6E10.0,10X,A10).
// The netname is the only non-numeric token on the card.
bool parseNetworkCard(const Entry& e, int& nm, int& nn, std::string& netname) {
    if (e.values.size() < 3)
        return false;
    if (!e.values[0].numeric || !e.values[1].numeric)
        return false;
    if (e.values.back().numeric)
        return false;

    nm = static_cast<int>(std::lround(e.values[0].number));
    nn = static_cast<int>(std::lround(e.values[1].number));
    netname = e.values.back().text;
    return nm > 0 && nn > 0 && !netname.empty();
}

} // namespace

void AeroModel::extractNetworks() {
    for (const auto& sec : sections) {
        if (sec.keywordNorm == "points")
            extractPointsSection(sec);
        else if (sec.keywordNorm == "circular" || sec.keywordNorm == "cir")
            extractCircularSection(sec);
        else if (sec.keywordNorm == "quadrilateral" ||
                 sec.keywordNorm == "qua")
            extractQuadrilateralSection(sec);
        else if (sec.keywordNorm == "elliptic" ||
                 sec.keywordNorm == "ell")
            extractEllipticalSection(sec);
        else if (sec.keywordNorm == "trailing" || sec.keywordNorm == "tra")
            extractTrailingSection(sec);
    }
}

void AeroModel::extractPointsSection(const Section& sec) {
    // The original reader (meshp / inbc) consumes a $points section as:
    //
    //   kn card      : 6E10.0 (field 1 = kn, field 6 = cpnorm)
    //   kt card      : 6E10.0 (field 1 = kt, field 5 = amnsw, field 6 = dnsmsh)
    //   network card : "nm nn ... netname"
    //   coordinates  : ceil(nm/2) cards per column, 2 points per card
    //
    // Modern decks spell the first two cards as "=kn" / "=kt" label cards
    // followed by their value cards.  Older decks (ellip.inp) put the value
    // cards directly.  A single $points section may contain several networks;
    // kn/kt then apply to all of them.
    {
        int    sectionKn = 0;
        int    sectionKt = 0;
        double sectionCpnorm = 2.0;

        enum class Header { NeedKn, NeedKt, Done };
        Header header = Header::NeedKn;
        bool pendingKn = false;
        bool pendingKt = false;

        for (std::size_t ei = 0; ei < sec.entries.size(); ++ei) {
            const Entry& e = sec.entries[ei];

            if (e.kind == EntryKind::Variable) {
                if (entryHasName(e, "kn"))
                    pendingKn = true;
                else if (entryHasName(e, "kt"))
                    pendingKt = true;
                // Other "=..." labels (nm/nn/netname, bet1, ...) are skipped.
                continue;
            }
            if (e.kind != EntryKind::FreeData)
                continue;

            // ---- kn card -------------------------------------------------
            if (header == Header::NeedKn || pendingKn) {
                if (!entryAllNumeric(e))
                    continue;
                sectionKn = static_cast<int>(
                    std::lround(e.values.front().number));
                // The cpnorm value sits in the sixth fixed-width field.  After
                // tokenization the last number is the best available proxy.
                if (e.values.size() >= 2)
                    sectionCpnorm = e.values.back().number;
                header = Header::NeedKt;
                pendingKn = false;
                continue;
            }

            // ---- kt card -------------------------------------------------
            if (header == Header::NeedKt || pendingKt) {
                if (!e.values.empty() && e.values.front().numeric)
                    sectionKt = static_cast<int>(
                        std::lround(e.values.front().number));
                header = Header::Done;
                pendingKt = false;
                continue;
            }

            // ---- network card -------------------------------------------
            int nm = 0, nn = 0;
            std::string netname;
            if (!parseNetworkCard(e, nm, nn, netname))
                continue;

            Network net;
            net.netname = netname;
            net.nm = nm;
            net.nn = nn;
            net.kn = sectionKn;
            net.kt = sectionKt;
            net.cpnorm = sectionCpnorm;

            const std::size_t needed =
                static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn) * 3;
            std::size_t got = 0;

            // Read exactly nm*nn points.  Cards normally hold two points
            // (six numbers); the last card of an odd column holds one point.
            // Label cards ("=x y z ...") may sit between the network card and
            // its coordinates, so skip Variable/Star entries.
            std::size_t cj = ei + 1;
            while (cj < sec.entries.size() && got < needed) {
                const Entry& ce = sec.entries[cj];
                if (ce.kind == EntryKind::Variable ||
                    ce.kind == EntryKind::Star) {
                    ++cj;
                    continue;
                }
                if (ce.kind != EntryKind::FreeData || !entryAllNumeric(ce))
                    break; // next network card or a control card
                for (const auto& v : ce.values) {
                    net.coordinates.push_back(v.number);
                    ++got;
                    if (got == needed)
                        break;
                }
                ++cj;
            }

            if (got != needed) {
                BLOG(warning) << "network '" << net.netname << "' (nm=" << nm
                              << ", nn=" << nn << ") has " << (got / 3)
                              << " coordinate points, expected " << nm * nn;
            }

            networks.push_back(std::move(net));
            ei = cj - 1; // continue after the coordinate cards
        }
    }
}

void AeroModel::extractCircularSection(const Section& sec) {
    // Port of circ(): a surface of revolution generated from a meridional
    // profile (xs, r) and a list of angular positions th (degrees):
    //   x = xs(i), y = r(i)*cos(th), z = z(i) + r(i)*sin(th).
    // A section may contain kn network blocks, each introduced by an
    // "=nopt netname" card.  nopt > 0 expects an "=z" card with an axial
    // offset per station; the raw-card form without a "=z" label is not
    // handled yet.
    struct Block {
        bool started{false};
        int nopt{0};
        std::string name;
        int nm{0};
        int nn{0};
        std::vector<double> xs, r, z, th;
    };

    int sectionKn = 0;
    int sectionKt = 0;
    double sectionCpnorm = 2.0;
    bool pendingKn = false;
    bool pendingKt = false;

    auto groupValues = [&](std::size_t index) {
        std::vector<Value> values = sec.entries[index].values;
        for (std::size_t j = index + 1;
             j < sec.entries.size() &&
             sec.entries[j].kind == EntryKind::FreeData;
             ++j) {
            values.insert(values.end(),
                          sec.entries[j].values.begin(),
                          sec.entries[j].values.end());
        }
        return values;
    };
    auto numbers = [](const std::vector<Value>& values) {
        std::vector<double> out;
        for (const Value& v : values)
            if (v.numeric)
                out.push_back(v.number);
        return out;
    };

    Block current;
    auto finishBlock = [&]() {
        if (!current.started)
            return;
        const bool complete =
            !current.name.empty() && current.nm > 0 && current.nn > 0 &&
            static_cast<int>(current.xs.size()) >= current.nm &&
            static_cast<int>(current.r.size()) >= current.nm &&
            static_cast<int>(current.th.size()) >= current.nn;
        if (!complete) {
            BLOG(warning) << "circular network '" << current.name
                          << "' is incomplete";
            current = Block{};
            return;
        }

        Network net;
        net.netname = current.name;
        net.nm = current.nm;
        net.nn = current.nn;
        net.kn = sectionKn;
        net.kt = sectionKt;
        net.cpnorm = sectionCpnorm;
        if (static_cast<int>(current.z.size()) < current.nm)
            current.z.assign(static_cast<std::size_t>(current.nm), 0.0);

        constexpr double pi = 3.14159265358979323846;
        const double deg = pi / 180.0;
        net.coordinates.reserve(
            static_cast<std::size_t>(3) * current.nm * current.nn);
        for (int j = 0; j < current.nn; ++j) {
            const double st = std::sin(current.th[static_cast<std::size_t>(j)] * deg);
            const double ct = std::cos(current.th[static_cast<std::size_t>(j)] * deg);
            for (int i = 0; i < current.nm; ++i) {
                const std::size_t si = static_cast<std::size_t>(i);
                net.coordinates.push_back(current.xs[si]);
                net.coordinates.push_back(current.r[si] * ct);
                net.coordinates.push_back(current.z[si] +
                                          current.r[si] * st);
            }
        }
        networks.push_back(std::move(net));
        current = Block{};
    };

    for (std::size_t ei = 0; ei < sec.entries.size(); ++ei) {
        const Entry& e = sec.entries[ei];
        if (e.kind == EntryKind::Variable && !e.names.empty()) {
            const std::string& first = e.names.front();
            if (first == "kn") {
                pendingKn = true;
            } else if (first == "kt") {
                pendingKt = true;
            } else if (first == "nopt") {
                finishBlock();
                const std::vector<Value> values = groupValues(ei);
                current.started = true;
                for (const Value& v : values) {
                    if (v.numeric) {
                        current.nopt = static_cast<int>(
                            std::lround(v.number));
                        break;
                    }
                }
                for (const Value& v : values) {
                    if (!v.numeric) {
                        current.name = v.text;
                        break;
                    }
                }
            } else if (first == "nm") {
                const std::vector<double> values = numbers(groupValues(ei));
                if (!values.empty())
                    current.nm = static_cast<int>(
                        std::lround(values.front()));
            } else if (first == "xs") {
                const std::vector<double> values = numbers(groupValues(ei));
                current.xs.clear();
                current.r.clear();
                for (std::size_t k = 0; k + 1 < values.size(); k += 2) {
                    current.xs.push_back(values[k]);
                    current.r.push_back(values[k + 1]);
                }
            } else if (first == "z") {
                current.z = numbers(groupValues(ei));
            } else if (first == "nn") {
                const std::vector<double> values = numbers(groupValues(ei));
                if (!values.empty())
                    current.nn = static_cast<int>(
                        std::lround(values.front()));
            } else if (first == "th") {
                current.th = numbers(groupValues(ei));
            }
            continue;
        }
        if (e.kind != EntryKind::FreeData)
            continue;

        if (pendingKn) {
            const std::vector<double> values = numbers(e.values);
            if (!values.empty()) {
                sectionKn = static_cast<int>(
                    std::lround(values.front()));
                if (values.size() >= 2)
                    sectionCpnorm = values.back();
            }
            pendingKn = false;
        } else if (pendingKt) {
            const std::vector<double> values = numbers(e.values);
            if (!values.empty())
                sectionKt = static_cast<int>(
                    std::lround(values.front()));
            pendingKt = false;
        }
    }
    finishBlock();
}

void AeroModel::extractQuadrilateralSection(const Section& sec) {
    // Port of quadnt(): bilinear map from the unit square to the four given
    // corner points, with row/column division fractions ypc/xpc.  The corner
    // cards follow the kt card directly ("6E10.0,10X,A" for corners 1-2 and
    // the network name; the second card holds corners 3-4).  Labeled division
    // groups are "=nrow", "=ypc", "=ncol", "=xpc".
    int sectionKn = 0;
    int sectionKt = 0;
    double sectionCpnorm = 2.0;
    bool pendingKn = false;
    bool pendingKt = false;
    bool afterKt = false;
    std::array<Eigen::Vector3d, 4> corners{};
    int cornerCount = 0;
    std::string netname;
    int nrow = 0;
    int ncol = 0;
    std::vector<double> ypc;
    std::vector<double> xpc;

    auto groupValues = [&](std::size_t index) {
        std::vector<Value> values = sec.entries[index].values;
        for (std::size_t j = index + 1;
             j < sec.entries.size() &&
             sec.entries[j].kind == EntryKind::FreeData;
             ++j) {
            values.insert(values.end(),
                          sec.entries[j].values.begin(),
                          sec.entries[j].values.end());
        }
        return values;
    };
    auto numbers = [](const std::vector<Value>& values) {
        std::vector<double> out;
        for (const Value& v : values)
            if (v.numeric)
                out.push_back(v.number);
        return out;
    };

    for (std::size_t ei = 0; ei < sec.entries.size(); ++ei) {
        const Entry& e = sec.entries[ei];
        if (e.kind == EntryKind::Variable && !e.names.empty()) {
            const std::string& first = e.names.front();
            if (first == "kn") {
                pendingKn = true;
            } else if (first == "kt") {
                pendingKt = true;
            } else if (first == "nrow") {
                afterKt = false;
                const std::vector<double> values = numbers(groupValues(ei));
                if (!values.empty())
                    nrow = static_cast<int>(
                        std::lround(values.front()));
            } else if (first == "ypc") {
                afterKt = false;
                ypc = numbers(groupValues(ei));
            } else if (first == "ncol") {
                afterKt = false;
                const std::vector<double> values = numbers(groupValues(ei));
                if (!values.empty())
                    ncol = static_cast<int>(
                        std::lround(values.front()));
            } else if (first == "xpc") {
                afterKt = false;
                xpc = numbers(groupValues(ei));
            }
            continue;
        }
        if (e.kind != EntryKind::FreeData)
            continue;

        if (pendingKn) {
            const std::vector<double> values = numbers(e.values);
            if (!values.empty()) {
                sectionKn = static_cast<int>(
                    std::lround(values.front()));
                if (values.size() >= 2)
                    sectionCpnorm = values.back();
            }
            pendingKn = false;
        } else if (pendingKt) {
            const std::vector<double> values = numbers(e.values);
            if (!values.empty())
                sectionKt = static_cast<int>(
                    std::lround(values.front()));
            pendingKt = false;
            afterKt = true;
        } else if (afterKt && cornerCount < 4) {
            const std::vector<double> values = numbers(e.values);
            if (values.size() < 6)
                continue;
            corners[static_cast<std::size_t>(cornerCount)] =
                Eigen::Vector3d(values[0], values[1], values[2]);
            corners[static_cast<std::size_t>(cornerCount + 1)] =
                Eigen::Vector3d(values[3], values[4], values[5]);
            cornerCount += 2;
            for (const Value& v : e.values) {
                if (!v.numeric) {
                    netname = v.text;
                    break;
                }
            }
        }
    }

    if (cornerCount != 4 || nrow < 2 || ncol < 2 ||
        static_cast<int>(ypc.size()) < nrow ||
        static_cast<int>(xpc.size()) < ncol || netname.empty()) {
        BLOG(warning) << "$quadrilateral: incomplete network definition";
        return;
    }

    Network net;
    net.netname = netname;
    net.nm = nrow;
    net.nn = ncol;
    net.kn = sectionKn;
    net.kt = sectionKt;
    net.cpnorm = sectionCpnorm;
    net.coordinates.assign(
        static_cast<std::size_t>(3) * nrow * ncol, 0.0);

    const Eigen::Vector3d& c1 = corners[0];
    const Eigen::Vector3d& c2 = corners[1];
    const Eigen::Vector3d& c3 = corners[2];
    const Eigen::Vector3d& c4 = corners[3];

    for (int j = 0; j < ncol; ++j) {
        const double eps = xpc[static_cast<std::size_t>(j)];
        const double ceps = 1.0 - eps;
        for (int i = 0; i < nrow; ++i) {
            const double eta = ypc[static_cast<std::size_t>(i)];
            const double ceta = 1.0 - eta;
            const double ea = ceps * ceta;
            const double eb = eps * ceta;
            const double ec = eps * eta;
            const double ed = ceps * eta;
            const std::size_t base =
                static_cast<std::size_t>(3) *
                static_cast<std::size_t>(i + nrow * j);
            net.coordinates[base] =
                c1.x() * ea + c2.x() * eb + c3.x() * ec + c4.x() * ed;
            net.coordinates[base + 1] =
                c1.y() * ea + c2.y() * eb + c3.y() * ec + c4.y() * ed;
        }
    }

    const bool anyZ = c1.z() != 0.0 || c2.z() != 0.0 ||
                      c3.z() != 0.0 || c4.z() != 0.0;
    if (anyZ) {
        // Slant-wing z interpolation from quadnt().
        const double zas = c2.z() - c1.z();
        const double zbs = c3.z() - c4.z();
        const double da =
            std::hypot(c2.x() - c1.x(), c2.y() - c1.y());
        const double db =
            std::hypot(c3.x() - c4.x(), c3.y() - c4.y());
        const double zda = (da != 0.0) ? zas / da : 0.0;
        const double zdb = (db != 0.0) ? zbs / db : 0.0;
        const int nrm1 = nrow - 1;
        for (int j = 0; j < ncol; ++j) {
            const std::size_t k0 =
                static_cast<std::size_t>(nrow) *
                static_cast<std::size_t>(j);
            const double xa = net.coordinates[3 * k0];
            const double ya = net.coordinates[3 * k0 + 1];
            const double dr = std::hypot(xa - c1.x(), ya - c1.y());
            const double za = c1.z() + dr * zda;
            const std::size_t kb = k0 + static_cast<std::size_t>(nrm1);
            const double xb = net.coordinates[3 * kb];
            const double yb = net.coordinates[3 * kb + 1];
            const double ds = std::hypot(xb - c4.x(), yb - c4.y());
            const double zb = c4.z() + ds * zdb;
            const double dd = std::hypot(xb - xa, yb - ya);
            const double zd = (dd != 0.0) ? (zb - za) / dd : 0.0;
            for (int i = 0; i < nrow; ++i) {
                const std::size_t index =
                    k0 + static_cast<std::size_t>(i);
                const double xc = net.coordinates[3 * index];
                const double yc = net.coordinates[3 * index + 1];
                const double dt = std::hypot(xc - xa, yc - ya);
                net.coordinates[3 * index + 2] = za + dt * zd;
            }
        }
    }

    networks.push_back(std::move(net));
}

void AeroModel::extractEllipticalSection(const Section& sec) {
    // Port of ellpt(): a body of revolution with elliptical cross-sections.
    // The profile cards hold (xs, rY, rZ) per station; the angular positions
    // th (degrees) map to
    //   y = rY*rZ*cos(th)/den, z = z0 + rY*rZ*sin(th)/den,
    //   den = sqrt(rZ^2*cos(th)^2 + rY^2*sin(th)^2).
    struct Block {
        bool started{false};
        int nopt{0};
        std::string name;
        int nm{0};
        int nn{0};
        std::vector<double> xs, ry, rz, z, th;
    };

    int sectionKn = 0;
    int sectionKt = 0;
    double sectionCpnorm = 2.0;
    bool pendingKn = false;
    bool pendingKt = false;

    auto groupValues = [&](std::size_t index) {
        std::vector<Value> values = sec.entries[index].values;
        for (std::size_t j = index + 1;
             j < sec.entries.size() &&
             sec.entries[j].kind == EntryKind::FreeData;
             ++j) {
            values.insert(values.end(),
                          sec.entries[j].values.begin(),
                          sec.entries[j].values.end());
        }
        return values;
    };
    auto numbers = [](const std::vector<Value>& values) {
        std::vector<double> out;
        for (const Value& v : values)
            if (v.numeric)
                out.push_back(v.number);
        return out;
    };

    Block current;
    auto finishBlock = [&]() {
        if (!current.started)
            return;
        const bool complete =
            !current.name.empty() && current.nm > 0 && current.nn > 0 &&
            static_cast<int>(current.xs.size()) >= current.nm &&
            static_cast<int>(current.ry.size()) >= current.nm &&
            static_cast<int>(current.rz.size()) >= current.nm &&
            static_cast<int>(current.th.size()) >= current.nn;
        if (!complete) {
            BLOG(warning) << "elliptical network '" << current.name
                          << "' is incomplete";
            current = Block{};
            return;
        }

        Network net;
        net.netname = current.name;
        net.nm = current.nm;
        net.nn = current.nn;
        net.kn = sectionKn;
        net.kt = sectionKt;
        net.cpnorm = sectionCpnorm;
        if (static_cast<int>(current.z.size()) < current.nm)
            current.z.assign(static_cast<std::size_t>(current.nm), 0.0);

        constexpr double pi = 3.14159265358979323846;
        const double deg = pi / 180.0;
        net.coordinates.reserve(
            static_cast<std::size_t>(3) * current.nm * current.nn);
        for (int j = 0; j < current.nn; ++j) {
            const double st =
                std::sin(current.th[static_cast<std::size_t>(j)] * deg);
            const double ct =
                std::cos(current.th[static_cast<std::size_t>(j)] * deg);
            for (int i = 0; i < current.nm; ++i) {
                const std::size_t si = static_cast<std::size_t>(i);
                const double ry = current.ry[si];
                const double rz = current.rz[si];
                const double den =
                    std::sqrt(rz * rz * ct * ct + ry * ry * st * st);
                net.coordinates.push_back(current.xs[si]);
                if (den > 0.0) {
                    net.coordinates.push_back(ry * rz * ct / den);
                    net.coordinates.push_back(current.z[si] +
                                              ry * rz * st / den);
                } else {
                    net.coordinates.push_back(0.0);
                    net.coordinates.push_back(current.z[si]);
                }
            }
        }
        networks.push_back(std::move(net));
        current = Block{};
    };

    for (std::size_t ei = 0; ei < sec.entries.size(); ++ei) {
        const Entry& e = sec.entries[ei];
        if (e.kind == EntryKind::Variable && !e.names.empty()) {
            const std::string& first = e.names.front();
            if (first == "kn") {
                pendingKn = true;
            } else if (first == "kt") {
                pendingKt = true;
            } else if (first == "nopt") {
                finishBlock();
                const std::vector<Value> values = groupValues(ei);
                current.started = true;
                for (const Value& v : values) {
                    if (v.numeric) {
                        current.nopt = static_cast<int>(
                            std::lround(v.number));
                        break;
                    }
                }
                for (const Value& v : values) {
                    if (!v.numeric) {
                        current.name = v.text;
                        break;
                    }
                }
            } else if (first == "nm") {
                const std::vector<double> values = numbers(groupValues(ei));
                if (!values.empty())
                    current.nm = static_cast<int>(
                        std::lround(values.front()));
            } else if (first == "xs") {
                const std::vector<double> values = numbers(groupValues(ei));
                current.xs.clear();
                current.ry.clear();
                current.rz.clear();
                for (std::size_t k = 0; k + 2 < values.size(); k += 3) {
                    current.xs.push_back(values[k]);
                    current.ry.push_back(values[k + 1]);
                    current.rz.push_back(values[k + 2]);
                }
            } else if (first == "z") {
                current.z = numbers(groupValues(ei));
            } else if (first == "nn") {
                const std::vector<double> values = numbers(groupValues(ei));
                if (!values.empty())
                    current.nn = static_cast<int>(
                        std::lround(values.front()));
            } else if (first == "th") {
                current.th = numbers(groupValues(ei));
            }
            continue;
        }
        if (e.kind != EntryKind::FreeData)
            continue;

        if (pendingKn) {
            const std::vector<double> values = numbers(e.values);
            if (!values.empty()) {
                sectionKn = static_cast<int>(
                    std::lround(values.front()));
                if (values.size() >= 2)
                    sectionCpnorm = values.back();
            }
            pendingKn = false;
        } else if (pendingKt) {
            const std::vector<double> values = numbers(e.values);
            if (!values.empty())
                sectionKt = static_cast<int>(
                    std::lround(values.front()));
            pendingKt = false;
        }
    }
    finishBlock();
}

void AeroModel::extractTrailingSection(const Section& sec) {
    // Port of trwake(): generate a two-row wake network attached to one
    // side of an existing network.  The "=inat insd xwake twake netname"
    // group provides one row per generated network.
    int sectionKn = 0;
    int sectionKt = 0;
    double sectionCpnorm = 2.0;
    bool pendingKn = false;
    bool pendingKt = false;
    std::vector<Value> trailingRows;

    auto groupValues = [&](std::size_t index) {
        std::vector<Value> values = sec.entries[index].values;
        for (std::size_t j = index + 1;
             j < sec.entries.size() &&
             sec.entries[j].kind == EntryKind::FreeData;
             ++j) {
            values.insert(values.end(),
                          sec.entries[j].values.begin(),
                          sec.entries[j].values.end());
        }
        return values;
    };
    auto numbers = [](const std::vector<Value>& values) {
        std::vector<double> out;
        for (const Value& v : values)
            if (v.numeric)
                out.push_back(v.number);
        return out;
    };

    for (std::size_t ei = 0; ei < sec.entries.size(); ++ei) {
        const Entry& e = sec.entries[ei];
        if (e.kind == EntryKind::Variable && !e.names.empty()) {
            const std::string& first = e.names.front();
            if (first == "kn")
                pendingKn = true;
            else if (first == "kt")
                pendingKt = true;
            else if (first == "inat")
                trailingRows = groupValues(ei);
            continue;
        }
        if (e.kind != EntryKind::FreeData)
            continue;
        if (pendingKn) {
            const std::vector<double> values = numbers(e.values);
            if (!values.empty()) {
                sectionKn = static_cast<int>(
                    std::lround(values.front()));
                if (values.size() >= 2)
                    sectionCpnorm = values.back();
            }
            pendingKn = false;
        } else if (pendingKt) {
            const std::vector<double> values = numbers(e.values);
            if (!values.empty())
                sectionKt = static_cast<int>(
                    std::lround(values.front()));
            pendingKt = false;
        }
    }

    if (trailingRows.empty())
        return;

    constexpr double pi = 3.14159265358979323846;
    const double aarg = alpc * pi / 180.0;
    const double barg = -betc * pi / 180.0;
    const Eigen::Vector3d freestream(
        std::cos(aarg) * std::cos(barg),
        std::cos(aarg) * std::sin(barg), std::sin(aarg));

    for (std::size_t i = 0; i + 4 < trailingRows.size(); i += 5) {
        const Value& attachValue = trailingRows[i];
        const Value& sideValue = trailingRows[i + 1];
        const Value& xwakeValue = trailingRows[i + 2];
        const Value& twakeValue = trailingRows[i + 3];
        const Value& nameValue = trailingRows[i + 4];
        if (attachValue.numeric || !sideValue.numeric ||
            !xwakeValue.numeric || !twakeValue.numeric ||
            nameValue.numeric) {
            BLOG(warning) << "$trailing: malformed =inat row";
            break;
        }

        const std::string attachName = lowerCopy(attachValue.text);
        int attached = -1;
        for (std::size_t k = 0; k < networks.size(); ++k)
            if (lowerCopy(networks[k].netname) == attachName) {
                attached = static_cast<int>(k);
                break;
            }
        if (attached < 0) {
            BLOG(warning) << "$trailing: network '" << attachValue.text
                          << "' not found";
            continue;
        }

        const Network& parent =
            networks[static_cast<std::size_t>(attached)];
        const int nrow = parent.nm;
        const int ncol = parent.nn;
        const int insd = static_cast<int>(std::lround(sideValue.number));
        int nte = 0;
        int sideIndex = 0;
        bool sideIsRow = false;
        if (insd == 1 || insd == 3) {
            nte = ncol;
            sideIndex = (insd == 1) ? 0 : nrow - 1;
            sideIsRow = true;
        } else if (insd == 2 || insd == 4) {
            nte = nrow;
            sideIndex = (insd == 2) ? ncol - 1 : 0;
            sideIsRow = false;
        } else {
            BLOG(warning) << "$trailing: invalid insd=" << insd;
            continue;
        }
        if (nte <= 0)
            continue;

        const Eigen::Vector3d gen =
            (twakeValue.number == 1.0)
                ? freestream
                : Eigen::Vector3d(1.0, 0.0, 0.0);
        if (!(std::abs(gen.x()) > 1.0e-14)) {
            BLOG(warning) << "$trailing: zero streamwise generator";
            continue;
        }

        Network wake;
        wake.netname = nameValue.text;
        wake.nm = 2;
        wake.nn = nte;
        wake.kn = sectionKn;
        wake.kt = sectionKt;
        wake.cpnorm = sectionCpnorm;
        wake.coordinates.reserve(
            static_cast<std::size_t>(6) * static_cast<std::size_t>(nte));

        bool complete = true;
        for (int j = 0; j < nte; ++j) {
            const int pIndex = sideIsRow
                                   ? j * nrow + sideIndex
                                   : sideIndex * nrow + j;
            const std::size_t base =
                static_cast<std::size_t>(3 * pIndex);
            if (base + 2 >= parent.coordinates.size()) {
                complete = false;
                break;
            }
            const Eigen::Vector3d p(parent.coordinates[base],
                                    parent.coordinates[base + 1],
                                    parent.coordinates[base + 2]);
            const double alf = (xwakeValue.number - p.x()) / gen.x();
            const Eigen::Vector3d q = p + alf * gen;
            wake.coordinates.push_back(p.x());
            wake.coordinates.push_back(p.y());
            wake.coordinates.push_back(p.z());
            wake.coordinates.push_back(q.x());
            wake.coordinates.push_back(q.y());
            wake.coordinates.push_back(q.z());
        }
        if (!complete) {
            BLOG(warning) << "$trailing: network '" << parent.netname
                          << "' has too few coordinates";
            continue;
        }
        networks.push_back(std::move(wake));
    }
}

void AeroModel::extractMaterials() {
    materials.assign(1, MaterialProperties{});
    surfaceMaterials.clear();

    const Section* section = findSection("material");
    if (section == nullptr)
        return;

    // A "=name" declaration owns the ordinary data cards that follow it
    // until the next "=..." group.  The parser keeps those cards as separate
    // FreeData entries, so collect them here.
    auto groupValues = [&](std::size_t index) {
        std::vector<Value> values = section->entries[index].values;
        for (std::size_t j = index + 1;
             j < section->entries.size() &&
             section->entries[j].kind == EntryKind::FreeData;
             ++j) {
            values.insert(values.end(),
                          section->entries[j].values.begin(),
                          section->entries[j].values.end());
        }
        return values;
    };

    int materialCount = 0;
    for (std::size_t i = 0; i < section->entries.size(); ++i) {
        const Entry& e = section->entries[i];
        if (e.kind != EntryKind::Variable || e.names.empty())
            continue;
        if (e.names.front() != "nmat")
            continue;
        const std::vector<Value> values = groupValues(i);
        for (const Value& v : values) {
            if (v.numeric) {
                materialCount =
                    static_cast<int>(std::lround(v.number));
                break;
            }
        }
        break;
    }
    materialCount = std::min(std::max(materialCount, 0), 10);

    // Material definitions: "tratio pratio name" per material.
    for (std::size_t i = 0; i < section->entries.size(); ++i) {
        const Entry& e = section->entries[i];
        if (e.kind != EntryKind::Variable || e.names.empty())
            continue;
        if (e.names.front().find("tratio") == std::string::npos)
            continue;
        const std::vector<Value> values = groupValues(i);
        for (std::size_t k = 0; k + 2 < values.size(); ++k) {
            if (materialCount > 0 &&
                static_cast<int>(materials.size()) - 1 >= materialCount)
                break;
            if (!values[k].numeric || !values[k + 1].numeric ||
                values[k + 2].numeric)
                continue;
            MaterialProperties material;
            material.temperatureRatio = values[k].number;
            material.pressureRatio = values[k + 1].number;
            material.name = values[k + 2].text;
            if (material.temperatureRatio == 0.0)
                material.temperatureRatio = 1.0;
            if (material.pressureRatio == 0.0)
                material.pressureRatio = 1.0;
            materials.push_back(std::move(material));
            k += 2;
        }
    }

    // Same clamps as matprp: air properties at very low Mach, and pressure
    // ratios not below the critical value.
    constexpr double gamma = 1.4;
    if (amach < 0.01) {
        for (std::size_t i = 1; i < materials.size(); ++i) {
            materials[i].temperatureRatio = 1.0;
            materials[i].pressureRatio = 1.0;
        }
    } else {
        const double tpexp = gamma / (1.0 - gamma);
        double tpcrit =
            std::pow(1.0 + 0.5 * (gamma - 1.0) * amach * amach, tpexp) +
            1.0e-12;
        tpcrit = std::min(tpcrit, 1.0);
        for (std::size_t i = 1; i < materials.size(); ++i)
            if (materials[i].pressureRatio < tpcrit)
                materials[i].pressureRatio = tpcrit;
    }

    // Network/surface -> material assignments.
    for (std::size_t ei = 0; ei < section->entries.size(); ++ei) {
        const Entry& e = section->entries[ei];
        if (e.kind != EntryKind::Variable || e.names.empty())
            continue;
        if (e.names.front() != "netname")
            continue;
        const std::vector<Value> values = groupValues(ei);
        for (std::size_t i = 0; i + 2 < values.size(); i += 3) {
            SurfaceMaterial assignment;
            assignment.network = values[i].text;
            const std::string surface = lowerCopy(values[i + 1].text);
            if (surface.find("upp") != std::string::npos)
                assignment.surface = 1;
            else if (surface.find("low") != std::string::npos)
                assignment.surface = -1;
            else
                continue;

            int material = 0;
            if (values[i + 2].numeric) {
                material = static_cast<int>(
                    std::lround(values[i + 2].number));
            } else {
                const std::string wanted = lowerCopy(values[i + 2].text);
                for (std::size_t k = 0; k < materials.size(); ++k)
                    if (lowerCopy(materials[k].name) == wanted) {
                        material = static_cast<int>(k);
                        break;
                    }
            }
            if (material < 0 ||
                material >= static_cast<int>(materials.size()))
                material = 0;
            assignment.material = material;
            surfaceMaterials.push_back(std::move(assignment));
        }
    }
}

int AeroModel::materialIndex(const std::string& netname, int surface) const {
    const std::string wanted = lowerCopy(netname);
    int material = 0;
    for (const SurfaceMaterial& assignment : surfaceMaterials) {
        if (assignment.surface != surface)
            continue;
        if (lowerCopy(assignment.network) == wanted)
            material = assignment.material; // later entries override
    }
    return material;
}

bool AeroModel::isSolved() const {
    return solved_;
}

bool AeroModel::solve() {
    // Geometric precursor only: build the nine canonical points and the eight
    // triangular subpanels for every quad cell.  Full higher-order AIC solve
    // lives in the Fortran A502 code and is not ported here yet.
    //
    // Panel order is the A502 order (geomc): columns outer, rows inner:
    //   panel index = i + (nm-1)*j.
    panels.clear();
    controlPoints.clear();
    totalSurfaceArea = 0.0;
    totalPanels = 0;
    std::size_t skipped = 0;
    solved_ = false;
    flowSolved_ = false;

    flow = Compressibility::fromAngles(amach, alpc, betc);

    for (std::size_t netIdx = 0; netIdx < networks.size(); ++netIdx) {
        const Network& n = networks[netIdx];
        const int nm = n.nm;
        const int nn = n.nn;
        if (nm < 2 || nn < 2)
            continue;

        const std::size_t nPoints =
            static_cast<std::size_t>(nm) * static_cast<std::size_t>(nn);
        if (n.coordinates.size() < nPoints * 3) {
            BLOG(warning) << "solve(): network '" << n.netname
                          << "' has too few coordinates ("
                          << n.coordinates.size() / 3 << " of " << nPoints
                          << ")";
            continue;
        }

        for (int j = 0; j + 1 < nn; ++j) {
            for (int i = 0; i + 1 < nm; ++i) {
                PanelGeometry geom =
                    buildPanelGeometry(n.coordinates, nm, nn, i, j);
                if (geom.degenerate) {
                    ++skipped;
                    continue;
                }
                applyCompressibleFrames(geom, flow);

                Panel pan;
                pan.geom = std::move(geom);
                pan.network = static_cast<int>(netIdx);
                pan.panelI = i;
                pan.panelJ = j;
                // Synthetic Cp so ParaView has a visible scalar (placeholder).
                pan.cp = pan.geom.normal.x();
                totalSurfaceArea += pan.geom.area;
                panels.push_back(std::move(pan));
            }
        }

        std::vector<ControlPoint> cps = buildNetworkControlPoints(n);
        for (auto& cp : cps) {
            cp.network = static_cast<int>(netIdx);
            controlPoints.push_back(std::move(cp));
        }
    }

    totalPanels = panels.size();
    solved_ = !panels.empty();
    if (solved_) {
        BLOG(info) << "solve(): " << totalPanels << " panels, "
                   << skipped << " degenerate skipped, total area="
                   << totalSurfaceArea << ", control points="
                   << controlPoints.size();
    } else {
        BLOG(warning) << "solve(): no usable panels from networks";
    }
    return solved_;
}

bool AeroModel::isFlowSolved() const {
    return flowSolved_;
}

bool AeroModel::solveFlow() {
    flowSolved_ = false;
    if (!solved_ || networks.empty())
        return false;

    const std::vector<Network>& nets = networks;
    const GlobalSingularityLayout layout =
        buildGlobalSingularityLayout(nets);
    if (layout.network.size() != nets.size()) {
        BLOG(warning) << "solveFlow(): singularity layout failed";
        return false;
    }

    // Source-only networks have one unknown per panel and use panel-centre
    // rows.  Doublet-only networks use the staggered parameter grid, whose
    // cardinality equals the number of control points.  Composite networks
    // are still under-determined.
    std::vector<bool> includeEdgeAndCorner(nets.size(), false);
    for (std::size_t k = 0; k < nets.size(); ++k) {
        const NetworkSingularityLayout& netLayout = layout.network[k];
        const bool sourceOnly = netLayout.hasSource && !netLayout.hasDoublet;
        const bool doubletOnly = netLayout.hasDoublet && !netLayout.hasSource;
        if (!sourceOnly && !doubletOnly) {
            BLOG(warning) << "solveFlow(): network '" << nets[k].netname
                          << "' is composite";
            return false;
        }
        if (nets[k].nm < 2 || nets[k].nn < 2)
            return false;
        includeEdgeAndCorner[k] = doubletOnly;
    }

    const GlobalBoundarySystem system = buildGlobalNormalFluxSystem(
        nets, flow, layout, 1.0, 1.0, 1.0e-10, includeEdgeAndCorner);
    if (!system.valid) {
        BLOG(warning) << "solveFlow(): boundary system assembly failed";
        return false;
    }

    const GlobalSystemSolution solution =
        solveGlobalBoundarySystem(system);
    if (!solution.valid) {
        const auto diagnostics = analyzeGlobalBoundarySystem(system);
        BLOG(warning) << "solveFlow(): global solve failed (rows="
                      << system.nRows << ", cols=" << system.nColumns
                      << ", rank=" << diagnostics.rank << ")";
        return false;
    }

    const GlobalParameterValues values =
        expandGlobalSolution(solution, layout);
    if (!values.valid) {
        BLOG(warning) << "solveFlow(): solution expansion failed";
        return false;
    }

    const SurfaceFlowResult surface = evaluateSurfaceFlow(
        nets, flow, layout, values, true);
    if (!surface.valid) {
        BLOG(warning) << "solveFlow(): surface flow evaluation failed";
        return false;
    }

    // evaluateSurfaceFlow walks the same all-control-point order as
    // AeroModel::controlPoints, so record each network's panel-centre row.
    std::vector<std::vector<int>> centerIndex(nets.size());
    for (std::size_t k = 0; k < nets.size(); ++k)
        centerIndex[k].assign(
            static_cast<std::size_t>((nets[k].nm - 1) *
                                     (nets[k].nn - 1)),
            -1);
    std::size_t surfaceIndex = 0;
    for (const ControlPoint& cp : controlPoints) {
        if (cp.kind == ControlPointKind::PanelCenter &&
            cp.network >= 0 &&
            cp.network < static_cast<int>(centerIndex.size())) {
            std::vector<int>& index =
                centerIndex[static_cast<std::size_t>(cp.network)];
            if (cp.panelIndex >= 0 &&
                cp.panelIndex < static_cast<int>(index.size()))
                index[static_cast<std::size_t>(cp.panelIndex)] =
                    static_cast<int>(surfaceIndex);
        }
        ++surfaceIndex;
    }
    if (surfaceIndex != surface.cp.size()) {
        BLOG(warning) << "solveFlow(): control-point index mismatch ("
                      << surfaceIndex << " vs " << surface.cp.size()
                      << ")";
        return false;
    }

    for (Panel& panel : panels) {
        if (panel.network < 0 ||
            panel.network >= static_cast<int>(nets.size()))
            continue;
        const Network& net = nets[static_cast<std::size_t>(panel.network)];
        const std::vector<int>& index =
            centerIndex[static_cast<std::size_t>(panel.network)];
        const int panelIndex = panel.panelI + (net.nm - 1) * panel.panelJ;
        if (panelIndex < 0 || panelIndex >= static_cast<int>(index.size()))
            continue;
        const int si = index[static_cast<std::size_t>(panelIndex)];
        if (si >= 0)
            panel.cp = surface.cp[static_cast<std::size_t>(si)];
    }

    flowSolved_ = true;
    BLOG(info) << "solveFlow(): " << nets.size() << " network(s), unknowns="
               << solution.unknowns.size()
               << ", residual=" << solution.residualNorm;
    return true;
}

} // namespace a502
