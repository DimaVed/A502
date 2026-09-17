// AeroModel.cpp
// -------------
// Implementation of the A502 input deck parser / model (see AeroModel.h).

#include "AeroModel.h"

#include <algorithm>
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
}

void AeroModel::parse(const std::string& inpPath) {
    clear();

    std::ifstream in(inpPath);
    if (!in.is_open())
        throw std::runtime_error("AeroModel::parse: cannot open input file \"" +
                                 inpPath + "\"");

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
}
void AeroModel::parse(const std::string& inpPath, const std::string& outPath) {
    parse(inpPath);
    writeOut(outPath);
}

void AeroModel::writeOut(const std::string& outPath) const {
    std::ofstream out(outPath);
    if (!out.is_open())
        throw std::runtime_error("AeroModel::writeOut: cannot create file \"" +
                                 outPath + "\"");

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
        out << "  net '" << n.netname << "'  kn=" << n.kn << "  kt=" << n.kt
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
}

void AeroModel::PrintParaview(const std::string& vtkPath) const {
#ifdef WITH_VTK
    if (networks.empty())
        throw std::runtime_error(
            "AeroModel::PrintParaview: model contains no grid networks");

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
    if (writer->Write() == 0)
        throw std::runtime_error(
            "AeroModel::PrintParaview: VTK writer failed for \"" + vtkPath +
            "\"");
#else
    (void)vtkPath;
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

void AeroModel::extractNetworks() {
    for (const auto& sec : sections) {
        if (sec.keywordNorm != "points")
            continue;

        // Flatten every data token of the section into a single sequence,
        // keeping the physical order in which the cards appeared.
        std::vector<Value> flat;
        for (const auto& e : sec.entries)
            flat.insert(flat.end(), e.values.begin(), e.values.end());

        // Network names (non-numeric tokens such as "ab1", "wingta", "kt5")
        // delimit the "nm nn netname" triplets.  Every such occurrence starts
        // the description of one network.
        std::vector<std::size_t> nameIdx;
        for (std::size_t i = 0; i < flat.size(); ++i)
            if (!flat[i].numeric)
                nameIdx.push_back(i);

        if (nameIdx.empty())
            continue;

        // Section-level kn / kt / cpnorm taken from the leading cards that
        // precede the first "nm nn netname" triplet.
        int    kn = 0, kt = 0;
        double cpnorm = 0.0;
        {
            std::size_t first = nameIdx.front();
            std::vector<double> lead;
            // The "nm nn" pair (index first-1, first-2) belongs to the first
            // network, not to the section header, so exclude it here.
            for (std::size_t i = 0; i + 2 < first; ++i)
                if (flat[i].numeric)
                    lead.push_back(flat[i].number);
            if (!lead.empty()) {
                kn = static_cast<int>(std::lround(lead.front()));
                if (lead.size() >= 3) {
                    cpnorm = lead[1];
                    kt = static_cast<int>(std::lround(lead.back()));
                } else {
                    kt = static_cast<int>(std::lround(lead.back()));
                }
            }
        }

        for (std::size_t k = 0; k < nameIdx.size(); ++k) {
            std::size_t a = nameIdx[k];
            if (a < 2 || !flat[a - 1].numeric || !flat[a - 2].numeric) {
                // A name without its nm/nn pair before it (legacy layout);
                // skip to keep the model well formed.
                continue;
            }
            Network n;
            n.netname = flat[a].text;
            n.nm = static_cast<int>(std::lround(flat[a - 2].number));
            n.nn = static_cast<int>(std::lround(flat[a - 1].number));
            n.kn = kn;
            n.kt = kt;
            n.cpnorm = cpnorm;

            // Coordinates of this network run from just after its netname up
            // to (but excluding) the nm/nn pair of the following network.
            std::size_t b = (k + 1 < nameIdx.size()) ? nameIdx[k + 1] : flat.size();
            std::size_t end = (k + 1 < nameIdx.size() && b >= 2) ? b - 2 : b;
            for (std::size_t i = a + 1; i < end; ++i)
                if (flat[i].numeric)
                    n.coordinates.push_back(flat[i].number);

            networks.push_back(std::move(n));
        }
    }
}

} // namespace a502
