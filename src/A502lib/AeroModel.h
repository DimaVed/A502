#pragma once

// AeroModel.h
// -----------
// AeroModel is a first step of the ongoing rewrite of the PanAir / A502
// (Higher Order Panel Aerodynamics) program from Fortran to C++.
//
// The class holds a structured representation of an A502 *input deck*
// (a "namelist"-style text file, normally with the ".inp" extension) and is
// able to:
//   * parse such a file (see AeroModel::parse), preserving every section,
//     variable declaration, "*" marker and raw data card it contains, and
//   * write a human readable report of the parsed model to an output file
//     (see AeroModel::writeOut).
//
// The parser is deliberately generic: it does not hard-code every possible
// section name.  Everything found in the input is stored in an ordered list
// of AeroModel::Section objects.  In addition a small set of well known
// scalar parameters (title, mach number, angles of attack, reference data,
// printout options, ...) and the simple grid networks ($points sections) are
// extracted into typed members as a convenience.
//
// Only the C++ standard library and Eigen (for geometry vectors) are used; no
// Boost / other external dependency is required by this translation unit.
#define EIGEN_NO_DEBUG           // Eigen full speed for the solver prototype

#include <ostream>
#include <string>
#include <vector>
#include <cstddef>

#include <Eigen/Dense>

namespace a502 {

// A single token read from a data card.  A token is either a number (that is
// how coordinates / coefficients are encoded) or a text word (network names,
// material names and so on).
struct Value {
    std::string text;   // original token verbatim
    bool        numeric{false};
    double      number{0.0};
};

enum class EntryKind {
    Variable,   // a "=name1 name2 ..." declaration followed by its values
    Star,       // a "*marker" card (e.g. *pressure, *cut) followed by values
    FreeData    // a plain data card not attached to any named group
};

// One assignment group stored in a Section.
struct Entry {
    EntryKind                kind{EntryKind::FreeData};
    std::string              text;   // raw trimmed card text (or the "=..." part)
    std::vector<std::string> names;  // Variable: the declared names
    std::vector<Value>       values; // the tokens that belong to this group
};

// One input section: everything between two '$' cards.
struct Section {
    std::string keyword;      // first word after '$' (kept as written)
    std::string keywordNorm;  // lower-cased, non-alphanumeric stripped
    std::string description;  // the rest of the '$' line
    std::vector<Entry> entries;
};

// A simple grid network reconstructed from a $points section.
// PanAir generates the full grid from a set of spanwise "columns"; we keep
// the metadata plus the flat sequence of raw coordinate values.
struct Network {
    std::string          netname;
    int                  nm{0};       // points per column
    int                  nn{0};       // number of columns
    int                  kn{0};       // boundary condition / network type
    int                  kt{0};       // panel type
    double               cpnorm{0.0}; // normal velocity boundary condition
    std::vector<double>  coordinates; // x,y,z triplets as read
};

// Reference quantities used for forces & moments accumulation.
struct ReferenceData {
    double xref{0.0}, yref{0.0}, zref{0.0};
    int    nref{1},   nprcof{1};
    double sref{0.0}, bref{0.0}, cref{0.0}, dref{0.0};
};

// One panel resulting from solve().  It carries the geometric data a
// higher-order panel method needs (centroid, oriented unit normal, area) plus
// a placeholder surface pressure coefficient (cp) so the result can be shaded
// in a visualizer.  Vectors are stored as Eigen::Vector3d.
struct Panel {
    Eigen::Vector3d center;  // panel centroid
    Eigen::Vector3d normal;  // oriented unit normal
    double          area{0.0};
    double          cp{0.0};
    int             network{0};  // owning $points network index
};

class AeroModel {
public:
    // Reads and parses the given A502 input deck.
    // Throws std::runtime_error if the file can not be opened.
    void parse(const std::string& inpPath);

    // Convenience: parse "inpPath" and immediately write the resulting model
    // report into "outPath" (shortcut for parse() + writeOut()).
    void parse(const std::string& inpPath, const std::string& outPath);

    // Writes a human readable report of the parsed model to "outPath".
    // Throws std::runtime_error if the file can not be created.
    void writeOut(const std::string& outPath) const;

    // Short summary via BLOG(info): title, flow parameters, network table.
    void printSummary() const;

    // Writes the parsed model as a VTK unstructured grid (.vtu) so it can be
    // inspected in ParaView.  Every $points network is turned into a quad mesh
    // (netname recorded as per-cell label, network index as per-point index).
    // Requires the project to be built with VTK support (WITH_VTK).  Throws
    // std::runtime_error if VTK is unavailable or the file can not be written.
    void PrintParaview(const std::string& vtkPath) const;

    // Removes all previously parsed data.
    void clear();

    // ----- typed, well-known scalar parameters -----------------------------
    std::vector<std::string> titleLines;          // $title
    int    ndtchk{0};                             // $datacheck
    int    misym{0}, mjsym{0};                    // $symmetric flow / $sym
    double amach{0.0}, betc{0.0};                 // $mach number, $yaw angles
    int    nacase{1};                             // $cases
    double alpc{0.0};                             // $angles-of-attack
    std::vector<double> alpha;                    // $angles-of-attack
    std::vector<double> beta;                     // $yaw angles / $sideslip
    std::vector<double> printoutOptions;          // $printout options
    ReferenceData refs;                           // $reference(s)
    double eat{0.0};                              // $eat
    int    ivcorr{0};                             // $velocity correction
    double nflowv{0.0}, tpoff{0.0};               // $flowfield properties
    std::vector<std::vector<double>> offBodyPoints; // $flow-fields / $xyz ...

    // ----- generic, complete model ----------------------------------------
    std::vector<Section>  sections;   // every input section, in file order
    std::vector<Network>  networks;   // networks reconstructed from $points

    // ----- prototype solution ---------------------------------------------
    // solve() panelizes every $points network into a flat list of Panels in
    // i-major cell order (the same order as the quads emitted by
    // PrintParaview, so the Cp field lines up with the cells).  It is the
    // geometric precursor of the real higher-order panel solve and stamps a
    // synthetic Cp field so the result can be inspected in a visualizer.
    std::vector<Panel>  panels;                // all cells, in i-major order
    double              totalSurfaceArea{0.0}; // sum of panel areas
    std::size_t         totalPanels{0};        // panels.size() after solve()

    bool solve();                  // false if the model has no usable grid
    bool isSolved() const;         // true after a successful solve()

    // Finds the first section whose normalized keyword contains "kwFragment"
    // (sub-string, case-insensitive).  Returns nullptr when not found.
    const Section* findSection(const std::string& kwFragment) const;

private:
    void extractParameters();
    void extractNetworks();
    void writeHeader(std::ostream& out) const;
    bool solved_{false};
};

} // namespace a502
