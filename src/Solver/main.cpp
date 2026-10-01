// Solver
// ------
// Parse an A502 / PanAir input deck and write a text model report.
// Full panel AIC solve is not ported yet; this entry point exercises the
// C++ deck reader and geometry precursor.
//
// usage:
//   Solver <input.inp> [report.out]

#include <filesystem>
#include <string>

#include "AeroModel.h"
#include "Logger.h"

namespace fs = std::filesystem;

static fs::path defaultReportPath(const fs::path& inp) {
    fs::path out = inp;
    out.replace_extension(".report.txt");
    return out;
}

int main(int argc, char** argv) {
    InitBoostLogFilter("logs/solver", boost::log::trivial::info);

    if (argc < 2) {
        BLOG(error) << "usage: " << (argc > 0 ? argv[0] : "Solver")
                    << " <input.inp> [report.out]";
        return 2;
    }

    const fs::path inpPath = argv[1];
    const fs::path reportPath =
        (argc >= 3) ? fs::path(argv[2]) : defaultReportPath(inpPath);

    try {
        a502::AeroModel model;
        BLOG(info) << "Reading A502 input: " << inpPath.string();
        model.parse(inpPath.string());

        if (model.networks.empty()) {
            BLOG(warning) << "No $points networks found in input";
        }

        model.printSummary();

        if (!model.solve()) {
            BLOG(warning) << "Geometric solve produced no panels";
        }

        model.writeOut(reportPath.string());
        BLOG(info) << "Done";
        return 0;
    } catch (const std::exception& e) {
        BLOG(error) << e.what();
        return 1;
    }
}
