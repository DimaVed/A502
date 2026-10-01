// Visualizer
// ----------
// Reads an A502 / PanAir input deck (.inp), prints a short mesh summary,
// and writes a VTK UnstructuredGrid (.vtu) for ParaView.
//
// usage:
//   Visualizer <input.inp> [output.vtu]
//
// If the output path is omitted, "<input-stem>.vtu" is written next to the
// input file (or in the current working directory if the stem alone is used).

#include <filesystem>
#include <string>

#include "AeroModel.h"
#include "Logger.h"

namespace fs = std::filesystem;

static fs::path defaultVtuPath(const fs::path& inp) {
    fs::path out = inp;
    out.replace_extension(".vtu");
    return out;
}

int main(int argc, char** argv) {
    InitBoostLogFilter("logs/visualizer", boost::log::trivial::info);

    if (argc < 2) {
        BLOG(error) << "usage: " << (argc > 0 ? argv[0] : "Visualizer")
                    << " <input.inp> [output.vtu]";
        return 2;
    }

    const fs::path inpPath = argv[1];
    const fs::path vtuPath =
        (argc >= 3) ? fs::path(argv[2]) : defaultVtuPath(inpPath);

    try {
        a502::AeroModel model;
        BLOG(info) << "Reading A502 input: " << inpPath.string();
        model.parse(inpPath.string());

        if (model.networks.empty()) {
            BLOG(error) << "No $points networks found in input";
            return 1;
        }

        model.printSummary();

        BLOG(info) << "Exporting geometry to ParaView";
        model.PrintParaview(vtuPath.string());
        BLOG(info) << "Done";
        return 0;
    } catch (const std::exception& e) {
        BLOG(error) << e.what();
        return 1;
    }
}
