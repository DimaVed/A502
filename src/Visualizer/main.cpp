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
#include <iostream>
#include <string>

#include "AeroModel.h"

namespace fs = std::filesystem;

static fs::path defaultVtuPath(const fs::path& inp) {
    fs::path out = inp;
    out.replace_extension(".vtu");
    return out;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "usage: " << (argc > 0 ? argv[0] : "Visualizer")
                  << " <input.inp> [output.vtu]\n";
        return 2;
    }

    const fs::path inpPath = argv[1];
    const fs::path vtuPath =
        (argc >= 3) ? fs::path(argv[2]) : defaultVtuPath(inpPath);

    try {
        a502::AeroModel model;
        std::cout << "Reading A502 input: " << inpPath.string() << "\n";
        model.parse(inpPath.string());

        if (model.networks.empty()) {
            std::cerr << "error: no $points networks found in input\n";
            return 1;
        }

        model.printSummary(std::cout);

        std::cout << "Writing ParaView VTU: " << vtuPath.string() << "\n";
        model.PrintParaview(vtuPath.string());
        std::cout << "Done.\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
