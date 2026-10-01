// Visualizer
// ----------
// Model / geometry viewer: reads an A502 / PanAir input deck (.inp) and
// writes a VTK UnstructuredGrid (.vtu) of the panel networks for ParaView.
// For solution fields (Cp, velocity) use ResultVisualizer on a *.out file.
//
// usage:
//   Visualizer [options] <input.inp> [output.vtu]
//   Visualizer --help
//   Visualizer --version

#include <filesystem>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "AeroModel.h"
#include "Logger.h"
#include "Version.h"

namespace fs = std::filesystem;
namespace po = boost::program_options;

static fs::path defaultVtuPath(const fs::path& inp) {
    fs::path out = inp;
    out.replace_extension(".vtu");
    return out;
}

int main(int argc, char** argv) {
    try {
        po::options_description desc(
            "Visualizer — A502 model/geometry (.inp) to ParaView VTU");
        desc.add_options()
            ("help,h", "show this help message and exit")
            ("version,V", "print version and exit")
            ("input,i", po::value<std::string>(),
             "A502 / PanAir input deck (.inp)")
            ("output,o", po::value<std::string>(),
             "output VTU path (default: <input>.vtu)");

        po::positional_options_description pos;
        pos.add("input", 1);
        pos.add("output", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv)
                      .options(desc)
                      .positional(pos)
                      .run(),
                  vm);
        po::notify(vm);

        if (vm.count("help")) {
            std::cout << desc << "\n"
                      << "This tool exports input-deck geometry only.\n"
                      << "For Cp / velocity from a solution use ResultVisualizer.\n\n"
                      << "Examples:\n"
                      << "  Visualizer samples/swb.inp\n"
                      << "  Visualizer -i samples/swb.inp -o swb.vtu\n";
            return 0;
        }
        if (vm.count("version")) {
            std::cout << "Visualizer " << a502::kVersion
                      << " (" << a502::kVersionDate << ")\n"
                      << a502::kProductName << "\n";
            return 0;
        }
        if (!vm.count("input")) {
            std::cerr << "error: input file is required\n\n" << desc;
            return 2;
        }

        InitBoostLogFilter("logs/visualizer", boost::log::trivial::info);

        const fs::path inpPath = vm["input"].as<std::string>();
        const fs::path vtuPath =
            vm.count("output") ? fs::path(vm["output"].as<std::string>())
                               : defaultVtuPath(inpPath);

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
    } catch (const po::error& e) {
        std::cerr << "error: " << e.what() << "\n"
                  << "Try --help for usage.\n";
        return 2;
    } catch (const std::exception& e) {
        try {
            BLOG(error) << e.what();
        } catch (...) {
            std::cerr << "error: " << e.what() << "\n";
        }
        return 1;
    }
}
