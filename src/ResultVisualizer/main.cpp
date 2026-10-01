// ResultVisualizer
// ----------------
// Reads a PanAir / A502 text solution (*.out) and writes a ParaView VTU with
// panel Cp / velocity fields (and optional off-body samples).
// Optionally paints Cp onto surface quads from a companion *.inp mesh.
//
// usage:
//   ResultVisualizer [options] <input.out> [output.vtu]
//   ResultVisualizer --help
//   ResultVisualizer --version

#include <filesystem>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "AeroModel.h"
#include "AeroResult.h"
#include "Logger.h"
#include "Version.h"

namespace fs = std::filesystem;
namespace po = boost::program_options;

static fs::path defaultVtuPath(const fs::path& outFile) {
    fs::path p = outFile;
    p.replace_extension("");
    return fs::path(p.string() + "_result.vtu");
}

int main(int argc, char** argv) {
    try {
        po::options_description desc(
            "ResultVisualizer — A502 .out solution to ParaView VTU");
        desc.add_options()
            ("help,h", "show this help message and exit")
            ("version,V", "print version and exit")
            ("input,i", po::value<std::string>(),
             "A502 / PanAir solution report (.out)")
            ("output,o", po::value<std::string>(),
             "output VTU path (default: <input>_result.vtu)")
            ("mesh,m", po::value<std::string>(),
             "optional .inp mesh to paint Cp on surface quads");

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
                      << "Examples:\n"
                      << "  ResultVisualizer samples/swb.out\n"
                      << "  ResultVisualizer -i samples/swb.out -o swb_result.vtu\n"
                      << "  ResultVisualizer -i samples/swb.out --mesh samples/swb.inp\n";
            return 0;
        }
        if (vm.count("version")) {
            std::cout << "ResultVisualizer " << a502::kVersion
                      << " (" << a502::kVersionDate << ")\n"
                      << a502::kProductName << "\n";
            return 0;
        }
        if (!vm.count("input")) {
            std::cerr << "error: input .out file is required\n\n" << desc;
            return 2;
        }

        InitBoostLogFilter("logs/result_visualizer",
                           boost::log::trivial::info);

        const fs::path outPath = vm["input"].as<std::string>();
        const fs::path vtuPath =
            vm.count("output") ? fs::path(vm["output"].as<std::string>())
                               : defaultVtuPath(outPath);

        a502::AeroResult result;
        BLOG(info) << "Reading A502 result: " << outPath.string();
        result.parse(outPath.string());
        result.printSummary();

        if (result.panels.empty() && result.offBody.empty()) {
            BLOG(error) << "No panel or off-body data found in result file";
            return 1;
        }

        if (vm.count("mesh")) {
            const fs::path meshPath = vm["mesh"].as<std::string>();
            BLOG(info) << "Loading companion mesh: " << meshPath.string();
            a502::AeroModel model;
            model.parse(meshPath.string());
            result.PrintParaview(vtuPath.string(), model);
        } else {
            BLOG(info) << "Exporting panel/off-body point cloud to ParaView";
            result.PrintParaview(vtuPath.string());
        }

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
