// Solver
// ------
// Parse an A502 / PanAir input deck and write a text model report.
// Full panel AIC solve is not ported yet; this entry point exercises the
// C++ deck reader and geometry precursor.
//
// usage:
//   Solver [options] <input.inp> [report.out]
//   Solver --help
//   Solver --version

#include <filesystem>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "AeroModel.h"
#include "Logger.h"
#include "Version.h"

namespace fs = std::filesystem;
namespace po = boost::program_options;

static fs::path defaultReportPath(const fs::path& inp) {
    fs::path out = inp;
    out.replace_extension(".report.txt");
    return out;
}

int main(int argc, char** argv) {
    try {
        po::options_description desc("Solver — A502 input parse / geometry report");
        desc.add_options()
            ("help,h", "show this help message and exit")
            ("version,V", "print version and exit")
            ("input,i", po::value<std::string>(),
             "A502 / PanAir input deck (.inp)")
            ("output,o", po::value<std::string>(),
             "text report path (default: <input>.report.txt)");

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
                      << "  Solver samples/swb.inp\n"
                      << "  Solver -i samples/swb.inp -o swb.report.txt\n";
            return 0;
        }
        if (vm.count("version")) {
            std::cout << "Solver " << a502::kVersion
                      << " (" << a502::kVersionDate << ")\n"
                      << a502::kProductName << "\n";
            return 0;
        }
        if (!vm.count("input")) {
            std::cerr << "error: input file is required\n\n" << desc;
            return 2;
        }

        InitBoostLogFilter("logs/solver", boost::log::trivial::info);

        const fs::path inpPath = vm["input"].as<std::string>();
        const fs::path reportPath =
            vm.count("output") ? fs::path(vm["output"].as<std::string>())
                               : defaultReportPath(inpPath);

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
