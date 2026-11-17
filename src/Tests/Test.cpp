// Tests for the C++ A502 port.
//
// The first block covers the $points network reader: it must follow the same
// card sequence as meshp/inbc in panair.f90 (kn card, kt card, "nm nn
// netname" cards, then exactly nm*nn coordinate points).

#include <catch2/catch_test_macros.hpp>

#include "AeroModel.h"

#include <filesystem>
#include <fstream>
#include <string>

namespace fs = std::filesystem;

namespace {

fs::path samplePath(const char* name) {
    return fs::path(A502_SOURCE_DIR) / "samples" / name;
}

struct ExpectedNetwork {
    const char* name;
    int         kt;
    int         nm;
    int         nn;
};

void checkNetworks(const fs::path& deck,
                   const std::vector<ExpectedNetwork>& expected) {
    a502::AeroModel model;
    model.parse(deck.string());

    INFO("deck: " << deck.string());
    REQUIRE(model.networks.size() == expected.size());
    for (std::size_t i = 0; i < expected.size(); ++i) {
        const auto& n = model.networks[i];
        INFO("network #" << i);
        CHECK(n.netname == expected[i].name);
        CHECK(n.kt == expected[i].kt);
        CHECK(n.nm == expected[i].nm);
        CHECK(n.nn == expected[i].nn);
        CHECK(n.coordinates.size() ==
              static_cast<std::size_t>(n.nm) * n.nn * 3);
    }
}

} // namespace

TEST_CASE("$points: swb.inp", "[parser][points]") {
    checkNetworks(samplePath("swb.inp"),
                  {{"winga", 1, 11, 3},
                   {"wingta", 1, 6, 2},
                   {"bodyl", 1, 11, 3},
                   {"bodyu", 1, 11, 3},
                   {"bodyb", 5, 5, 2},
                   {"awbw", 20, 4, 2}});
}

TEST_CASE("$points: ellip.inp (legacy card layout)", "[parser][points]") {
    checkNetworks(samplePath("ellip.inp"),
                  {{"wingup", 1, 21, 11},
                   {"winglo", 1, 21, 11},
                   {"wingtip", 1, 21, 2},
                   {"wake", 18, 11, 11}});
}

TEST_CASE("$points: vbc/vepd/vss small networks", "[parser][points]") {
    checkNetworks(samplePath("vbc.inp"),
                  {{"kt11", 11, 4, 4},
                   {"kt15", 15, 4, 4},
                   {"kt17", 17, 4, 4},
                   {"kt12a", 12, 4, 4},
                   {"kt12b", 12, 4, 4},
                   {"kt14", 14, 4, 4}});

    checkNetworks(samplePath("vepd.inp"),
                  {{"kt11", 11, 4, 4},
                   {"kt15", 15, 4, 4},
                   {"kt17", 17, 4, 4},
                   {"kt12a", 12, 4, 4},
                   {"kt12b", 12, 4, 4},
                   {"kt14", 14, 4, 4},
                   {"kt5", 5, 4, 4},
                   {"kt21", 21, 4, 4}});

    checkNetworks(samplePath("vss.inp"),
                  {{"kt11", 11, 4, 4},
                   {"kt15", 15, 4, 4},
                   {"kt17", 17, 4, 4},
                   {"kt12a", 12, 4, 4},
                   {"kt12b", 12, 4, 4},
                   {"kt13", 13, 4, 4},
                   {"kt14", 14, 4, 4},
                   {"kt5", 5, 4, 4},
                   {"kt21", 21, 4, 4}});
}

TEST_CASE("$points: ppbc/btac", "[parser][points]") {
    checkNetworks(samplePath("ppbc.inp"),
                  {{"kt5", 5, 4, 4}, {"kt21", 21, 4, 4}});

    checkNetworks(samplePath("btac.inp"),
                  {{"ab1", 1, 51, 7},
                   {"ab2", 1, 51, 7},
                   {"akl", 2, 11, 11},
                   {"a2", 18, 41, 11}});
}

TEST_CASE("$points: synthetic deck and coordinate order", "[parser][points]") {
    const fs::path deck = fs::temp_directory_path() / "a502_parser_test.inp";
    {
        std::ofstream out(deck);
        REQUIRE(out.is_open());
        out << "$title\n"
               "synthetic parser test\n"
               "second title line\n"
               "$points\n"
               "=kn\n"
               "1.\n"
               "=kt\n"
               "1.\n"
               "=nm       nn                                                          netname\n"
               "2.        2.                                                          plate\n"
               "=x y z\n"
               "0. 0. 0. 1. 0. 0.\n"
               "0. 1. 0. 1. 1. 0.\n"
               "$end\n";
    }

    a502::AeroModel model;
    model.parse(deck.string());
    REQUIRE(model.networks.size() == 1);
    const auto& n = model.networks.front();
    CHECK(n.netname == "plate");
    CHECK(n.nm == 2);
    CHECK(n.nn == 2);
    CHECK(n.kt == 1);
    REQUIRE(n.coordinates.size() == 12);
    // Column-major: (m=1,n=1), (m=2,n=1), (m=1,n=2), (m=2,n=2).
    CHECK(n.coordinates[0] == 0.0);
    CHECK(n.coordinates[1] == 0.0);
    CHECK(n.coordinates[2] == 0.0);
    CHECK(n.coordinates[3] == 1.0);
    CHECK(n.coordinates[4] == 0.0);
    CHECK(n.coordinates[5] == 0.0);
    CHECK(n.coordinates[6] == 0.0);
    CHECK(n.coordinates[7] == 1.0);
    CHECK(n.coordinates[8] == 0.0);
    CHECK(n.coordinates[9] == 1.0);
    CHECK(n.coordinates[10] == 1.0);
    CHECK(n.coordinates[11] == 0.0);

    std::error_code ec;
    fs::remove(deck, ec);
}

TEST_CASE("$circular is not a $points section", "[parser][points]") {
    // nac6.inp uses $circular; the $points reader must not invent networks.
    a502::AeroModel model;
    model.parse(samplePath("nac6.inp").string());
    CHECK(model.networks.empty());
}
