// Tests for the C++ A502 port.
//
// The first block covers the $points network reader: it must follow the same
// card sequence as meshp/inbc in panair.f90 (kn card, kt card, "nm nn
// netname" cards, then exactly nm*nn coordinate points).

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "AeroModel.h"
#include "PanelGeometry.h"

#include <cmath>
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

TEST_CASE("panel geometry: unit square", "[geometry]") {
    // Column-major grid:
    //   (m=0,n=0)=(0,0,0)  (m=1,n=0)=(1,0,0)
    //   (m=0,n=1)=(0,1,0)  (m=1,n=1)=(1,1,0)
    const std::vector<double> coords = {
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 1.0, 1.0, 0.0};

    const a502::PanelGeometry g =
        a502::buildPanelGeometry(coords, 2, 2, 0, 0);

    REQUIRE_FALSE(g.degenerate);
    CHECK(g.cp[0].isApprox(Eigen::Vector3d(0, 0, 0)));
    CHECK(g.cp[1].isApprox(Eigen::Vector3d(0, 1, 0)));
    CHECK(g.cp[2].isApprox(Eigen::Vector3d(1, 1, 0)));
    CHECK(g.cp[3].isApprox(Eigen::Vector3d(1, 0, 0)));
    CHECK(g.cp[8].isApprox(Eigen::Vector3d(0.5, 0.5, 0)));
    CHECK(g.center.isApprox(g.cp[8]));
    // A502 normal = (+column) x (+row) = (0,1,0) x (1,0,0).
    CHECK(g.normal.isApprox(Eigen::Vector3d(0, 0, -1)));
    CHECK(g.area == Catch::Approx(1.0).epsilon(1e-14));

    double subArea = 0.0;
    for (const auto& s : g.subpanels) {
        CHECK_FALSE(s.degenerate);
        CHECK(s.area > 0.0);
        CHECK(std::abs(s.normal.norm() - 1.0) < 1e-14);
        CHECK(s.normal.dot(g.normal) > 1.0 - 1e-14);
        CHECK((s.toGlobal * s.toLocal -
               Eigen::Matrix3d::Identity()).norm() < 1e-14);
        CHECK((s.toLocal * (s.vertex[0] - s.origin)).norm() < 1e-14);
        CHECK(std::abs((s.toLocal * (s.vertex[1] - s.origin)).z()) < 1e-14);
        CHECK(std::abs((s.toLocal * (s.vertex[2] - s.origin)).z()) < 1e-14);
        subArea += s.area;
    }
    CHECK(subArea == Catch::Approx(g.area).epsilon(1e-14));
}

TEST_CASE("panel geometry: collapsed edge becomes a triangle", "[geometry]") {
    // P1 == P2 (network points (m=0,n=0) and (m=0,n=1)), so edge 1 is
    // collapsed and the panel is the triangle P1,P3,P4.
    const std::vector<double> coords = {
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 1.0, 0.0};

    const a502::PanelGeometry g =
        a502::buildPanelGeometry(coords, 2, 2, 0, 0);

    CHECK_FALSE(g.degenerate);
    CHECK(g.collapsedEdge == 1);
    CHECK(g.area == Catch::Approx(0.5).epsilon(1e-14));
    CHECK(g.normal.isApprox(Eigen::Vector3d(0, 0, -1)));
}

TEST_CASE("panel geometry: sample decks", "[geometry]") {
    SECTION("swb") {
        a502::AeroModel model;
        model.parse(samplePath("swb.inp").string());
        REQUIRE(model.solve());
        CHECK(model.totalPanels == 72);

        const auto& first = model.panels.front();
        CHECK(first.network == 0);
        CHECK(first.geom.i == 0);
        CHECK(first.geom.j == 0);
        // A502 panel order is column-major: first panel is i=0, j=0.
        const auto& net = model.networks[0];
        REQUIRE(net.coordinates.size() >= 6);
        CHECK(first.geom.cp[0].isApprox(Eigen::Vector3d(
            net.coordinates[0], net.coordinates[1], net.coordinates[2])));
        const std::size_t p2 = static_cast<std::size_t>(net.nm) * 3;
        CHECK(first.geom.cp[1].isApprox(Eigen::Vector3d(
            net.coordinates[p2], net.coordinates[p2 + 1],
            net.coordinates[p2 + 2])));

        for (const auto& pan : model.panels) {
            const auto& g = pan.geom;
            INFO("network " << pan.network << " panel " << g.i << "," << g.j);
            CHECK_FALSE(g.degenerate);
            CHECK(g.area > 0.0);
            CHECK(std::abs(g.normal.norm() - 1.0) < 1e-12);

            double subArea = 0.0;
            for (const auto& s : g.subpanels) {
                subArea += s.area;
                CHECK(std::abs(s.normal.norm() - 1.0) < 1e-10);
            }
            CHECK(subArea == Catch::Approx(g.area).epsilon(1e-10));

            // P5..P9 are the planar inner parallelogram; P9 is its center.
            const Eigen::Vector3d& n = g.normal;
            const double scale = std::max(1.0, g.cp[0].norm());
            for (int k = 4; k < 9; ++k)
                CHECK(std::abs((g.cp[k] - g.cp[4]).dot(n)) <
                      1e-9 * scale);
        }
    }

    SECTION("ellip and vss") {
        for (const char* name : {"ellip.inp", "vss.inp"}) {
            a502::AeroModel model;
            model.parse(samplePath(name).string());
            REQUIRE(model.solve());
            CHECK(model.totalPanels > 0);
            INFO("deck " << name);
            for (const auto& pan : model.panels) {
                double subArea = 0.0;
                for (const auto& s : pan.geom.subpanels)
                    subArea += s.area;
                CHECK(subArea ==
                      Catch::Approx(pan.geom.area).epsilon(1e-10));
            }
        }
    }
}
