// Tests for the C++ A502 port.
//
// The first block covers the $points network reader: it must follow the same
// card sequence as meshp/inbc in panair.f90 (kn card, kt card, "nm nn
// netname" cards, then exactly nm*nn coordinate points).

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "AeroModel.h"
#include "Compressibility.h"
#include "BoundaryConditions.h"
#include "FlowSolution.h"
#include "GlobalBoundarySystem.h"
#include "GlobalInfluence.h"
#include "Influence.h"
#include "Matching.h"
#include "NearField.h"
#include "NetworkInfluence.h"
#include "PanelGeometry.h"
#include "PanelParameters.h"
#include "OuterSplines.h"
#include "EdgeSplines.h"
#include "SubpanelSplines.h"
#include "TriangleInterpolation.h"

#include <algorithm>
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


TEST_CASE("compressibility: freestream metric", "[compressibility]") {
    SECTION("M = 0") {
        const auto c = a502::Compressibility::fromAngles(0.0, 0.0, 0.0);
        CHECK(c.betams == Catch::Approx(1.0));
        CHECK(c.abetms == Catch::Approx(1.0));
        CHECK(c.betam == Catch::Approx(1.0));
        CHECK(c.sbetam == Catch::Approx(1.0));
        CHECK(c.btsqi == Catch::Approx(1.0));
        CHECK(c.compd.isApprox(Eigen::Vector3d(1, 0, 0)));
        CHECK(c.arotc.isApprox(Eigen::Matrix3d::Identity()));
        CHECK(c.ggcp.isApprox(Eigen::Matrix3d::Identity()));
        CHECK(c.ggcpit.isApprox(Eigen::Matrix3d::Identity()));
    }
    SECTION("M = 0.6, alpha_c = 4 deg") {
        const auto c = a502::Compressibility::fromAngles(0.6, 4.0, 0.0);
        CHECK(c.betams == Catch::Approx(0.64));
        CHECK(c.betam == Catch::Approx(0.8));
        CHECK(c.btsqi == Catch::Approx(1.0 / 0.64));
        const double a = 4.0 * 3.14159265358979323846 / 180.0;
        CHECK(c.compd.isApprox(Eigen::Vector3d(std::cos(a), 0.0,
                                                std::sin(a))));
        CHECK((c.arotc * c.arotc.transpose() -
               Eigen::Matrix3d::Identity()).norm() < 1e-14);
        CHECK((c.arotci - c.arotc.transpose()).norm() < 1e-14);
    }
    SECTION("M = 1.2") {
        const auto c = a502::Compressibility::fromAngles(1.2, 0.0, 0.0);
        CHECK(c.betams == Catch::Approx(-0.44));
        CHECK(c.abetms == Catch::Approx(0.44));
        CHECK(c.betam == Catch::Approx(std::sqrt(0.44)));
        CHECK(c.sbetam == Catch::Approx(-1.0));
        CHECK(c.btsqi == Catch::Approx(1.0 / 0.44));
    }
    SECTION("M = 1 is rejected") {
        REQUIRE_THROWS_AS(a502::Compressibility::fromAngles(1.0, 0.0, 0.0),
                          std::invalid_argument);
    }
}

TEST_CASE("refloc: scaled local frame", "[compressibility][refloc]") {
    const Eigen::Vector3d n(0, 0, 1);
    const Eigen::Vector3d gen(1, 0, 0);

    SECTION("subsonic") {
        const auto c = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
        const auto f = a502::buildReflocFrame(c, n, gen);
        REQUIRE(f.valid);
        CHECK(f.inclination == 1);
        CHECK((f.toLocal * f.toGlobal -
               Eigen::Matrix3d::Identity()).norm() < 1e-12);
        CHECK((f.toGlobal * f.toLocal -
               Eigen::Matrix3d::Identity()).norm() < 1e-12);
        // The subpanel plane z = 0 stays at local z = 0.
        for (const auto& p : {Eigen::Vector3d(0, 0, 0),
                              Eigen::Vector3d(1, 0, 0),
                              Eigen::Vector3d(0, 1, 0),
                              Eigen::Vector3d(2, 3, 0)}) {
            CHECK(std::abs((f.toLocal * p).z()) < 1e-12);
        }
        CHECK(f.toLocal(1, 1) == Catch::Approx(c.betam));
        CHECK(f.jacobian == Catch::Approx(c.btsqi * c.betam));
    }

    SECTION("supersonic superinclined") {
        const auto c = a502::Compressibility::fromAngles(1.2, 0.0, 0.0);
        const auto f = a502::buildReflocFrame(c, Eigen::Vector3d(1, 0, 0),
                                              Eigen::Vector3d(0, 1, 0));
        REQUIRE(f.valid);
        CHECK(f.inclination == -1);
        CHECK((f.toLocal * f.toGlobal -
               Eigen::Matrix3d::Identity()).norm() < 1e-12);
        for (const auto& p : {Eigen::Vector3d(0, 0, 0),
                              Eigen::Vector3d(0, 1, 0),
                              Eigen::Vector3d(0, 0, 1)}) {
            CHECK(std::abs((f.toLocal * p).z()) < 1e-12);
        }
    }
}

TEST_CASE("panel scaled frames", "[compressibility][refloc]") {
    for (const char* name : {"swb.inp", "vss.inp"}) {
        a502::AeroModel model;
        model.parse(samplePath(name).string());
        REQUIRE(model.solve());
        INFO("deck " << name << ", M=" << model.flow.mach);
        for (const auto& pan : model.panels) {
            for (const auto& s : pan.geom.subpanels) {
                if (s.degenerate)
                    continue;
                CHECK(s.scaledValid);
                if (!s.scaledValid)
                    continue;
                CHECK((s.scaledToLocal * s.scaledToGlobal -
                       Eigen::Matrix3d::Identity()).norm() < 1e-10);
                for (const auto& v : s.vertex) {
                    const Eigen::Vector3d local =
                        s.scaledToLocal * (v - s.origin);
                    CHECK(std::abs(local.z()) < 1e-9);
                }
            }
        }
    }
}

TEST_CASE("influence: compressible panel diameter and PIC level",
          "[influence]") {
    a502::Network net;
    net.netname = "square";
    net.nm = 2;
    net.nn = 2;
    net.coordinates = {0.0, 0.0, 0.0,
                       1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0,
                       1.0, 1.0, 0.0};

    const auto panel =
        a502::buildPanelGeometry(net.coordinates, 2, 2, 0, 0);
    REQUIRE_FALSE(panel.degenerate);

    const auto m0 = a502::Compressibility::fromAngles(0.0, 0.0, 0.0);
    CHECK(a502::compressedPanelDiameter(panel, m0) ==
          Catch::Approx(std::sqrt(2.0)).epsilon(1e-12));

    const Eigen::Vector3d center = panel.cp[8];
    auto atHeight = [&](double h) -> Eigen::Vector3d {
        return center + Eigen::Vector3d(0.0, 0.0, h);
    };

    // eps = sqrt(2)/(2 h) for this panel at M = 0.
    CHECK(a502::influenceLevel(panel, m0, center) == 6);
    CHECK(a502::influenceLevel(panel, m0, atHeight(0.5)) == 6);
    CHECK(a502::influenceLevel(panel, m0, atHeight(1.0)) == 5);
    CHECK(a502::influenceLevel(panel, m0, atHeight(2.0)) == 4);
    CHECK(a502::influenceLevel(panel, m0, atHeight(5.0)) == 3);
    CHECK(a502::influenceLevel(panel, m0, atHeight(10.0)) == 2);
    CHECK(a502::influenceLevel(panel, m0, atHeight(100.0)) == 1);

    // Compressibility changes the diameter but not the near/far ordering.
    const auto m06 = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    CHECK(a502::compressedPanelDiameter(panel, m06) ==
          Catch::Approx(2.0 * std::sqrt(0.41)).epsilon(1e-12));
    CHECK(a502::influenceLevel(panel, m06, atHeight(0.1)) == 6);
    CHECK(a502::influenceLevel(panel, m06, atHeight(100.0)) == 1);
}

TEST_CASE("subsonic near-field influence: constant source and doublet",
          "[nearfield]") {
    const std::array<Eigen::Vector3d, 4> p = {
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0)};
    const Eigen::Vector3d x(0.5, 0.5, 1.0);
    const auto infl =
        a502::subsonicSubpanelInfluence(p, 4, 0, 3, x, 1.0, 6);
    REQUIRE(infl.valid);

    // 8-point Gauss-Legendre quadrature on the unit square.
    constexpr double nodes[8] = {
        -0.9602898564975363, -0.7966664774136267,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363};
    constexpr double weights[8] = {
        0.1012285362903763, 0.2223810344533745,
        0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763};
    constexpr double pi = 3.14159265358979323846;
    const double fourPi = 4.0 * pi;

    double sourcePotential = 0.0;
    double sourceVz = 0.0;
    double doubletPotential = 0.0;
    double doubletX3 = 0.0;

    for (int j = 0; j < 8; ++j) {
        for (int i = 0; i < 8; ++i) {
            const double xp = 0.5 * (1.0 + nodes[i]);
            const double yp = 0.5 * (1.0 + nodes[j]);
            const double w = 0.25 * weights[i] * weights[j];
            const double rx = x.x() - xp;
            const double ry = x.y() - yp;
            const double rz = x.z();
            const double r2 = rx * rx + ry * ry + rz * rz;
            const double r = std::sqrt(r2);
            sourcePotential += w / (fourPi * r);
            sourceVz += w * (-rz) / (fourPi * r2 * r);
            doubletPotential += w * rz / (fourPi * r2 * r);
            doubletX3 += w * (xp * xp * xp / 6.0) * rz /
                         (fourPi * r2 * r);
        }
    }

    // PAN AIR's source potential is minus the 1/(4*pi*r) sheet integral.
    CHECK(infl.source(0, 0) ==
          Catch::Approx(-sourcePotential).epsilon(1e-9));
    CHECK(infl.source(1, 0) == 0.0);
    CHECK(infl.source(2, 0) == 0.0);
    CHECK(infl.source(3, 0) == Catch::Approx(-sourceVz).epsilon(1e-9));
    CHECK(infl.doublet(0, 0) ==
          Catch::Approx(doubletPotential).epsilon(1e-9));
    // subsbi returns no velocity influence for the constant doublet term;
    // its derivative is represented through the source-like terms.
    CHECK(infl.doublet(3, 0) == 0.0);

    // Cubic terms: the lower-order coefficients must agree with the
    // quadratic fit, and the x^3/6 coefficients match direct quadrature.
    const auto infl10 =
        a502::subsonicSubpanelInfluence(p, 4, 0, 3, x, 1.0, 10);
    REQUIRE(infl10.valid);
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 3; ++j) {
            CHECK(infl10.source(i, j) ==
                  Catch::Approx(infl.source(i, j)).epsilon(1e-9));
        }
        for (int j = 0; j < 6; ++j) {
            CHECK(infl10.doublet(i, j) ==
                  Catch::Approx(infl.doublet(i, j)).epsilon(1e-9));
        }
    }
    CHECK(infl10.doublet(0, 6) ==
          Catch::Approx(doubletX3).epsilon(1e-8));

    // The source-only and doublet-only modes leave the other block zero.
    const auto sourceOnly =
        a502::subsonicSubpanelInfluence(p, 4, 0, 1, x, 1.0, 6);
    REQUIRE(sourceOnly.valid);
    CHECK(sourceOnly.doublet.isZero(0.0));
    const auto doubletOnly =
        a502::subsonicSubpanelInfluence(p, 4, 0, 2, x, 1.0, 6);
    REQUIRE(doubletOnly.valid);
    CHECK(doubletOnly.doublet(0, 0) ==
          Catch::Approx(doubletPotential).epsilon(1e-9));
}

TEST_CASE("subsonic near-field influence: triangular subpanel",
          "[nearfield]") {
    const std::array<Eigen::Vector3d, 4> p = {
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 0.0)};
    const Eigen::Vector3d x(1.0 / 3.0, 1.0 / 3.0, 1.0);
    const auto infl =
        a502::subsonicSubpanelInfluence(p, 3, 0, 3, x, 1.0, 6);
    REQUIRE(infl.valid);

    constexpr double nodes[8] = {
        -0.9602898564975363, -0.7966664774136267,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363};
    constexpr double weights[8] = {
        0.1012285362903763, 0.2223810344533745,
        0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763};
    constexpr double pi = 3.14159265358979323846;
    const double fourPi = 4.0 * pi;

    double sourcePotential = 0.0;
    double sourceVz = 0.0;
    double doubletPotential = 0.0;
    for (int j = 0; j < 8; ++j) {
        for (int i = 0; i < 8; ++i) {
            const double u = 0.5 * (1.0 + nodes[i]);
            const double v = 0.5 * (1.0 + nodes[j]);
            const double xp = u;
            const double yp = v * (1.0 - u);
            const double jac = 1.0 - u;
            const double w = 0.25 * weights[i] * weights[j] * jac;
            const double rx = x.x() - xp;
            const double ry = x.y() - yp;
            const double rz = x.z();
            const double r2 = rx * rx + ry * ry + rz * rz;
            const double r = std::sqrt(r2);
            sourcePotential += w / (fourPi * r);
            sourceVz += w * (-rz) / (fourPi * r2 * r);
            doubletPotential += w * rz / (fourPi * r2 * r);
        }
    }

    CHECK(infl.source(0, 0) ==
          Catch::Approx(-sourcePotential).epsilon(1e-9));
    CHECK(infl.source(3, 0) ==
          Catch::Approx(-sourceVz).epsilon(1e-8));
    CHECK(infl.doublet(0, 0) ==
          Catch::Approx(doubletPotential).epsilon(1e-8));
}

TEST_CASE("subsonic near-field influence: panel-centre self term",
          "[nearfield]") {
    a502::Network net;
    net.netname = "square";
    net.nm = 2;
    net.nn = 2;
    net.coordinates = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0, 1.0, 1.0, 0.0};
    const auto flow = a502::Compressibility::fromAngles(0.0, 0.0, 0.0);
    auto panel = a502::buildPanelGeometry(net.coordinates, 2, 2, 0, 0);
    REQUIRE_FALSE(panel.degenerate);
    a502::applyCompressibleFrames(panel, flow);

    const auto infl = a502::subsonicPanelNearFieldInfluence(
        panel, flow, panel.cp[8], 6);
    REQUIRE(infl.valid);
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 5; ++c)
            CHECK(std::isfinite(infl.source(r, c)));
    for (int r = 0; r < 4; ++r)
        for (int c = 0; c < 9; ++c)
            CHECK(std::isfinite(infl.doublet(r, c)));

    // subsbi uses sign(h) for the solid-angle term.  At h = 0 the raw
    // kernel therefore returns a zero constant-doublet potential; it does
    // not include the surface jump.
    const double jump =
        infl.doublet.row(0).dot(a502::PanelDoubletValues::Ones());
    CHECK(jump == Catch::Approx(0.0).margin(1e-12));
}

TEST_CASE("subsonic near-field influence: whole panel", "[nearfield]") {
    a502::Network net;
    net.netname = "square";
    net.nm = 2;
    net.nn = 2;
    net.coordinates = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0, 1.0, 1.0, 0.0};
    const auto flow = a502::Compressibility::fromAngles(0.0, 0.0, 0.0);
    auto panel = a502::buildPanelGeometry(net.coordinates, 2, 2, 0, 0);
    REQUIRE_FALSE(panel.degenerate);
    a502::applyCompressibleFrames(panel, flow);

    const Eigen::Vector3d x(0.3, 0.7, 1.0);
    const auto infl =
        a502::subsonicPanelNearFieldInfluence(panel, flow, x, 6);
    REQUIRE(infl.valid);

    // Gauss-Legendre quadrature over the unit square.
    constexpr double nodes[8] = {
        -0.9602898564975363, -0.7966664774136267,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363};
    constexpr double weights[8] = {
        0.1012285362903763, 0.2223810344533745,
        0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763};
    constexpr double pi = 3.14159265358979323846;
    const double fourPi = 4.0 * pi;
    const double normalZ = panel.normal.z();

    auto mu = [](const Eigen::Vector3d& p) {
        return 1.0 + 2.0 * p.x() + 3.0 * p.y() +
               4.0 * p.x() * p.x() / 2.0 + 5.0 * p.x() * p.y() +
               6.0 * p.y() * p.y() / 2.0;
    };
    auto sigma = [](const Eigen::Vector3d& p) {
        return 1.0 + 2.0 * p.x() + 3.0 * p.y();
    };

    double sourcePotential = 0.0; // -int sigma/(4 pi r)
    double sourceVx = 0.0;
    double sourceVy = 0.0;
    double sourceVz = 0.0;
    double constantSourcePotential = 0.0;
    double doubletPotential = 0.0;
    double constantDoubletPotential = 0.0;

    for (int j = 0; j < 8; ++j) {
        for (int i = 0; i < 8; ++i) {
            const double xp = 0.5 * (1.0 + nodes[i]);
            const double yp = 0.5 * (1.0 + nodes[j]);
            const double w = 0.25 * weights[i] * weights[j];
            const Eigen::Vector3d p(xp, yp, 0.0);
            const double rx = x.x() - xp;
            const double ry = x.y() - yp;
            const double rz = x.z();
            const double r2 = rx * rx + ry * ry + rz * rz;
            const double r = std::sqrt(r2);
            const double s = sigma(p);
            sourcePotential -= w * s / (fourPi * r);
            sourceVx += w * s * rx / (fourPi * r2 * r);
            sourceVy += w * s * ry / (fourPi * r2 * r);
            sourceVz += w * s * rz / (fourPi * r2 * r);
            constantSourcePotential -= w / (fourPi * r);
            doubletPotential +=
                w * mu(p) * normalZ * rz / (fourPi * r2 * r);
            constantDoubletPotential +=
                w * normalZ * rz / (fourPi * r2 * r);
        }
    }

    a502::PanelSourceValues sigmaValues;
    for (int i = 0; i < 4; ++i)
        sigmaValues(i) = sigma(panel.cp[static_cast<std::size_t>(i)]);
    sigmaValues(4) = sigma(panel.cp[8]);

    a502::PanelDoubletValues muValues;
    for (int i = 0; i < 9; ++i)
        muValues(i) = mu(panel.cp[static_cast<std::size_t>(i)]);

    const Eigen::VectorXd sourceHit = infl.source * sigmaValues;
    CHECK(sourceHit(0) ==
          Catch::Approx(sourcePotential).epsilon(1e-7));
    CHECK(sourceHit(1) == Catch::Approx(sourceVx).epsilon(1e-7));
    CHECK(sourceHit(2) == Catch::Approx(sourceVy).epsilon(1e-7));
    CHECK(sourceHit(3) == Catch::Approx(sourceVz).epsilon(1e-7));

    const double constantSource =
        infl.source.row(0).dot(a502::PanelSourceValues::Ones());
    CHECK(constantSource ==
          Catch::Approx(constantSourcePotential).epsilon(1e-7));

    const double doubletHit = (infl.doublet * muValues)(0);
    CHECK(doubletHit ==
          Catch::Approx(doubletPotential).epsilon(1e-7));

    // The intermediate (quasi-far) single-panel evaluation must also be
    // exact for constant source and doublet distributions.
    const auto intermediate =
        a502::subsonicPanelIntermediateInfluence(panel, flow, x, 6);
    REQUIRE(intermediate.valid);
    CHECK(intermediate.source(0, 0) ==
          Catch::Approx(constantSourcePotential).epsilon(1e-7));
    CHECK(intermediate.doublet(0, 0) ==
          Catch::Approx(constantDoubletPotential).epsilon(1e-7));
}

TEST_CASE("subsonic near-field influence: collapsed-edge panel",
          "[nearfield]") {
    const std::vector<double> coords = {
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 1.0, 0.0};
    const auto flow = a502::Compressibility::fromAngles(0.0, 0.0, 0.0);
    auto panel = a502::buildPanelGeometry(coords, 2, 2, 0, 0);
    REQUIRE_FALSE(panel.degenerate);
    REQUIRE(panel.collapsedEdge == 1);
    a502::applyCompressibleFrames(panel, flow);

    const Eigen::Vector3d x(0.5, 0.25, 1.0);
    const auto infl =
        a502::subsonicPanelNearFieldInfluence(panel, flow, x, 6);
    REQUIRE(infl.valid);

    // Triangle (0,0), (1,0), (1,1): x=u, y=u*v with Jacobian u.
    constexpr double nodes[8] = {
        -0.9602898564975363, -0.7966664774136267,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363};
    constexpr double weights[8] = {
        0.1012285362903763, 0.2223810344533745,
        0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763};
    constexpr double pi = 3.14159265358979323846;
    const double fourPi = 4.0 * pi;

    double sourcePotential = 0.0;
    double doubletPotential = 0.0;
    for (int j = 0; j < 8; ++j) {
        for (int i = 0; i < 8; ++i) {
            const double u = 0.5 * (1.0 + nodes[i]);
            const double v = 0.5 * (1.0 + nodes[j]);
            const double xp = u;
            const double yp = u * v;
            const double w = 0.25 * weights[i] * weights[j] * u;
            const double rx = x.x() - xp;
            const double ry = x.y() - yp;
            const double rz = x.z();
            const double r2 = rx * rx + ry * ry + rz * rz;
            const double r = std::sqrt(r2);
            sourcePotential -= w / (fourPi * r);
            doubletPotential +=
                w * panel.normal.z() * rz / (fourPi * r2 * r);
        }
    }

    CHECK(infl.source.row(0).dot(a502::PanelSourceValues::Ones()) ==
          Catch::Approx(sourcePotential).epsilon(1e-7));
    CHECK(infl.doublet.row(0).dot(a502::PanelDoubletValues::Ones()) ==
          Catch::Approx(doubletPotential).epsilon(1e-7));
}

TEST_CASE("network influence: outer splines and near-field PIC",
          "[networkinfluence]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto flow = a502::Compressibility::fromAngles(0.0, 0.0, 0.0);
    const auto layout = a502::buildNetworkSingularityLayout(net);
    REQUIRE(layout.source.count == 49);
    REQUIRE(layout.doublet.count == 81);

    const Eigen::Vector3d x(3.5, 3.5, 1.0);
    const auto infl =
        a502::evaluateNetworkNearFieldInfluence(net, flow, x, layout);
    REQUIRE(infl.valid);
    REQUIRE(infl.source.cols() == 49);
    REQUIRE(infl.doublet.cols() == 81);

    auto sigma = [](const Eigen::Vector3d& p) {
        return 1.0 + 0.2 * p.x() - 0.3 * p.y();
    };
    auto mu = [](const Eigen::Vector3d& p) {
        return 1.0 + 2.0 * p.x() + 3.0 * p.y() +
               4.0 * p.x() * p.x() / 2.0 + 5.0 * p.x() * p.y() +
               6.0 * p.y() * p.y() / 2.0;
    };

    Eigen::VectorXd sourceValues(layout.source.count);
    for (int k = 0; k < layout.source.count; ++k)
        sourceValues(k) = sigma(
            layout.source.location[static_cast<std::size_t>(k)]);

    Eigen::VectorXd doubletValues(layout.doublet.count);
    for (int k = 0; k < layout.doublet.count; ++k)
        doubletValues(k) = mu(
            layout.doublet.location[static_cast<std::size_t>(k)]);

    // Direct panel-by-panel Gauss-Legendre quadrature.
    constexpr double nodes[8] = {
        -0.9602898564975363, -0.7966664774136267,
        -0.5255324099163290, -0.1834346424956498,
         0.1834346424956498,  0.5255324099163290,
         0.7966664774136267,  0.9602898564975363};
    constexpr double weights[8] = {
        0.1012285362903763, 0.2223810344533745,
        0.3137066458778873, 0.3626837833783620,
        0.3626837833783620, 0.3137066458778873,
        0.2223810344533745, 0.1012285362903763};
    constexpr double pi = 3.14159265358979323846;
    const double fourPi = 4.0 * pi;

    double sourcePotential = 0.0;
    double sourceVx = 0.0;
    double sourceVy = 0.0;
    double sourceVz = 0.0;
    double doubletPotential = 0.0;

    const auto basePanel =
        a502::buildPanelGeometry(net.coordinates, 8, 8, 0, 0);
    const double normalZ = basePanel.normal.z();

    for (int pj = 0; pj < 7; ++pj) {
        for (int pi = 0; pi < 7; ++pi) {
            for (int j = 0; j < 8; ++j) {
                for (int i = 0; i < 8; ++i) {
                    const double xp =
                        pi + 0.5 * (1.0 + nodes[i]);
                    const double yp =
                        pj + 0.5 * (1.0 + nodes[j]);
                    const double w = 0.25 * weights[i] * weights[j];
                    const Eigen::Vector3d p(xp, yp, 0.0);
                    const double rx = x.x() - xp;
                    const double ry = x.y() - yp;
                    const double rz = x.z();
                    const double r2 = rx * rx + ry * ry + rz * rz;
                    const double r = std::sqrt(r2);
                    const double s = sigma(p);
                    const double m = mu(p);
                    sourcePotential -= w * s / (fourPi * r);
                    sourceVx += w * s * rx / (fourPi * r2 * r);
                    sourceVy += w * s * ry / (fourPi * r2 * r);
                    sourceVz += w * s * rz / (fourPi * r2 * r);
                    doubletPotential +=
                        w * m * normalZ * rz / (fourPi * r2 * r);
                }
            }
        }
    }

    const Eigen::Vector4d sourceHit = infl.source * sourceValues;
    CHECK(sourceHit(0) ==
          Catch::Approx(sourcePotential).epsilon(1e-6));
    CHECK(sourceHit(1) == Catch::Approx(sourceVx).epsilon(1e-6));
    CHECK(sourceHit(2) == Catch::Approx(sourceVy).epsilon(1e-6));
    CHECK(sourceHit(3) == Catch::Approx(sourceVz).epsilon(1e-6));

    const double doubletHit = (infl.doublet * doubletValues)(0);
    CHECK(doubletHit ==
          Catch::Approx(doubletPotential).epsilon(1e-6));
}

TEST_CASE("network AIC: rows at all control points", "[networkinfluence]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto flow = a502::Compressibility::fromAngles(0.0, 0.0, 0.0);
    const auto layout = a502::buildNetworkSingularityLayout(net);
    const auto controlPoints = a502::buildNetworkControlPoints(net);
    REQUIRE(controlPoints.size() == 81);

    std::vector<Eigen::Vector3d> positions;
    positions.reserve(controlPoints.size());
    for (const auto& cp : controlPoints)
        positions.push_back(cp.position);

    const auto aic =
        a502::buildNetworkNearFieldAic(net, flow, layout, positions);
    REQUIRE(aic.valid);
    CHECK(aic.nControlPoints == 81);
    CHECK(aic.nSource == 49);
    CHECK(aic.nDoublet == 81);
    REQUIRE(aic.source.rows() == 4 * 81);
    REQUIRE(aic.source.cols() == 49);
    REQUIRE(aic.doublet.rows() == 4 * 81);
    REQUIRE(aic.doublet.cols() == 81);

    CHECK(aic.source.allFinite());
    CHECK(aic.doublet.allFinite());
}

TEST_CASE("network singularity types follow kt", "[networkinfluence]") {
    CHECK(a502::networkSingularityTypes(1).sourceType == 1);
    CHECK(a502::networkSingularityTypes(1).doubletType == 12);
    CHECK(a502::networkSingularityTypes(7).sourceType == 1);
    CHECK(a502::networkSingularityTypes(7).doubletType == 0);
    CHECK(a502::networkSingularityTypes(12).sourceType == 0);
    CHECK(a502::networkSingularityTypes(12).doubletType == 12);
    CHECK(a502::networkSingularityTypes(18).doubletType == 18);
    CHECK(a502::networkSingularityTypes(19).doubletType == 6);
    CHECK(a502::networkSingularityTypes(20).doubletType == 20);

    a502::Network net;
    net.netname = "grid";
    net.nm = 3;
    net.nn = 3;
    for (int n = 0; n < 3; ++n) {
        for (int m = 0; m < 3; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    net.kt = 12;
    const auto doubletsOnly =
        a502::buildNetworkSingularityLayout(net);
    CHECK_FALSE(doubletsOnly.hasSource);
    CHECK(doubletsOnly.hasDoublet);

    net.kt = 7;
    const auto sourceOnly =
        a502::buildNetworkSingularityLayout(net);
    CHECK(sourceOnly.hasSource);
    CHECK_FALSE(sourceOnly.hasDoublet);

    net.kt = 1;
    const auto composite =
        a502::buildNetworkSingularityLayout(net);
    CHECK(composite.hasSource);
    CHECK(composite.hasDoublet);

    a502::Network doubletNet = net;
    doubletNet.kt = 12;
    a502::Network compositeNet = net;
    compositeNet.kt = 1;
    const std::vector<a502::Network> mixed = {doubletNet, compositeNet};
    const auto global = a502::buildGlobalSingularityLayout(mixed);
    REQUIRE(global.network.size() == 2);
    CHECK(global.activeSourceOffset[0] == -1);
    CHECK(global.activeSourceOffset[1] == 0);
    CHECK(global.activeDoubletOffset[0] == 0);
    CHECK(global.activeDoubletOffset[1] == 16);
    CHECK(global.activeSource == 4);
    CHECK(global.activeDoublet == 32);
}

TEST_CASE("network type defaults follow inbc", "[networkinfluence]") {
    const auto kt1 = a502::networkTypeDefaults(1);
    CHECK(kt1.sourceType == 1);
    CHECK(kt1.doubletType == 12);
    CHECK(kt1.ipot == 2);
    CHECK(kt1.nlopt1 == 5);
    CHECK(kt1.nropt1 == 3);
    CHECK(kt1.nlopt2 == 7);
    CHECK(kt1.nropt2 == 2);

    const auto kt2 = a502::networkTypeDefaults(2);
    CHECK(kt2.sourceType == 0);
    CHECK(kt2.doubletType == 12);
    CHECK(kt2.ipot == 0);
    CHECK(kt2.nlopt1 == 0);
    CHECK(kt2.nropt1 == 0);
    CHECK(kt2.nlopt2 == 4);
    CHECK(kt2.nropt2 == 3);

    const auto kt3 = a502::networkTypeDefaults(3);
    CHECK(kt3.nlopt1 == 3);
    CHECK(kt3.nropt1 == 2);

    const auto kt4 = a502::networkTypeDefaults(4);
    CHECK(kt4.nlopt1 == 2);
    CHECK(kt4.nropt1 == 1);

    const auto baseSubsonic = a502::networkTypeDefaults(5, 0.6);
    CHECK(baseSubsonic.nlopt1 == 6);
    CHECK(baseSubsonic.nropt1 == 9);
    CHECK(baseSubsonic.nlopt2 == 7);
    const auto baseSupersonic = a502::networkTypeDefaults(5, 1.2);
    CHECK(baseSupersonic.nlopt1 == 2);
    CHECK(baseSupersonic.nropt1 == 2);
    CHECK(baseSupersonic.nlopt2 == 6);
    CHECK(baseSupersonic.nropt2 == 2);

    const auto kt6 = a502::networkTypeDefaults(6);
    CHECK(kt6.sourceType == 1);
    CHECK(kt6.doubletType == 0);
    CHECK(kt6.nlopt1 == 5);
    CHECK(kt6.nropt1 == 2);

    const auto kt7 = a502::networkTypeDefaults(7);
    CHECK(kt7.sourceType == 1);
    CHECK(kt7.doubletType == 0);
    CHECK(kt7.nlopt1 == 2);
    CHECK(kt7.nropt1 == 3);

    const auto kt8 = a502::networkTypeDefaults(8);
    CHECK(kt8.sourceType == 1);
    CHECK(kt8.doubletType == 12);
    CHECK(kt8.nlopt1 == 2);

    const auto kt9 = a502::networkTypeDefaults(9);
    CHECK(kt9.nropt1 == 7);

    const auto kt11 = a502::networkTypeDefaults(11);
    CHECK(kt11.nlopt1 == 14);

    const auto kt12 = a502::networkTypeDefaults(12);
    CHECK(kt12.sourceType == 0);
    CHECK(kt12.doubletType == 12);
    CHECK(kt12.ipot == 0);
    CHECK(kt12.nlopt1 == 0);
    CHECK(kt12.nropt1 == 3);
    CHECK(kt12.nlopt2 == 13);
    CHECK(kt12.nropt2 == 3);

    const auto kt13 = a502::networkTypeDefaults(13, 1.5);
    CHECK(kt13.nlopt1 == 12);
    CHECK(kt13.nropt1 == 2);

    const auto kt14 = a502::networkTypeDefaults(14);
    CHECK(kt14.nlopt1 == 11);
    CHECK(kt14.nropt1 == 1);

    const auto kt15 = a502::networkTypeDefaults(15);
    CHECK(kt15.nlopt1 == 11);
    CHECK(kt15.nropt1 == 3);

    const auto kt17 = a502::networkTypeDefaults(17);
    CHECK(kt17.sourceType == 1);
    CHECK(kt17.doubletType == 0);
    CHECK(kt17.nlopt1 == 11);

    const auto kt18 = a502::networkTypeDefaults(18);
    CHECK(kt18.sourceType == 0);
    CHECK(kt18.doubletType == 18);
    CHECK(kt18.ipot == 2);
    CHECK(kt18.nlopt1 == 0);
    CHECK(kt18.nlopt2 == 15);
    CHECK(kt18.nropt2 == 2);

    const auto kt19 = a502::networkTypeDefaults(19);
    CHECK(kt19.sourceType == 0);
    CHECK(kt19.doubletType == 6);
    CHECK(kt19.nlopt2 == 18);

    const auto kt20 = a502::networkTypeDefaults(20);
    CHECK(kt20.sourceType == 0);
    CHECK(kt20.doubletType == 20);
    CHECK(kt20.nlopt2 == 6);

    const auto kt21 = a502::networkTypeDefaults(21);
    CHECK(kt21.sourceType == 0);
    CHECK(kt21.doubletType == 12);
    CHECK(kt21.ipot == 4);
    CHECK(kt21.nropt2 == 4);

    // Thin wakes are doublet-only by default.
    CHECK(a502::networkSingularityTypes(2).sourceType == 0);
    CHECK(a502::networkSingularityTypes(2).doubletType == 12);

    // Unknown kt codes keep the permissive composite defaults.
    const auto unknown = a502::networkTypeDefaults(99);
    CHECK(unknown.sourceType == 1);
    CHECK(unknown.doubletType == 12);
    CHECK(unknown.nlopt1 == 5);
}

TEST_CASE("network influence: inactive singularity types are skipped",
          "[networkinfluence]") {
    auto makeGrid = [](int kt) {
        a502::Network net;
        net.netname = "grid";
        net.nm = 4;
        net.nn = 4;
        net.kt = kt;
        for (int n = 0; n < 4; ++n) {
            for (int m = 0; m < 4; ++m) {
                net.coordinates.push_back(static_cast<double>(m));
                net.coordinates.push_back(static_cast<double>(n));
                net.coordinates.push_back(0.0);
            }
        }
        return net;
    };

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const Eigen::Vector3d x(1.5, 1.5, 1.0);

    const a502::Network sourceOnly = makeGrid(7);
    const auto sourceLayout =
        a502::buildNetworkSingularityLayout(sourceOnly);
    REQUIRE(sourceLayout.hasSource);
    REQUIRE_FALSE(sourceLayout.hasDoublet);
    const auto sourceInfl = a502::evaluateNetworkNearFieldInfluence(
        sourceOnly, flow, x, sourceLayout);
    REQUIRE(sourceInfl.valid);
    REQUIRE(sourceInfl.source.cols() == sourceLayout.source.count);
    REQUIRE(sourceInfl.doublet.cols() == sourceLayout.doublet.count);
    CHECK_FALSE(sourceInfl.source.isZero(1e-12));
    CHECK(sourceInfl.doublet.isZero(0.0));

    const a502::Network doubletOnly = makeGrid(12);
    const auto doubletLayout =
        a502::buildNetworkSingularityLayout(doubletOnly);
    REQUIRE_FALSE(doubletLayout.hasSource);
    REQUIRE(doubletLayout.hasDoublet);
    const auto doubletInfl = a502::evaluateNetworkNearFieldInfluence(
        doubletOnly, flow, x, doubletLayout);
    REQUIRE(doubletInfl.valid);
    CHECK(doubletInfl.source.isZero(0.0));
    CHECK_FALSE(doubletInfl.doublet.isZero(1e-12));

    // Skipping the inactive type must not change the active columns.
    const a502::Network composite = makeGrid(1);
    const auto compositeLayout =
        a502::buildNetworkSingularityLayout(composite);
    const auto compositeInfl = a502::evaluateNetworkNearFieldInfluence(
        composite, flow, x, compositeLayout);
    REQUIRE(compositeInfl.valid);
    CHECK((sourceInfl.source - compositeInfl.source).norm() < 1e-12);
    CHECK((doubletInfl.doublet - compositeInfl.doublet).norm() < 1e-12);

    // Local terms follow the same type flags.
    auto panelCenter = [](const a502::Network& net) {
        const auto panel = a502::buildPanelGeometry(
            net.coordinates, net.nm, net.nn, 1, 1);
        return panel.cp[8];
    };

    const auto sourceLocal = a502::buildPanelPointLocalTerms(
        sourceOnly, flow, 1, 1, panelCenter(sourceOnly),
        sourceLayout.doublet, sourceLayout.hasSource,
        sourceLayout.hasDoublet);
    REQUIRE(sourceLocal.valid);
    CHECK_FALSE(sourceLocal.sourceId.empty());
    CHECK(sourceLocal.sourceValue.size() ==
          static_cast<Eigen::Index>(sourceLocal.sourceId.size()));
    CHECK(sourceLocal.doubletId.empty());
    CHECK(sourceLocal.doubletValue.size() == 0);

    const auto doubletLocal = a502::buildPanelPointLocalTerms(
        doubletOnly, flow, 1, 1, panelCenter(doubletOnly),
        doubletLayout.doublet, doubletLayout.hasSource,
        doubletLayout.hasDoublet);
    REQUIRE(doubletLocal.valid);
    CHECK(doubletLocal.sourceId.empty());
    CHECK(doubletLocal.sourceValue.size() == 0);
    CHECK_FALSE(doubletLocal.doubletId.empty());
    CHECK(doubletLocal.doubletValue.cols() ==
          static_cast<Eigen::Index>(doubletLocal.doubletId.size()));

    // Explicit flags on a composite network select the requested terms.
    const auto onlyDoublet = a502::buildPanelPointLocalTerms(
        composite, flow, 1, 1, panelCenter(composite),
        compositeLayout.doublet, false, true);
    REQUIRE(onlyDoublet.valid);
    CHECK(onlyDoublet.sourceId.empty());
    CHECK_FALSE(onlyDoublet.doubletId.empty());
}

TEST_CASE("global influence: two networks and column offsets",
          "[networkinfluence]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    a502::Network shifted = net;
    shifted.netname = "grid2";
    for (std::size_t k = 0; k < shifted.coordinates.size(); k += 3)
        shifted.coordinates[k] += 20.0;

    const std::vector<a502::Network> networks = {net, shifted};
    const auto layout = a502::buildGlobalSingularityLayout(networks);
    REQUIRE(layout.network.size() == 2);
    REQUIRE(layout.totalSource == 98);
    REQUIRE(layout.totalDoublet == 162);
    CHECK(layout.sourceOffset[0] == 0);
    CHECK(layout.sourceOffset[1] == 49);
    CHECK(layout.doubletOffset[0] == 0);
    CHECK(layout.doubletOffset[1] == 81);

    const auto flow = a502::Compressibility::fromAngles(0.0, 0.0, 0.0);
    const Eigen::Vector3d x(3.5, 3.5, 1.0);
    const auto global =
        a502::evaluateGlobalNearFieldInfluence(networks, flow, x, layout);
    REQUIRE(global.valid);
    REQUIRE(global.source.cols() == 98);
    REQUIRE(global.doublet.cols() == 162);

    const auto local0 = a502::evaluateNetworkNearFieldInfluence(
        networks[0], flow, x, layout.network[0]);
    const auto local1 = a502::evaluateNetworkNearFieldInfluence(
        networks[1], flow, x, layout.network[1]);
    REQUIRE(local0.valid);
    REQUIRE(local1.valid);

    CHECK((global.source.block(0, 0, 4, 49) - local0.source).norm() <
          1e-12);
    CHECK((global.source.block(0, 49, 4, 49) - local1.source).norm() <
          1e-12);
    CHECK((global.doublet.block(0, 0, 4, 81) - local0.doublet).norm() <
          1e-12);
    CHECK((global.doublet.block(0, 81, 4, 81) - local1.doublet).norm() <
          1e-12);

    const std::vector<Eigen::Vector3d> fieldPoints = {
        x, Eigen::Vector3d(5.0, 5.0, 2.0)};
    const auto aic =
        a502::buildGlobalNearFieldAic(networks, flow, layout, fieldPoints);
    REQUIRE(aic.valid);
    CHECK(aic.nControlPoints == 2);
    CHECK(aic.nSource == 98);
    CHECK(aic.nDoublet == 162);
    REQUIRE(aic.source.rows() == 8);
    REQUIRE(aic.doublet.rows() == 8);
    CHECK((aic.source.block(0, 0, 4, 98) - global.source).norm() < 1e-12);
    CHECK((aic.doublet.block(0, 0, 4, 162) - global.doublet).norm() <
          1e-12);
    CHECK(aic.source.allFinite());
    CHECK(aic.doublet.allFinite());
}

TEST_CASE("boundary condition rows: nct=1 and nct=4",
          "[boundaryconditions]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const std::vector<a502::Network> networks = {net};
    const auto globalLayout =
        a502::buildGlobalSingularityLayout(networks);
    REQUIRE(globalLayout.totalSource == 49);
    REQUIRE(globalLayout.totalDoublet == 81);

    auto panel = a502::buildPanelGeometry(net.coordinates, 8, 8, 3, 3);
    REQUIRE_FALSE(panel.degenerate);
    a502::applyCompressibleFrames(panel, flow);

    const Eigen::Vector3d fieldPoint =
        panel.cp[8] + 0.01 * panel.normal;
    const auto influence = a502::evaluateGlobalNearFieldInfluence(
        networks, flow, fieldPoint, globalLayout);
    REQUIRE(influence.valid);

    const auto local = a502::buildPanelCenterLocalTerms(
        net, flow, 3, 3, globalLayout.network[0].doublet);
    REQUIRE(local.valid);

    Eigen::VectorXd sourceValues(49);
    Eigen::VectorXd doubletValues(81);
    for (int k = 0; k < 49; ++k)
        sourceValues(k) = 0.01 * (k + 1);
    for (int k = 0; k < 81; ++k)
        doubletValues(k) = -0.005 * (k + 1);

    const Eigen::Vector3d inducedVelocity =
        influence.source.block(1, 0, 3, 49) * sourceValues +
        influence.doublet.block(1, 0, 3, 81) * doubletValues;

    // nct=1: ca*n.v + cd*local source value.
    a502::BoundaryConditionCoefficients bc1;
    bc1.nct = 1;
    bc1.cu = 0.7;
    bc1.cl = 0.2;
    bc1.bet = 0.3;
    const auto row1 = a502::buildBoundaryRow(
        influence, panel.normal, flow, bc1, local);
    REQUIRE(row1.valid);
    CHECK(row1.rhs == Catch::Approx(0.3).epsilon(1e-14));

    const double ca = bc1.cu + bc1.cl;
    const double cd = 0.5 * (bc1.cu - bc1.cl);
    double localSourceValue = 0.0;
    for (std::size_t k = 0; k < local.sourceId.size(); ++k)
        localSourceValue +=
            local.sourceValue(static_cast<Eigen::Index>(k)) *
            sourceValues(local.sourceId[k]);

    const double got1 = row1.source.dot(sourceValues) +
                        row1.doublet.dot(doubletValues);
    const double want1 =
        ca * panel.normal.dot(inducedVelocity) + cd * localSourceValue;
    CHECK(got1 == Catch::Approx(want1).epsilon(1e-10));

    // nct=4: da*potential + dd*local doublet potential.
    a502::BoundaryConditionCoefficients bc4;
    bc4.nct = 4;
    bc4.du = 0.4;
    bc4.dl = 0.1;
    bc4.bet = -0.2;
    const auto row4 = a502::buildBoundaryRow(
        influence, panel.normal, flow, bc4, local);
    REQUIRE(row4.valid);
    CHECK(row4.rhs == Catch::Approx(-0.2).epsilon(1e-14));

    const double da = bc4.du + bc4.dl;
    const double dd = 0.5 * (bc4.du - bc4.dl);
    const double potential =
        influence.source.row(0).dot(sourceValues) +
        influence.doublet.row(0).dot(doubletValues);
    double localDoubletValue = 0.0;
    for (std::size_t k = 0; k < local.doubletId.size(); ++k)
        localDoubletValue +=
            local.doubletValue(0, static_cast<Eigen::Index>(k)) *
            doubletValues(local.doubletId[k]);

    const double got4 = row4.source.dot(sourceValues) +
                        row4.doublet.dot(doubletValues);
    const double want4 =
        da * potential + dd * localDoubletValue;
    CHECK(got4 == Catch::Approx(want4).epsilon(1e-10));

    // nct=2 reduces to nct=1 when cu == cl and the general/potential
    // coefficients vanish (td = 0).
    a502::BoundaryConditionCoefficients bc1eq;
    bc1eq.nct = 1;
    bc1eq.cu = 1.0;
    bc1eq.cl = 1.0;
    bc1eq.bet = 0.25;
    const auto row1eq = a502::buildBoundaryRow(
        influence, panel.normal, flow, bc1eq, local);
    REQUIRE(row1eq.valid);

    a502::BoundaryConditionCoefficients bc2eq;
    bc2eq.nct = 2;
    bc2eq.cu = 1.0;
    bc2eq.cl = 1.0;
    bc2eq.bet = 0.25;
    const auto row2eq = a502::buildBoundaryRow(
        influence, panel.normal, flow, bc2eq, local);
    REQUIRE(row2eq.valid);
    CHECK((row2eq.source - row1eq.source).norm() < 1e-10);
    CHECK((row2eq.doublet - row1eq.doublet).norm() < 1e-10);
    CHECK(row2eq.rhs == Catch::Approx(0.25).epsilon(1e-14));

    // nct=2 reduces to nct=4 when only the potential part is present.
    a502::BoundaryConditionCoefficients bc4eq;
    bc4eq.nct = 4;
    bc4eq.du = 1.0;
    bc4eq.dl = 1.0;
    bc4eq.bet = -0.1;
    const auto row4eq = a502::buildBoundaryRow(
        influence, panel.normal, flow, bc4eq, local);
    REQUIRE(row4eq.valid);

    a502::BoundaryConditionCoefficients bc2pot;
    bc2pot.nct = 2;
    bc2pot.du = 1.0;
    bc2pot.dl = 1.0;
    bc2pot.bet = -0.1;
    const auto row2pot = a502::buildBoundaryRow(
        influence, panel.normal, flow, bc2pot, local);
    REQUIRE(row2pot.valid);
    CHECK((row2pot.source - row4eq.source).norm() < 1e-10);
    CHECK((row2pot.doublet - row4eq.doublet).norm() < 1e-10);

    a502::BoundaryConditionCoefficients bad;
    bad.nct = 3;
    CHECK_FALSE(a502::buildBoundaryRow(
                    influence, panel.normal, flow, bad, local)
                    .valid);

    // Generic point evaluation must agree with the P9 specialization and
    // with the source outer spline at an edge midpoint and a corner.
    const auto genericCenter = a502::buildPanelPointLocalTerms(
        net, flow, 3, 3, panel.cp[8], globalLayout.network[0].doublet);
    REQUIRE(genericCenter.valid);
    REQUIRE(genericCenter.sourceId == local.sourceId);
    CHECK((genericCenter.sourceValue - local.sourceValue).norm() < 1e-10);
    CHECK((genericCenter.doubletValue - local.doubletValue).norm() < 1e-10);

    const auto edgeLocal = a502::buildPanelPointLocalTerms(
        net, flow, 3, 3, panel.cp[4], globalLayout.network[0].doublet);
    REQUIRE(edgeLocal.valid);
    const auto sourceSpline = a502::buildSourceOuterSpline(net, 3, 3);
    REQUIRE(sourceSpline.valid);
    REQUIRE(edgeLocal.sourceId.size() == 9);
    for (int k = 0; k < 9; ++k) {
        const double expected =
            0.5 * (sourceSpline.coefficients(0, k) +
                   sourceSpline.coefficients(1, k));
        CHECK(edgeLocal.sourceValue(k) ==
              Catch::Approx(expected).epsilon(1e-10));
    }
    CHECK(edgeLocal.doubletValue.allFinite());

    const auto cornerLocal = a502::buildPanelPointLocalTerms(
        net, flow, 3, 3, panel.cp[0], globalLayout.network[0].doublet);
    REQUIRE(cornerLocal.valid);
    REQUIRE(cornerLocal.sourceId.size() == 9);
    for (int k = 0; k < 9; ++k) {
        CHECK(cornerLocal.sourceValue(k) ==
              Catch::Approx(sourceSpline.coefficients(0, k))
                  .epsilon(1e-10));
    }
    CHECK(cornerLocal.doubletValue.allFinite());

    const auto controlPoints = a502::buildNetworkControlPoints(net);
    const auto system = a502::buildNetworkNormalFluxRows(
        net, flow, globalLayout.network[0], controlPoints, 1.0, 0.0);
    REQUIRE(system.valid);
    CHECK(system.nRows == 49);
    CHECK(system.nSource == 49);
    CHECK(system.nDoublet == 81);
    REQUIRE(system.source.rows() == 49);
    REQUIRE(system.doublet.rows() == 49);
    REQUIRE(system.rhs.size() == 49);
    CHECK(system.source.allFinite());
    CHECK(system.doublet.allFinite());
    CHECK(system.rhs.allFinite());

    const a502::ControlPoint& firstCenter =
        controlPoints[static_cast<std::size_t>(
            system.controlPointIndex[0])];
    const double expectedRhs =
        -firstCenter.normal.dot(flow.compd);
    CHECK(system.rhs(0) ==
          Catch::Approx(expectedRhs).epsilon(1e-12));
}

TEST_CASE("matching: shared full edge", "[matching]") {
    auto makeGrid = [](double y0, const char* name) {
        a502::Network net;
        net.netname = name;
        net.nm = 4;
        net.nn = 4;
        for (int n = 0; n < 4; ++n) {
            for (int m = 0; m < 4; ++m) {
                net.coordinates.push_back(static_cast<double>(m));
                net.coordinates.push_back(y0 + static_cast<double>(n));
                net.coordinates.push_back(0.0);
            }
        }
        return net;
    };

    const std::vector<a502::Network> networks = {
        makeGrid(0.0, "a"), makeGrid(3.0, "b")};
    const auto edges = a502::buildNetworkEdges(networks);
    REQUIRE(edges.size() == 8);

    const auto matches = a502::findFullEdgeMatches(edges, 1e-10);
    REQUIRE(matches.size() == 1);
    CHECK(edges[static_cast<std::size_t>(matches[0].edgeA)].network == 0);
    CHECK(edges[static_cast<std::size_t>(matches[0].edgeA)].edge == 2);
    CHECK(edges[static_cast<std::size_t>(matches[0].edgeB)].network == 1);
    CHECK(edges[static_cast<std::size_t>(matches[0].edgeB)].edge == 4);
    CHECK_FALSE(matches[0].reversed);

    const auto layout =
        a502::buildGlobalSingularityLayout(networks);
    REQUIRE(layout.totalDoublet == 50);
    const auto rows =
        a502::buildEdgeDoubletMatchingRows(networks, layout, matches);
    REQUIRE(rows.valid);
    REQUIRE(rows.nRows == 5); // 3 intervals + 2 endpoints
    REQUIRE(rows.doublet.rows() == 5);
    REQUIRE(rows.doublet.cols() == 50);
    CHECK(rows.rhs.isZero(0.0));

    for (int r = 0; r < rows.nRows; ++r)
        CHECK(rows.doublet.row(r).cwiseAbs().sum() ==
              Catch::Approx(2.0).epsilon(1e-14));

    const int baseA = layout.network[0].doublet.edgeMidBase[1];
    const int baseB = 25 + layout.network[1].doublet.edgeMidBase[3];
    CHECK(rows.doublet(0, baseA) == Catch::Approx(1.0).epsilon(1e-14));
    CHECK(rows.doublet(0, baseB) == Catch::Approx(-1.0).epsilon(1e-14));
}

TEST_CASE("matching: source-only networks do not match doublets",
          "[matching]") {
    auto makeGrid = [](double y0, int kt, const char* name) {
        a502::Network net;
        net.netname = name;
        net.nm = 4;
        net.nn = 4;
        net.kt = kt;
        for (int n = 0; n < 4; ++n) {
            for (int m = 0; m < 4; ++m) {
                net.coordinates.push_back(static_cast<double>(m));
                net.coordinates.push_back(y0 + static_cast<double>(n));
                net.coordinates.push_back(0.0);
            }
        }
        return net;
    };

    const std::vector<a502::Network> networks = {
        makeGrid(0.0, 7, "src"), makeGrid(3.0, 12, "dbl")};
    const auto edges = a502::buildNetworkEdges(networks);
    REQUIRE(edges.size() == 8);
    const auto matches = a502::findFullEdgeMatches(edges, 1e-10);
    REQUIRE(matches.size() == 1);

    const auto layout =
        a502::buildGlobalSingularityLayout(networks);
    const auto rows =
        a502::buildEdgeDoubletMatchingRows(networks, layout, matches);
    REQUIRE(rows.valid);
    CHECK(rows.nRows == 0);
    CHECK(rows.doublet.rows() == 0);
    CHECK(rows.doublet.cols() == layout.totalDoublet);
    CHECK(rows.rhs.size() == 0);
}

TEST_CASE("global boundary system: active columns and matching rows",
          "[globalsystem]") {
    auto makeGrid = [](double y0, int kt, const char* name) {
        a502::Network net;
        net.netname = name;
        net.nm = 4;
        net.nn = 4;
        net.kt = kt;
        for (int n = 0; n < 4; ++n) {
            for (int m = 0; m < 4; ++m) {
                net.coordinates.push_back(static_cast<double>(m));
                net.coordinates.push_back(y0 + static_cast<double>(n));
                net.coordinates.push_back(0.0);
            }
        }
        return net;
    };

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const std::vector<a502::Network> networks = {
        makeGrid(0.0, 1, "a"), makeGrid(3.0, 1, "b")};
    const auto layout = a502::buildGlobalSingularityLayout(networks);
    REQUIRE(layout.activeSource == 18);
    REQUIRE(layout.activeDoublet == 50);

    const auto system = a502::buildGlobalNormalFluxSystem(
        networks, flow, layout, 1.0, 1.0, 1e-10);
    REQUIRE(system.valid);
    CHECK(system.nBcRows == 18);
    CHECK(system.nMatchingRows == 5);
    CHECK(system.nRows == 23);
    CHECK(system.nActiveSource == 18);
    CHECK(system.nActiveDoublet == 50);
    CHECK(system.nColumns == 68);
    REQUIRE(system.matrix.rows() == 23);
    REQUIRE(system.matrix.cols() == 68);
    REQUIRE(system.rhs.size() == 23);
    CHECK(system.matrix.allFinite());
    CHECK(system.rhs.allFinite());

    Eigen::FullPivLU<Eigen::MatrixXd> lu(system.matrix);
    CHECK(lu.rank() == system.nRows);

    const auto systemDiagnostics =
        a502::analyzeGlobalBoundarySystem(system);
    REQUIRE(systemDiagnostics.valid);
    CHECK(systemDiagnostics.rank == system.nRows);
    CHECK(systemDiagnostics.dependentRows == 0);
    CHECK(systemDiagnostics.zeroRows == 0);
    // With panel-centre rows only, some doublet parameters (mostly edge and
    // corner values) are still untouched by the assembled equations.
    CHECK(systemDiagnostics.zeroColumns > 0);
    CHECK(std::isfinite(systemDiagnostics.conditionNumber));

    const std::vector<a502::Network> doubletOnly = {
        makeGrid(0.0, 12, "d")};
    const auto doubletLayout =
        a502::buildGlobalSingularityLayout(doubletOnly);
    CHECK(doubletLayout.activeSource == 0);
    CHECK(doubletLayout.activeDoublet == 25);
    const auto doubletSystem = a502::buildGlobalNormalFluxSystem(
        doubletOnly, flow, doubletLayout, 1.0, 1.0, 1e-10);
    REQUIRE(doubletSystem.valid);
    CHECK(doubletSystem.nBcRows == 9);
    CHECK(doubletSystem.nMatchingRows == 0);
    CHECK(doubletSystem.nColumns == 25);
    REQUIRE(doubletSystem.matrix.rows() == 9);
    REQUIRE(doubletSystem.matrix.cols() == 25);

    const auto allPoints = a502::buildGlobalNormalFluxSystem(
        doubletOnly, flow, doubletLayout, 1.0, 1.0, 1e-10, true);
    REQUIRE(allPoints.valid);
    CHECK(allPoints.nBcRows == 25);
    CHECK(allPoints.nRows == 25);
    CHECK(allPoints.nColumns == 25);
    Eigen::FullPivLU<Eigen::MatrixXd> allLu(allPoints.matrix);
    CHECK(allLu.rank() == 25);

    const auto flowAlpha =
        a502::Compressibility::fromAngles(0.6, 5.0, 0.0);
    const auto solveSystem = a502::buildGlobalNormalFluxSystem(
        doubletOnly, flowAlpha, doubletLayout, 1.0, 1.0, 1e-10, true);
    REQUIRE(solveSystem.valid);
    const auto solution =
        a502::solveGlobalBoundarySystem(solveSystem);
    REQUIRE(solution.valid);
    CHECK(solution.unknowns.size() == 25);
    CHECK(solution.unknowns.allFinite());
    CHECK(solution.unknowns.norm() > 1e-6);
    CHECK(solution.residualNorm < 1e-3);

    const auto values =
        a502::expandGlobalSolution(solution, doubletLayout);
    REQUIRE(values.valid);
    REQUIRE(values.source.size() == 9);
    CHECK(values.source.isZero(0.0));
    REQUIRE(values.doublet.size() == 25);
    CHECK((values.doublet - solution.unknowns).norm() < 1e-12);

    a502::GlobalParameterValues zeroValues;
    zeroValues.source = Eigen::VectorXd::Zero(9);
    zeroValues.doublet = Eigen::VectorXd::Zero(25);
    zeroValues.valid = true;
    const auto uniformFlow = a502::evaluateSurfaceFlow(
        doubletOnly, flowAlpha, doubletLayout, zeroValues, true);
    REQUIRE(uniformFlow.valid);
    REQUIRE(uniformFlow.velocity.size() == 25);
    REQUIRE(uniformFlow.cp.size() == 25);
    for (std::size_t k = 0; k < uniformFlow.velocity.size(); ++k) {
        CHECK((uniformFlow.velocity[k] - flowAlpha.compd).norm() < 1e-12);
        CHECK(std::abs(uniformFlow.cp[k]) < 1e-12);
    }

    const auto surfaceFlow = a502::evaluateSurfaceFlow(
        doubletOnly, flowAlpha, doubletLayout, values, true);
    REQUIRE(surfaceFlow.valid);
    REQUIRE(surfaceFlow.velocity.size() == 25);
    REQUIRE(surfaceFlow.cp.size() == 25);
    bool nonZeroCp = false;
    for (std::size_t k = 0; k < surfaceFlow.cp.size(); ++k) {
        CHECK(surfaceFlow.velocity[k].allFinite());
        CHECK(std::isfinite(surfaceFlow.cp[k]));
        if (std::abs(surfaceFlow.cp[k]) > 1e-8)
            nonZeroCp = true;
    }
    CHECK(nonZeroCp);

    // Solver sanity check with an exactly solvable system.
    a502::GlobalBoundarySystem exact;
    exact.nRows = 3;
    exact.nColumns = 3;
    exact.nActiveSource = 3;
    exact.nActiveDoublet = 0;
    exact.matrix.resize(3, 3);
    exact.matrix << 2.0, 1.0, 0.0,
                    1.0, 3.0, 1.0,
                    0.0, 1.0, 2.0;
    Eigen::VectorXd wanted(3);
    wanted << 1.0, -2.0, 3.0;
    exact.rhs = exact.matrix * wanted;
    exact.valid = true;
    const auto exactSolution =
        a502::solveGlobalBoundarySystem(exact);
    REQUIRE(exactSolution.valid);
    CHECK(exactSolution.rank == 3);
    CHECK((exactSolution.unknowns - wanted).norm() < 1e-12);
    CHECK(exactSolution.residualNorm < 1e-12);
    CHECK(exactSolution.scaledResidualNorm < 1e-12);

    // The under-determined composite system is not solved directly.
    CHECK_FALSE(a502::solveGlobalBoundarySystem(system).valid);
}

TEST_CASE("global system: mixed source-only and doublet-only networks",
          "[globalsystem]") {
    auto makeGrid = [](double y0, int kt, const char* name) {
        a502::Network net;
        net.netname = name;
        net.nm = 4;
        net.nn = 4;
        net.kt = kt;
        for (int n = 0; n < 4; ++n) {
            for (int m = 0; m < 4; ++m) {
                net.coordinates.push_back(static_cast<double>(m));
                net.coordinates.push_back(y0 + static_cast<double>(n));
                net.coordinates.push_back(0.0);
            }
        }
        return net;
    };

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const std::vector<a502::Network> networks = {
        makeGrid(0.0, 7, "src"), makeGrid(10.0, 12, "dbl")};
    const auto layout =
        a502::buildGlobalSingularityLayout(networks);
    REQUIRE(layout.activeSource == 9);
    REQUIRE(layout.activeDoublet == 25);
    CHECK(layout.activeSourceOffset[0] == 0);
    CHECK(layout.activeSourceOffset[1] == -1);
    CHECK(layout.activeDoubletOffset[0] == -1);
    CHECK(layout.activeDoubletOffset[1] == 0);
    CHECK(layout.totalSource == 18);
    CHECK(layout.totalDoublet == 50);

    // Inactive per-network blocks of a global influence row stay zero.
    const auto panel = a502::buildPanelGeometry(
        networks[0].coordinates, 4, 4, 1, 1);
    const Eigen::Vector3d point = panel.cp[8] + Eigen::Vector3d(0, 0, 0.5);
    const auto influence =
        a502::evaluateGlobalNearFieldInfluence(
            networks, flow, point, layout);
    REQUIRE(influence.valid);
    CHECK_FALSE(influence.source
                    .block(0, layout.sourceOffset[0], 4,
                           layout.network[0].source.count)
                    .isZero(0.0));
    CHECK(influence.source
              .block(0, layout.sourceOffset[1], 4,
                     layout.network[1].source.count)
              .isZero(0.0));
    CHECK(influence.doublet
              .block(0, layout.doubletOffset[0], 4,
                     layout.network[0].doublet.count)
              .isZero(0.0));
    CHECK_FALSE(influence.doublet
                    .block(0, layout.doubletOffset[1], 4,
                           layout.network[1].doublet.count)
                    .isZero(0.0));

    const auto system = a502::buildGlobalNormalFluxSystem(
        networks, flow, layout, 1.0, 1.0, 1e-10);
    REQUIRE(system.valid);
    CHECK(system.nBcRows == 18);
    CHECK(system.nMatchingRows == 0);
    CHECK(system.nRows == 18);
    CHECK(system.nActiveSource == 9);
    CHECK(system.nActiveDoublet == 25);
    CHECK(system.nColumns == 34);
    REQUIRE(system.matrix.rows() == 18);
    REQUIRE(system.matrix.cols() == 34);
    REQUIRE(system.rhs.size() == 18);
    CHECK(system.matrix.allFinite());
    CHECK(system.rhs.allFinite());

    Eigen::FullPivLU<Eigen::MatrixXd> lu(system.matrix);
    CHECK(lu.rank() == system.nRows);

    // Adjacent source-only and doublet-only networks share an edge.  The
    // source-alone edge is excluded from the doublet matching search, so no
    // matching rows are generated.
    const std::vector<a502::Network> adjacent = {
        makeGrid(0.0, 7, "src"), makeGrid(3.0, 12, "dbl")};
    const auto adjacentLayout =
        a502::buildGlobalSingularityLayout(adjacent);
    const auto adjacentSystem = a502::buildGlobalNormalFluxSystem(
        adjacent, flow, adjacentLayout, 1.0, 1.0, 1e-10);
    REQUIRE(adjacentSystem.valid);
    CHECK(adjacentSystem.nMatchingRows == 0);
    CHECK(adjacentSystem.nBcRows == 18);
    CHECK(adjacentSystem.nColumns == 34);
}

TEST_CASE("global system: equilibration and diagnostics", "[globalsystem]") {
    a502::GlobalBoundarySystem scaled;
    scaled.nRows = 3;
    scaled.nColumns = 3;
    scaled.nActiveSource = 3;
    scaled.nActiveDoublet = 0;
    scaled.matrix.resize(3, 3);
    scaled.matrix << 1.0e3, 1.0e3, 1.0e-3,
                     1.0e3, -1.0e3, 2.0e-3,
                     1.0e3, 0.0, 3.0e-3;
    Eigen::VectorXd wanted(3);
    wanted << 1.0, -2.0, 3.0;
    scaled.rhs = scaled.matrix * wanted;
    scaled.valid = true;

    const auto scaling =
        a502::buildGlobalSystemScaling(scaled);
    REQUIRE(scaling.valid);
    REQUIRE(scaling.row.size() == 3);
    REQUIRE(scaling.column.size() == 3);
    for (int i = 0; i < 3; ++i)
        CHECK(scaling.row(i) == Catch::Approx(1.0e-3).epsilon(1e-12));
    CHECK(scaling.column(0) == Catch::Approx(1.0).epsilon(1e-12));
    CHECK(scaling.column(1) == Catch::Approx(1.0).epsilon(1e-12));
    CHECK(scaling.column(2) ==
          Catch::Approx(1.0 / 3.0e-6).epsilon(1e-12));

    const auto solution =
        a502::solveGlobalBoundarySystem(scaled);
    REQUIRE(solution.valid);
    CHECK(solution.rank == 3);
    CHECK((solution.unknowns - wanted).norm() < 1.0e-8);
    CHECK(solution.residualNorm < 1.0e-8);
    CHECK(solution.scaledResidualNorm <= solution.residualNorm);

    const auto diagnostics =
        a502::analyzeGlobalBoundarySystem(scaled);
    REQUIRE(diagnostics.valid);
    CHECK(diagnostics.nRows == 3);
    CHECK(diagnostics.nColumns == 3);
    CHECK(diagnostics.rank == 3);
    CHECK(diagnostics.dependentRows == 0);
    CHECK(diagnostics.zeroRows == 0);
    CHECK(diagnostics.zeroColumns == 0);
    CHECK(diagnostics.maxAbs == 1.0e3);
    CHECK(std::isfinite(diagnostics.conditionNumber));
    CHECK(diagnostics.conditionNumber > 1.0);

    // A dependent row and an all-zero row are reported, and the matrix is
    // recognized as rank deficient.
    a502::GlobalBoundarySystem deficient;
    deficient.nRows = 4;
    deficient.nColumns = 3;
    deficient.nActiveSource = 3;
    deficient.nActiveDoublet = 0;
    deficient.matrix.resize(4, 3);
    deficient.matrix << 1.0, 2.0, 0.0,
                        2.0, 4.0, 0.0,
                        0.0, 0.0, 1.0,
                        0.0, 0.0, 0.0;
    deficient.rhs = Eigen::VectorXd::Zero(4);
    deficient.valid = true;

    const auto deficientDiagnostics =
        a502::analyzeGlobalBoundarySystem(deficient);
    REQUIRE(deficientDiagnostics.valid);
    CHECK(deficientDiagnostics.rank == 2);
    CHECK(deficientDiagnostics.dependentRows == 2);
    CHECK(deficientDiagnostics.zeroRows == 1);
    CHECK(deficientDiagnostics.zeroColumns == 0);
    CHECK(std::isinf(deficientDiagnostics.conditionNumber));

    const auto deficientScaling =
        a502::buildGlobalSystemScaling(deficient);
    REQUIRE(deficientScaling.valid);
    CHECK(deficientScaling.row(3) == 1.0);
    CHECK(deficientScaling.row(0) == Catch::Approx(0.5).epsilon(1e-14));
    CHECK(deficientScaling.row(1) == Catch::Approx(0.25).epsilon(1e-14));
    CHECK_FALSE(a502::solveGlobalBoundarySystem(deficient).valid);

    a502::GlobalBoundarySystem empty;
    CHECK_FALSE(a502::buildGlobalSystemScaling(empty).valid);
    CHECK_FALSE(a502::analyzeGlobalBoundarySystem(empty).valid);
    CHECK_FALSE(a502::solveGlobalBoundarySystem(empty).valid);
}

TEST_CASE("control points: single square network", "[controlpoints]") {
    a502::Network net;
    net.netname = "plate";
    net.nm = 2;
    net.nn = 2;
    net.coordinates = {0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                       0.0, 1.0, 0.0, 1.0, 1.0, 0.0};

    const auto cps = a502::buildNetworkControlPoints(net);
    REQUIRE(cps.size() == 9);

    int centers = 0, edges = 0, corners = 0;
    for (const auto& cp : cps) {
        switch (cp.kind) {
        case a502::ControlPointKind::PanelCenter: ++centers; break;
        case a502::ControlPointKind::EdgeMidpoint: ++edges; break;
        case a502::ControlPointKind::Corner: ++corners; break;
        }
        CHECK(std::abs(cp.normal.norm() - 1.0) < 1e-12);
        CHECK(cp.normal.isApprox(Eigen::Vector3d(0, 0, -1)));
    }
    CHECK(centers == 1);
    CHECK(edges == 4);
    CHECK(corners == 4);

    // Panel center: P9 withdrawn 1% along the bisector towards P5/P8.
    const auto center = std::find_if(
        cps.begin(), cps.end(), [](const a502::ControlPoint& c) {
            return c.kind == a502::ControlPointKind::PanelCenter;
        });
    REQUIRE(center != cps.end());
    CHECK(center->position.isApprox(Eigen::Vector3d(0.4975, 0.4975, 0.0),
                                    1e-14));

    // Corner P1 = (0,0), withdrawn 25% towards P5/P8.
    const auto corner = std::find_if(
        cps.begin(), cps.end(), [](const a502::ControlPoint& c) {
            return c.kind == a502::ControlPointKind::Corner;
        });
    REQUIRE(corner != cps.end());
    CHECK(corner->position.isApprox(Eigen::Vector3d(0.0625, 0.0625, 0.0),
                                    1e-14));
}

TEST_CASE("control points: sample networks", "[controlpoints]") {
    a502::AeroModel model;
    model.parse(samplePath("swb.inp").string());
    REQUIRE(model.solve());

    std::size_t wingaCenters = 0, wingaEdges = 0, wingaCorners = 0;
    for (const auto& cp : model.controlPoints) {
        if (cp.network != 0)
            continue;
        // winga: nm=11, nn=3 -> 20 centers, 24 edge midpoints, 4 corners.
        CHECK(std::abs(cp.normal.norm() - 1.0) < 1e-10);
        switch (cp.kind) {
        case a502::ControlPointKind::PanelCenter: ++wingaCenters; break;
        case a502::ControlPointKind::EdgeMidpoint: ++wingaEdges; break;
        case a502::ControlPointKind::Corner: ++wingaCorners; break;
        }
    }
    CHECK(wingaCenters == 20);
    CHECK(wingaEdges == 24);
    CHECK(wingaCorners == 4);

    // No coincident control points within a network.
    for (std::size_t i = 0; i < model.controlPoints.size(); ++i) {
        for (std::size_t j = i + 1; j < model.controlPoints.size(); ++j) {
            const auto& a = model.controlPoints[i];
            const auto& b = model.controlPoints[j];
            if (a.network != b.network)
                continue;
            CHECK((a.position - b.position).norm() > 1e-9);
        }
    }

    // bodyb contains collapsed rows/points; the generator must not crash and
    // must not emit coincident points.
    std::size_t bodybCount = 0;
    for (const auto& cp : model.controlPoints)
        if (cp.network == 4)
            ++bodybCount;
    CHECK(bodybCount > 0);
}

namespace {

double expectedPanelSource(int cpIndex, const a502::PanelSourceValues& v) {
    switch (cpIndex) {
    case 0: return v(0);
    case 1: return v(1);
    case 2: return v(2);
    case 3: return v(3);
    case 4: return 0.5 * (v(0) + v(1));
    case 5: return 0.5 * (v(1) + v(2));
    case 6: return 0.5 * (v(2) + v(3));
    case 7: return 0.5 * (v(3) + v(0));
    case 8: return v(4);
    default: return 0.0;
    }
}

double evaluateSourceAtGlobal(const a502::PanelGeometry& g, int s,
                              const Eigen::Vector3d& p,
                              const a502::PanelSourceValues& v) {
    const auto& sp = g.subpanels[static_cast<std::size_t>(s)];
    const Eigen::Matrix3d& frame =
        sp.scaledValid ? sp.scaledToLocal : sp.toLocal;
    const auto splines = a502::buildSourceSubpanelSplines(g);
    const auto& spline = splines[static_cast<std::size_t>(s)];
    const Eigen::Vector3d local = frame * (p - spline.origin);
    return a502::evaluateSource(spline, local.head<2>(), v);
}

} // namespace

TEST_CASE("source subpanel splines: vertex values and continuity",
          "[splines][source]") {
    const std::vector<double> coords = {
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 1.0, 1.0, 0.0};
    const a502::PanelGeometry g =
        a502::buildPanelGeometry(coords, 2, 2, 0, 0);
    REQUIRE_FALSE(g.degenerate);

    a502::PanelSourceValues vals;
    vals << 1.0, 2.0, 3.0, 4.0, 5.0;

    const auto splines = a502::buildSourceSubpanelSplines(g);
    for (int s = 0; s < 8; ++s) {
        INFO("subpanel " << s);
        const auto& spline = splines[static_cast<std::size_t>(s)];
        REQUIRE(spline.valid);
        const auto& sp = g.subpanels[static_cast<std::size_t>(s)];
        if (s >= 4)
            CHECK((spline.origin - g.cp[8]).norm() < 1e-14);
        else
            CHECK((spline.origin - sp.vertex[0]).norm() < 1e-14);
        for (int k = 0; k < 3; ++k) {
            const double got =
                evaluateSourceAtGlobal(g, s, sp.vertex[static_cast<std::size_t>(k)], vals);
            const double want =
                expectedPanelSource(sp.cpIndex[static_cast<std::size_t>(k)], vals);
            CHECK(got == Catch::Approx(want).epsilon(1e-12));
        }
    }

    // Shared edge P5-P8 between subpanels 0 and 4 (1 and 5 in A502 numbering).
    const Eigen::Vector3d p58 = 0.5 * (g.cp[4] + g.cp[7]);
    CHECK(evaluateSourceAtGlobal(g, 0, p58, vals) ==
          Catch::Approx(evaluateSourceAtGlobal(g, 4, p58, vals)).epsilon(1e-12));

    // Shared edge P5-P9 between subpanels 4 and 5 (5 and 6 in A502 numbering).
    const Eigen::Vector3d p59 = 0.5 * (g.cp[4] + g.cp[8]);
    CHECK(evaluateSourceAtGlobal(g, 4, p59, vals) ==
          Catch::Approx(evaluateSourceAtGlobal(g, 5, p59, vals)).epsilon(1e-12));
}

TEST_CASE("source subpanel splines: linear reproduction", "[splines][source]") {
    const std::vector<double> coords = {
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 1.0, 1.0, 0.0};
    const a502::PanelGeometry g =
        a502::buildPanelGeometry(coords, 2, 2, 0, 0);
    REQUIRE_FALSE(g.degenerate);

    auto f = [](const Eigen::Vector3d& p) {
        return 1.0 + 2.0 * p.x() + 3.0 * p.y();
    };
    a502::PanelSourceValues vals;
    vals << f(g.cp[0]), f(g.cp[1]), f(g.cp[2]), f(g.cp[3]), f(g.cp[8]);

    for (int s = 0; s < 8; ++s) {
        const auto& sp = g.subpanels[static_cast<std::size_t>(s)];
        const Eigen::Vector3d c =
            (sp.vertex[0] + sp.vertex[1] + sp.vertex[2]) / 3.0;
        INFO("subpanel " << s);
        CHECK(evaluateSourceAtGlobal(g, s, c, vals) ==
              Catch::Approx(f(c)).epsilon(1e-12));
    }
}

TEST_CASE("source subpanel splines: sample panel", "[splines][source]") {
    a502::AeroModel model;
    model.parse(samplePath("swb.inp").string());
    REQUIRE(model.solve());

    const auto& panel = model.panels.front().geom;
    const auto splines = a502::buildSourceSubpanelSplines(panel);

    a502::PanelSourceValues vals;
    vals << 0.1, -0.2, 0.3, -0.4, 0.5;

    for (int s = 0; s < 8; ++s) {
        INFO("subpanel " << s);
        REQUIRE(splines[static_cast<std::size_t>(s)].valid);
        const auto& sp = panel.subpanels[static_cast<std::size_t>(s)];
        for (int k = 0; k < 3; ++k) {
            CHECK(evaluateSourceAtGlobal(panel, s,
                                         sp.vertex[static_cast<std::size_t>(k)],
                                         vals) ==
                  Catch::Approx(expectedPanelSource(
                      sp.cpIndex[static_cast<std::size_t>(k)], vals))
                      .epsilon(1e-10));
        }
    }
}

TEST_CASE("biquadratic basis functions", "[splines][bqbfun]") {
    struct Pt { int index; double s; double t; };
    const Pt pts[9] = {{0, 1, 1},  {1, -1, 1},  {2, -1, -1},
                       {3, 1, -1}, {4, 0, 1},   {5, -1, 0},
                       {6, 0, -1}, {7, 1, 0},   {8, 0, 0}};

    for (const auto& pt : pts) {
        const auto b = a502::bqbfun(pt.s, pt.t);
        INFO("P" << (pt.index + 1) << " at (" << pt.s << "," << pt.t << ")");
        for (int i = 0; i < 9; ++i) {
            if (i == pt.index)
                CHECK(b.phi(i) == Catch::Approx(1.0).epsilon(1e-14));
            else
                CHECK(std::abs(b.phi(i)) < 1e-14);
        }
        CHECK(b.phi.sum() == Catch::Approx(1.0).epsilon(1e-14));
        CHECK(std::abs(b.phis.sum()) < 1e-12);
        CHECK(std::abs(b.phit.sum()) < 1e-12);
    }
}

TEST_CASE("triangle interpolation: linear map", "[splines][tcof]") {
    Eigen::Matrix3d p;
    p << 0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         0.0, 0.0, 0.0; // unused third row

    const auto map = a502::buildLinearTriangleMap(p);
    REQUIRE(map.valid);

    Eigen::Vector3d values(1.0, 5.0, 7.0);
    const Eigen::Vector3d c = map.b * values;
    CHECK(c(0) == Catch::Approx(1.0).epsilon(1e-14));
    CHECK(c(1) == Catch::Approx(4.0).epsilon(1e-14));
    CHECK(c(2) == Catch::Approx(6.0).epsilon(1e-14));
}

TEST_CASE("triangle interpolation: quadratic map", "[splines][tcof]") {
    Eigen::Matrix3d p;
    p << 0.0, 1.0, 0.0,
         0.0, 0.0, 1.0,
         0.0, 0.0, 0.0;

    const auto map = a502::buildQuadraticTriangleMap(p);
    REQUIRE(map.valid);

    auto f = [](double x, double y) {
        return 1.0 + 2.0 * x + 3.0 * y + 4.0 * x * x / 2.0 +
               5.0 * x * y + 6.0 * y * y / 2.0;
    };

    const double v0 = f(0.0, 0.0);
    const double v1 = f(1.0, 0.0);
    const double v2 = f(0.0, 1.0);
    const double m12 = f(0.5, 0.5);
    const double m20 = f(0.0, 0.5);
    const double m01 = f(0.5, 0.0);
    const double l1 = 2.0 * m12 - 0.5 * (v1 + v2);
    const double l2 = 2.0 * m20 - 0.5 * (v2 + v0);
    const double l3 = 2.0 * m01 - 0.5 * (v0 + v1);

    Eigen::Matrix<double, 6, 1> rhs;
    rhs << v0, v1, v2, l1, l2, l3;
    const Eigen::Matrix<double, 6, 1> c = map.c * rhs;
    const Eigen::Matrix<double, 6, 1> want =
        (Eigen::Matrix<double, 6, 1>() << 1, 2, 3, 4, 5, 6).finished();
    CHECK((c - want).norm() < 1e-12);

    const double x = 0.25;
    const double y = 0.75;
    const double evaluated = c(0) + c(1) * x + c(2) * y +
                             c(3) * x * x / 2.0 + c(4) * x * y +
                             c(5) * y * y / 2.0;
    CHECK(evaluated == Catch::Approx(f(x, y)).epsilon(1e-12));

    // A non-right triangle exercises the general inverse.
    Eigen::Matrix3d q;
    q << 0.0, 1.0, 0.3,
         0.0, 0.0, 0.8,
         0.0, 0.0, 0.0;
    const auto map2 = a502::buildQuadraticTriangleMap(q);
    REQUIRE(map2.valid);

    auto f2 = [](double x, double y) {
        return -2.0 + 0.5 * x - 1.5 * y + 2.0 * x * x / 2.0 -
               3.0 * x * y + 4.0 * y * y / 2.0;
    };
    const double w0 = f2(q(0, 0), q(1, 0));
    const double w1 = f2(q(0, 1), q(1, 1));
    const double w2 = f2(q(0, 2), q(1, 2));
    const auto mid = [&](int a, int b) -> Eigen::Vector2d {
        return 0.5 * (Eigen::Vector2d(q(0, a), q(1, a)) +
                      Eigen::Vector2d(q(0, b), q(1, b)));
    };
    const double n12 = f2(mid(1, 2).x(), mid(1, 2).y());
    const double n20 = f2(mid(2, 0).x(), mid(2, 0).y());
    const double n01 = f2(mid(0, 1).x(), mid(0, 1).y());
    Eigen::Matrix<double, 6, 1> rhs2;
    rhs2 << w0, w1, w2,
        2.0 * n12 - 0.5 * (w1 + w2),
        2.0 * n20 - 0.5 * (w2 + w0),
        2.0 * n01 - 0.5 * (w0 + w1);
    const Eigen::Matrix<double, 6, 1> c2 = map2.c * rhs2;
    const Eigen::Matrix<double, 6, 1> want2 =
        (Eigen::Matrix<double, 6, 1>() << -2.0, 0.5, -1.5, 2.0, -3.0, 4.0)
            .finished();
    CHECK((c2 - want2).norm() < 1e-12);
}

TEST_CASE("panel parameters: canonical points", "[splines][nrpthp]") {
    const std::vector<double> coords = {
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 1.0, 1.0, 0.0};
    const a502::PanelGeometry g =
        a502::buildPanelGeometry(coords, 2, 2, 0, 0);
    REQUIRE_FALSE(g.degenerate);

    const std::array<Eigen::Vector3d, 4> corners = {
        g.cp[0], g.cp[1], g.cp[2], g.cp[3]};

    // Canonical (s,t) values (same convention as bqbfun).
    const double expected[9][2] = {
        {1, 1},   {-1, 1},  {-1, -1}, {1, -1}, {0, 1},
        {-1, 0},  {0, -1},  {1, 0},   {0, 0}};

    for (int i = 0; i < 9; ++i) {
        const auto par = a502::closestPanelParameters(corners, g.cp[i]);
        INFO("P" << (i + 1));
        REQUIRE(par.valid);
        CHECK(par.s == Catch::Approx(expected[i][0]).margin(1e-10));
        CHECK(par.t == Catch::Approx(expected[i][1]).margin(1e-10));
    }

    // The same property must hold for a non-planar sample panel.
    a502::AeroModel model;
    model.parse(samplePath("swb.inp").string());
    REQUIRE(model.solve());
    const auto& gp = model.panels.front().geom;
    const std::array<Eigen::Vector3d, 4> c2 = {
        gp.cp[0], gp.cp[1], gp.cp[2], gp.cp[3]};
    for (int i = 0; i < 9; ++i) {
        const auto par = a502::closestPanelParameters(c2, gp.cp[i]);
        INFO("sample P" << (i + 1));
        REQUIRE(par.valid);
        CHECK(par.s == Catch::Approx(expected[i][0]).margin(1e-10));
        CHECK(par.t == Catch::Approx(expected[i][1]).margin(1e-10));
    }
}

namespace {

double evaluateDoubletAtGlobal(
    const a502::PanelGeometry& g,
    const std::array<a502::DoubletSubpanelSpline, 8>& splines, int s,
    const Eigen::Vector3d& p, const a502::PanelDoubletValues& v) {
    const auto& sp = g.subpanels[static_cast<std::size_t>(s)];
    const Eigen::Matrix3d& frame =
        sp.scaledValid ? sp.scaledToLocal : sp.toLocal;
    const auto& spline = splines[static_cast<std::size_t>(s)];
    const Eigen::Vector3d local = frame * (p - spline.origin);
    return a502::evaluateDoublet(spline, local.head<2>(), v);
}

} // namespace

TEST_CASE("doublet subpanel splines: constant and quadratic reproduction",
          "[splines][doublet]") {
    const std::vector<double> coords = {
        0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
        0.0, 1.0, 0.0, 1.0, 1.0, 0.0};
    const a502::PanelGeometry g =
        a502::buildPanelGeometry(coords, 2, 2, 0, 0);
    REQUIRE_FALSE(g.degenerate);

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const auto splines = a502::buildDoubletSubpanelSplines(g, flow);
    for (int s = 0; s < 8; ++s)
        REQUIRE(splines[static_cast<std::size_t>(s)].valid);

    // Constant distribution.
    a502::PanelDoubletValues ones;
    ones.setOnes();
    for (int s = 0; s < 8; ++s) {
        const auto c = a502::doubletTaylorCoefficients(
            splines[static_cast<std::size_t>(s)], ones);
        CHECK(c(0) == Catch::Approx(1.0).epsilon(1e-12));
        for (int k = 1; k < 6; ++k)
            CHECK(std::abs(c(k)) < 1e-12);
    }

    // Quadratic distribution in the inner-panel frame.
    const Eigen::Matrix3d& innerFrame = g.subpanels[4].scaledToLocal;
    auto innerLocal = [&](const Eigen::Vector3d& p) {
        return (innerFrame * (p - g.cp[8])).head<2>();
    };
    auto mu = [](const Eigen::Vector2d& q) {
        return 1.0 + 2.0 * q.x() + 3.0 * q.y() +
               4.0 * q.x() * q.x() / 2.0 +
               5.0 * q.x() * q.y() +
               6.0 * q.y() * q.y() / 2.0;
    };

    a502::PanelDoubletValues vals;
    for (int i = 0; i < 9; ++i)
        vals(i) = mu(innerLocal(g.cp[i]));

    for (int s = 0; s < 8; ++s) {
        const auto& sp = g.subpanels[static_cast<std::size_t>(s)];
        for (int k = 0; k < 3; ++k) {
            const int ci = sp.cpIndex[static_cast<std::size_t>(k)];
            INFO("subpanel " << s << " vertex " << k);
            CHECK(evaluateDoubletAtGlobal(g, splines, s,
                                          g.cp[ci], vals) ==
                  Catch::Approx(vals(ci)).epsilon(1e-10));
        }
        const Eigen::Vector3d centroid =
            (sp.vertex[0] + sp.vertex[1] + sp.vertex[2]) / 3.0;
        INFO("subpanel " << s << " centroid");
        CHECK(evaluateDoubletAtGlobal(g, splines, s, centroid, vals) ==
              Catch::Approx(mu(innerLocal(centroid))).epsilon(1e-10));
    }

    // Continuity on shared edges.
    const Eigen::Vector3d p58 = 0.5 * (g.cp[4] + g.cp[7]);
    CHECK(evaluateDoubletAtGlobal(g, splines, 0, p58, vals) ==
          Catch::Approx(evaluateDoubletAtGlobal(g, splines, 4, p58, vals))
              .epsilon(1e-10));
    const Eigen::Vector3d p59 = 0.5 * (g.cp[4] + g.cp[8]);
    CHECK(evaluateDoubletAtGlobal(g, splines, 4, p59, vals) ==
          Catch::Approx(evaluateDoubletAtGlobal(g, splines, 5, p59, vals))
              .epsilon(1e-10));
}

TEST_CASE("doublet subpanel splines: sample panel", "[splines][doublet]") {
    a502::AeroModel model;
    model.parse(samplePath("swb.inp").string());
    REQUIRE(model.solve());

    const auto& g = model.panels.front().geom;
    const auto splines = a502::buildDoubletSubpanelSplines(g, model.flow);
    for (int s = 0; s < 8; ++s) {
        INFO("subpanel " << s);
        REQUIRE(splines[static_cast<std::size_t>(s)].valid);
    }

    // Constant and a quadratic in the inner frame.
    const Eigen::Matrix3d& innerFrame = g.subpanels[4].scaledToLocal;
    auto innerLocal = [&](const Eigen::Vector3d& p) {
        return (innerFrame * (p - g.cp[8])).head<2>();
    };
    auto mu = [](const Eigen::Vector2d& q) {
        return -1.0 + 0.5 * q.x() - 2.0 * q.y() +
               3.0 * q.x() * q.x() / 2.0 -
               1.5 * q.x() * q.y() +
               2.5 * q.y() * q.y() / 2.0;
    };

    a502::PanelDoubletValues vals;
    for (int i = 0; i < 9; ++i)
        vals(i) = mu(innerLocal(g.cp[i]));

    for (int s = 0; s < 8; ++s) {
        const auto& sp = g.subpanels[static_cast<std::size_t>(s)];
        for (int k = 0; k < 3; ++k) {
            const int ci = sp.cpIndex[static_cast<std::size_t>(k)];
            CHECK(evaluateDoubletAtGlobal(g, splines, s, g.cp[ci], vals) ==
                  Catch::Approx(vals(ci)).epsilon(1e-9));
        }
        const Eigen::Vector3d centroid =
            (sp.vertex[0] + sp.vertex[1] + sp.vertex[2]) / 3.0;
        CHECK(evaluateDoubletAtGlobal(g, splines, s, centroid, vals) ==
              Catch::Approx(mu(innerLocal(centroid))).epsilon(1e-9));
    }
}

TEST_CASE("source outer spline: interior panel", "[splines][outer]") {
    // 4x4 grid of points -> 3x3 panels of unit size in the z=0 plane.
    a502::Network net;
    net.netname = "grid";
    net.nm = 4;
    net.nn = 4;
    net.coordinates.reserve(4 * 4 * 3);
    for (int n = 0; n < 4; ++n) {
        for (int m = 0; m < 4; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto spline = a502::buildSourceOuterSpline(net, 1, 1);
    REQUIRE(spline.valid);
    for (int k = 0; k < 9; ++k)
        CHECK(spline.neighborPanel[static_cast<std::size_t>(k)] >= 0);

    // Row sums must be one (constant source reproduction).
    for (int r = 0; r < 5; ++r)
        CHECK(spline.coefficients.row(r).sum() ==
              Catch::Approx(1.0).epsilon(1e-10));

    // A linear source function must be reproduced exactly at P1..P4,P9.
    auto f = [](double x, double y) { return 1.0 + 2.0 * x + 3.0 * y; };

    Eigen::Matrix<double, 9, 1> neighborValues;
    int k = 0;
    for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
            const double x = (1 + di) + 0.5;
            const double y = (1 + dj) + 0.5;
            neighborValues(k++) = f(x, y);
        }
    }

    const Eigen::Matrix<double, 5, 1> got = spline.coefficients * neighborValues;

    const a502::PanelGeometry panel =
        a502::buildPanelGeometry(net.coordinates, net.nm, net.nn, 1, 1);
    REQUIRE_FALSE(panel.degenerate);
    for (int r = 0; r < 4; ++r)
        CHECK(got(r) ==
              Catch::Approx(f(panel.cp[r].x(), panel.cp[r].y()))
                  .epsilon(1e-10));
    CHECK(got(4) ==
          Catch::Approx(f(panel.cp[8].x(), panel.cp[8].y())).epsilon(1e-10));

    // Edge and corner panels use the clipped 3x3 stencil of sing().
    const auto edge = a502::buildSourceOuterSpline(net, 0, 1);
    REQUIRE(edge.valid);
    int edgeMissing = 0;
    for (int q = 0; q < 9; ++q)
        if (edge.neighborPanel[static_cast<std::size_t>(q)] < 0)
            ++edgeMissing;
    CHECK(edgeMissing == 3);
    for (int r = 0; r < 5; ++r)
        CHECK(edge.coefficients.row(r).sum() ==
              Catch::Approx(1.0).epsilon(1e-10));

    Eigen::Matrix<double, 9, 1> edgeValues =
        Eigen::Matrix<double, 9, 1>::Zero();
    k = 0;
    for (int dj = -1; dj <= 1; ++dj) {
        for (int di = -1; di <= 1; ++di) {
            const int ni = 0 + di;
            const int nj = 1 + dj;
            if (ni >= 0 && ni < 3 && nj >= 0 && nj < 3)
                edgeValues(k) = f(ni + 0.5, nj + 0.5);
            ++k;
        }
    }
    const auto edgeGeom =
        a502::buildPanelGeometry(net.coordinates, net.nm, net.nn, 0, 1);
    const Eigen::Matrix<double, 5, 1> edgeGot = edge.coefficients * edgeValues;
    for (int r = 0; r < 4; ++r)
        CHECK(edgeGot(r) ==
              Catch::Approx(f(edgeGeom.cp[r].x(), edgeGeom.cp[r].y()))
                  .epsilon(1e-9));
    CHECK(edgeGot(4) ==
          Catch::Approx(f(edgeGeom.cp[8].x(), edgeGeom.cp[8].y()))
              .epsilon(1e-9));

    const auto corner = a502::buildSourceOuterSpline(net, 0, 0);
    REQUIRE(corner.valid);
    int cornerMissing = 0;
    for (int q = 0; q < 9; ++q)
        if (corner.neighborPanel[static_cast<std::size_t>(q)] < 0)
            ++cornerMissing;
    CHECK(cornerMissing == 5);
}

TEST_CASE("source outer spline: sample interior panel", "[splines][outer]") {
    a502::AeroModel model;
    model.parse(samplePath("btac.inp").string());
    REQUIRE(model.networks.size() > 2);
    // akl is an 11x11 network with interior panels.
    const auto& net = model.networks[2];
    REQUIRE(net.nm == 11);
    REQUIRE(net.nn == 11);

    const auto spline = a502::buildSourceOuterSpline(net, 1, 1);
    REQUIRE(spline.valid);
    for (int r = 0; r < 5; ++r)
        CHECK(spline.coefficients.row(r).sum() ==
              Catch::Approx(1.0).epsilon(1e-9));
}

TEST_CASE("doublet outer spline: deep interior panel", "[splines][outer]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    net.coordinates.reserve(8 * 8 * 3);
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const auto spline = a502::buildDoubletOuterSpline(net, flow, 3, 3);
    REQUIRE(spline.valid);

    // All 25 neighbour panel indices must be present and distinct.
    std::array<int, 25> idx = spline.neighborPanel;
    std::sort(idx.begin(), idx.end());
    for (int k = 0; k < 25; ++k) {
        CHECK(idx[static_cast<std::size_t>(k)] >= 0);
        if (k > 0)
            CHECK(idx[static_cast<std::size_t>(k)] !=
                  idx[static_cast<std::size_t>(k - 1)]);
    }

    // Constant doublet strength.
    Eigen::Matrix<double, 25, 1> ones =
        Eigen::Matrix<double, 25, 1>::Ones();
    const Eigen::Matrix<double, 9, 1> constant =
        spline.coefficients * ones;
    for (int r = 0; r < 9; ++r)
        CHECK(constant(r) == Catch::Approx(1.0).epsilon(1e-10));

    // Quadratic doublet distribution.
    auto mu = [](double x, double y) {
        return 1.0 + 2.0 * x + 3.0 * y + 4.0 * x * x / 2.0 +
               5.0 * x * y + 6.0 * y * y / 2.0;
    };
    Eigen::Matrix<double, 25, 1> values;
    for (int k = 0; k < 25; ++k) {
        const int p = spline.neighborPanel[static_cast<std::size_t>(k)];
        const int pi = p % (net.nm - 1);
        const int pj = p / (net.nm - 1);
        const auto g = a502::buildPanelGeometry(
            net.coordinates, net.nm, net.nn, pi, pj);
        REQUIRE_FALSE(g.degenerate);
        values(k) = mu(g.cp[8].x(), g.cp[8].y());
    }

    const a502::PanelGeometry panel =
        a502::buildPanelGeometry(net.coordinates, net.nm, net.nn, 3, 3);
    REQUIRE_FALSE(panel.degenerate);
    const Eigen::Matrix<double, 9, 1> got = spline.coefficients * values;
    for (int r = 0; r < 9; ++r)
        CHECK(got(r) ==
              Catch::Approx(mu(panel.cp[r].x(), panel.cp[r].y()))
                  .epsilon(1e-8));
}

TEST_CASE("doublet outer spline: sample deep interior panel",
          "[splines][outer]") {
    a502::AeroModel model;
    model.parse(samplePath("btac.inp").string());
    REQUIRE(model.solve());
    REQUIRE(model.networks.size() > 2);
    const auto& net = model.networks[2]; // akl, 11x11
    REQUIRE(net.nm == 11);
    REQUIRE(net.nn == 11);

    const auto spline = a502::buildDoubletOuterSpline(net, model.flow, 3, 3);
    REQUIRE(spline.valid);

    Eigen::Matrix<double, 25, 1> ones =
        Eigen::Matrix<double, 25, 1>::Ones();
    const Eigen::Matrix<double, 9, 1> constant =
        spline.coefficients * ones;
    for (int r = 0; r < 9; ++r)
        CHECK(constant(r) == Catch::Approx(1.0).epsilon(1e-8));
}

TEST_CASE("doublet network parameters: staggered parameter grid",
          "[splines][outer]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto parameters = a502::buildNetworkDoubletParameters(net);
    REQUIRE(parameters.nPanelM == 7);
    REQUIRE(parameters.nPanelN == 7);
    REQUIRE(parameters.nCenters == 49);
    // (nm+1)*(nn+1) = 81 parameters: 49 centres, 4 corners, 28 edge points.
    REQUIRE(parameters.count == 81);
    REQUIRE(parameters.location.size() == 81);

    CHECK(parameters.cornerId[0] == 49);
    CHECK(parameters.cornerId[1] == 50);
    CHECK(parameters.cornerId[2] == 51);
    CHECK(parameters.cornerId[3] == 52);
    CHECK(parameters.edgeMidBase[0] == 53); // top row
    CHECK(parameters.edgeMidBase[1] == 60); // right column
    CHECK(parameters.edgeMidBase[2] == 67); // bottom row
    CHECK(parameters.edgeMidBase[3] == 74); // left column

    // Parameter-grid mapping: 1-based (m,n) in (nm+1)x(nn+1).
    CHECK(parameters.parameterId(2, 2) == parameters.centerId(0, 0));
    CHECK(parameters.parameterId(8, 8) == parameters.centerId(6, 6));
    CHECK(parameters.parameterId(1, 1) == parameters.cornerId[0]);
    CHECK(parameters.parameterId(1, 9) == parameters.cornerId[1]);
    CHECK(parameters.parameterId(9, 9) == parameters.cornerId[2]);
    CHECK(parameters.parameterId(9, 1) == parameters.cornerId[3]);
    CHECK(parameters.parameterId(1, 2) == parameters.edgeMidBase[0]);
    CHECK(parameters.parameterId(2, 9) == parameters.edgeMidBase[1]);
    CHECK(parameters.parameterId(9, 2) == parameters.edgeMidBase[2]);
    CHECK(parameters.parameterId(2, 1) == parameters.edgeMidBase[3]);
    CHECK(parameters.parameterId(0, 1) == -1);
    CHECK(parameters.parameterId(10, 1) == -1);

    // Locations follow gcpcal(): panel centres, edge midpoints, corners.
    auto close = [](const Eigen::Vector3d& a, double x, double y) {
        return a.x() == Catch::Approx(x).epsilon(1e-14) &&
               a.y() == Catch::Approx(y).epsilon(1e-14) &&
               a.z() == Catch::Approx(0.0).epsilon(1e-14);
    };
    CHECK(close(parameters.location[parameters.centerId(0, 0)], 0.5, 0.5));
    CHECK(close(parameters.location[parameters.centerId(6, 6)], 6.5, 6.5));
    CHECK(close(parameters.location[parameters.cornerId[0]], 0.0, 0.0));
    CHECK(close(parameters.location[parameters.cornerId[2]], 7.0, 7.0));
    CHECK(close(parameters.location[parameters.edgeMidBase[0]], 0.0, 0.5));
    CHECK(close(parameters.location[parameters.edgeMidBase[1]], 0.5, 7.0));
    CHECK(close(parameters.location[parameters.edgeMidBase[2]], 7.0, 0.5));
    CHECK(close(parameters.location[parameters.edgeMidBase[3]], 0.5, 0.0));
}

TEST_CASE("doublet panel spline: edge panel", "[splines][outer]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const auto parameters = a502::buildNetworkDoubletParameters(net);
    REQUIRE(parameters.count == 81);

    // Panel (0,3) lies on the first network row, away from the corners.
    const auto spline =
        a502::buildDoubletPanelSpline(net, flow, 0, 3, parameters);
    REQUIRE(spline.valid);

    // P5 is the midpoint of the first panel edge and is an explicit edge
    // parameter.
    const int p5Id = parameters.edgeMidBase[0] + 3;
    REQUIRE(spline.rows[4].index.size() == 1);
    CHECK(spline.rows[4].index[0] == p5Id);
    CHECK(spline.rows[4].coefficient[0] == Catch::Approx(1.0).epsilon(1e-14));

    // P9 is the panel-centre parameter.
    const int p9Id = parameters.centerId(0, 3);
    REQUIRE(spline.rows[8].index.size() == 1);
    CHECK(spline.rows[8].index[0] == p9Id);
    CHECK(spline.rows[8].coefficient[0] == Catch::Approx(1.0).epsilon(1e-14));

    auto mu = [](double x, double y) {
        return 1.0 + 2.0 * x + 3.0 * y + 4.0 * x * x / 2.0 +
               5.0 * x * y + 6.0 * y * y / 2.0;
    };

    const a502::PanelGeometry panel =
        a502::buildPanelGeometry(net.coordinates, net.nm, net.nn, 0, 3);
    REQUIRE_FALSE(panel.degenerate);

    for (int r = 0; r < 9; ++r) {
        double got = 0.0;
        double rowSum = 0.0;
        for (std::size_t k = 0; k < spline.rows[static_cast<std::size_t>(r)]
                                         .index.size();
             ++k) {
            const int id =
                spline.rows[static_cast<std::size_t>(r)].index[k];
            const double c =
                spline.rows[static_cast<std::size_t>(r)].coefficient[k];
            const Eigen::Vector3d& p =
                parameters.location[static_cast<std::size_t>(id)];
            got += c * mu(p.x(), p.y());
            rowSum += c;
        }
        INFO("canonical point " << r);
        CHECK(rowSum == Catch::Approx(1.0).epsilon(1e-9));
        CHECK(got == Catch::Approx(mu(panel.cp[r].x(), panel.cp[r].y()))
                       .epsilon(1e-8));
    }
}

TEST_CASE("doublet panel spline: all synthetic panels", "[splines][outer]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const auto parameters = a502::buildNetworkDoubletParameters(net);

    auto mu = [](double x, double y) {
        return 1.0 + 2.0 * x + 3.0 * y + 4.0 * x * x / 2.0 +
               5.0 * x * y + 6.0 * y * y / 2.0;
    };

    for (int panelJ = 0; panelJ < 7; ++panelJ) {
        for (int panelI = 0; panelI < 7; ++panelI) {
            INFO("panel " << panelI << "," << panelJ);
            const auto spline = a502::buildDoubletPanelSpline(
                net, flow, panelI, panelJ, parameters);
            REQUIRE(spline.valid);

            const a502::PanelGeometry panel = a502::buildPanelGeometry(
                net.coordinates, net.nm, net.nn, panelI, panelJ);
            REQUIRE_FALSE(panel.degenerate);

            for (int r = 0; r < 9; ++r) {
                double got = 0.0;
                for (std::size_t k = 0;
                     k < spline.rows[static_cast<std::size_t>(r)]
                             .index.size();
                     ++k) {
                    const int id =
                        spline.rows[static_cast<std::size_t>(r)].index[k];
                    const double c =
                        spline.rows[static_cast<std::size_t>(r)]
                            .coefficient[k];
                    const Eigen::Vector3d& p =
                        parameters.location[static_cast<std::size_t>(id)];
                    got += c * mu(p.x(), p.y());
                }
                INFO("canonical point " << r);
                CHECK(got ==
                      Catch::Approx(
                          mu(panel.cp[r].x(), panel.cp[r].y()))
                          .epsilon(1e-8));
            }
        }
    }
}

TEST_CASE("doublet network spline: all synthetic panels",
          "[splines][outer]") {
    a502::Network net;
    net.netname = "grid";
    net.nm = 8;
    net.nn = 8;
    for (int n = 0; n < 8; ++n) {
        for (int m = 0; m < 8; ++m) {
            net.coordinates.push_back(static_cast<double>(m));
            net.coordinates.push_back(static_cast<double>(n));
            net.coordinates.push_back(0.0);
        }
    }

    const auto flow = a502::Compressibility::fromAngles(0.6, 0.0, 0.0);
    const auto parameters = a502::buildNetworkDoubletParameters(net);
    const auto spline =
        a502::buildNetworkDoubletSpline(net, flow, parameters);
    REQUIRE(spline.valid);
    CHECK(spline.nPanelM == 7);
    CHECK(spline.nPanelN == 7);
    CHECK(spline.nParameters == 81);
    REQUIRE(spline.rows.size() == 9 * 7 * 7);

    // Every row is a convex-like constant-reproducing stencil in exact
    // arithmetic; only the fit coefficients carry roundoff.
    for (const auto& row : spline.rows) {
        REQUIRE_FALSE(row.index.empty());
        double sum = 0.0;
        for (double c : row.coefficient)
            sum += c;
        CHECK(sum == Catch::Approx(1.0).epsilon(1e-9));
    }
}

TEST_CASE("doublet panel spline: sample edge panel", "[splines][outer]") {
    a502::AeroModel model;
    model.parse(samplePath("btac.inp").string());
    REQUIRE(model.solve());
    REQUIRE(model.networks.size() > 2);
    const auto& net = model.networks[2]; // akl, 11x11
    REQUIRE(net.nm == 11);
    REQUIRE(net.nn == 11);

    const auto parameters = a502::buildNetworkDoubletParameters(net);
    REQUIRE(parameters.count == 12 * 12);

    const auto spline =
        a502::buildDoubletPanelSpline(net, model.flow, 0, 4, parameters);
    REQUIRE(spline.valid);
    for (int r = 0; r < 9; ++r) {
        double rowSum = 0.0;
        for (double c : spline.rows[static_cast<std::size_t>(r)].coefficient)
            rowSum += c;
        CHECK(rowSum == Catch::Approx(1.0).epsilon(1e-8));
    }
}

TEST_CASE("doublet panel spline: all sample network panels",
          "[splines][outer]") {
    a502::AeroModel model;
    model.parse(samplePath("btac.inp").string());
    REQUIRE(model.solve());
    REQUIRE(model.networks.size() > 2);
    const auto& net = model.networks[2]; // akl, 11x11
    REQUIRE(net.nm == 11);
    REQUIRE(net.nn == 11);

    const auto parameters = a502::buildNetworkDoubletParameters(net);
    for (int panelJ = 0; panelJ < 10; ++panelJ) {
        for (int panelI = 0; panelI < 10; ++panelI) {
            INFO("panel " << panelI << "," << panelJ);
            const auto spline = a502::buildDoubletPanelSpline(
                net, model.flow, panelI, panelJ, parameters);
            REQUIRE(spline.valid);
            for (int r = 0; r < 9; ++r) {
                double rowSum = 0.0;
                for (double c :
                     spline.rows[static_cast<std::size_t>(r)].coefficient)
                    rowSum += c;
                CHECK(rowSum == Catch::Approx(1.0).epsilon(1e-7));
            }
        }
    }
}

TEST_CASE("doublet network spline: sample decks", "[splines][outer]") {
    const char* decks[] = {"vss.inp", "swb.inp", "ellip.inp"};
    for (const char* deck : decks) {
        INFO("deck " << deck);
        a502::AeroModel model;
        model.parse(samplePath(deck).string());
        REQUIRE(model.solve());

        for (std::size_t k = 0; k < model.networks.size(); ++k) {
            INFO("network " << k);
            const auto& net = model.networks[k];
            const auto parameters =
                a502::buildNetworkDoubletParameters(net);
            REQUIRE(parameters.count > 0);
            const auto spline = a502::buildNetworkDoubletSpline(
                net, model.flow, parameters);
            REQUIRE(spline.valid);
            CHECK(spline.rows.size() ==
                  static_cast<std::size_t>(9) *
                      static_cast<std::size_t>(spline.nPanelM) *
                      static_cast<std::size_t>(spline.nPanelN));
            for (const auto& row : spline.rows) {
                REQUIRE_FALSE(row.index.empty());
                double sum = 0.0;
                for (double c : row.coefficient)
                    sum += c;
                CHECK(sum == Catch::Approx(1.0).epsilon(1e-6));
            }
        }
    }
}

TEST_CASE("edgls: edge quadratic reproduction", "[splines][edge]") {
    const std::array<double, 4> d = {1.0, 1.0, 1.0, 1.0};
    const std::array<double, 4> w = {1.0, 1.0, 1.0, 1.0};
    const auto ls = a502::edgls(d, w);
    REQUIRE(ls.valid);

    const double x1 = -1.5;
    const double x2 = -0.5;
    const double x3 = 0.5;
    const double x4 = 1.5;

    auto quad = [](double x) { return 1.0 + 2.0 * x + 3.0 * x * x; };
    CHECK(ls.a[0] * quad(x1) + ls.a[1] * quad(x2) +
          ls.a[2] * quad(x3) + ls.a[3] * quad(x4) ==
          Catch::Approx(quad(0.0)).epsilon(1e-12));

    auto lin = [](double x) { return 4.0 - 1.5 * x; };
    CHECK(ls.a[0] * lin(x1) + ls.a[1] * lin(x2) +
          ls.a[2] * lin(x3) + ls.a[3] * lin(x4) ==
          Catch::Approx(lin(0.0)).epsilon(1e-12));

    // Non-uniform spacing and weights: a quadratic is still reproduced exactly.
    const std::array<double, 4> du = {0.5, 1.2, 0.8, 0.3};
    const std::array<double, 4> wu = {0.7, 1.0, 1.0, 2.3};
    const auto lsu = a502::edgls(du, wu);
    REQUIRE(lsu.valid);
    const double xu1 = -1.2 - 0.25;
    const double xu2 = -0.6;
    const double xu3 = 0.4;
    const double xu4 = 0.8 + 0.15;
    CHECK(lsu.a[0] * quad(xu1) + lsu.a[1] * quad(xu2) +
          lsu.a[2] * quad(xu3) + lsu.a[3] * quad(xu4) ==
          Catch::Approx(quad(0.0)).epsilon(1e-12));

    // Coincident inner points.
    const std::array<double, 4> dz = {1.0, 0.0, 0.0, 1.0};
    const auto lsz = a502::edgls(dz, w);
    CHECK(lsz.a[0] == 0.0);
    CHECK(lsz.a[1] == 0.5);
    CHECK(lsz.a[2] == 0.5);
    CHECK(lsz.a[3] == 0.0);
}
