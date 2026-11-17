// Tests for the C++ A502 port.
//
// The first block covers the $points network reader: it must follow the same
// card sequence as meshp/inbc in panair.f90 (kn card, kt card, "nm nn
// netname" cards, then exactly nm*nn coordinate points).

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "AeroModel.h"
#include "Compressibility.h"
#include "PanelGeometry.h"

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
