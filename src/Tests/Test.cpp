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
#include "PanelParameters.h"
#include "OuterSplines.h"
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
    const Eigen::Vector3d local = frame * (p - sp.origin);
    return a502::evaluateSource(
        a502::buildSourceSubpanelSplines(g)[static_cast<std::size_t>(s)],
        local.head<2>(), v);
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
        REQUIRE(splines[static_cast<std::size_t>(s)].valid);
        const auto& sp = g.subpanels[static_cast<std::size_t>(s)];
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
