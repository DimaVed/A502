// main.cpp
// Small demonstration: parse a PANAIR (A502) input deck and pretty-print
// the highlights of the model, including per-network mesh statistics.
//
// usage: panair_parser_bench <input-file>
#include <cstdio>
#include <string>

#include "panair_input.h"
#include "panair_parser.h"

using panair::PanairInput;
using panair::PanairParser;

static void dump(const PanairInput& in) {
    std::printf("title1 : %s\n", in.title1.c_str());
    std::printf("title2 : %s\n", in.title2.c_str());
    std::printf("ndtchk : %d      misym : %d   mjsym : %d\n",
                in.ndtchk, in.misym, in.mjsym);
    std::printf("amach  : %g      nacase: %d\n", in.amach, in.nacase);
    std::printf("alpc   : %g   alpha: ", in.alpc);
    for (double a : in.alpha) std::printf(" %g", a);
    std::printf("\nbetc   : %g   beta : ", in.betc);
    for (double b : in.beta) std::printf(" %g", b);
    std::printf("\nrefs   : xref=%g yref=%g zref=%g nref=%d nprcof=%d\n",
                in.refs.xref, in.refs.yref, in.refs.zref,
                in.refs.nref, in.refs.nprcof);
    std::printf("refs   : sref=%g bref=%g cref=%g dref=%g\n",
                in.refs.sref, in.refs.bref, in.refs.cref, in.refs.dref);

    std::printf("\n-- networks (%zu) --\n", in.networks.size());
    for (const auto& n : in.networks) {
        std::printf("  net %-10s kn=%d kt=%d nm=%d nn=%d cpnorm=%g amnsw=%g dnsmsh=%g\n",
                    n.netname.c_str(), n.kn, n.kt, n.nm, n.nn,
                    n.cpnorm, n.amnsw, n.dnsmsh);
        if (!n.points.empty()) {
            // bounding box + a sanity check through Eigen
            double xmin=1e30, xmax=-1e30, ymin=1e30, ymax=-1e30;
            for (const auto& p : n.points) {
                xmin = std::min(xmin, p.x()); xmax = std::max(xmax, p.x());
                ymin = std::min(ymin, p.y()); ymax = std::max(ymax, p.y());
            }
            std::printf("       bbox x[%g,%g] y[%g,%g]  spanext=(%g,%g)\n",
                        xmin, xmax, ymin, ymax, xmax - xmin, ymax - ymin);
        }
    }

    std::printf("\n-- trailing wakes (%zu) --\n", in.wakes.size());
    for (const auto& w : in.wakes)
        std::printf("  %-10s insd=%g xwake=%g twake=%g  -> %-10s\n",
                    w.refNet.c_str(), w.insd, w.xwake, w.twake, w.netname.c_str());

    std::printf("\n-- raw (undeclared) sections (%zu) --\n", in.rawSections.size());
    for (const auto& r : in.rawSections) {
        std::printf("  %-5s (%zu data cards) %s\n",
                    r.keyword.c_str(), r.cards.size(), r.description.c_str());
    }
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::fprintf(stderr, "usage: %s <input-file>\n", argv[0]);
        return 2;
    }
    try {
        PanairParser parser;
        PanairInput in = parser.parse(argv[1]);
        dump(in);
        return 0;
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        return 1;
    }
}