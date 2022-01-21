#include "firpm.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using pm::firpm;
using namespace pybind11::literals;

PYBIND11_MODULE(pyfirpm, m){

    py::enum_<pm::init_t>(m, "Strategy")
        .value("UNIFORM", pm::init_t::UNIFORM)
        .value("SCALING", pm::init_t::SCALING)
        .value("AFP", pm::init_t::AFP)
        .export_values();

    m.def("firpm", [](std::size_t n,
              std::vector<double> f, std::vector<double> a, std::vector<double> w,
              double eps,
              std::size_t nmax,
              pm::init_t strategy,
              std::size_t depth,
              pm::init_t rstrategy,
              unsigned long prec) {
            return firpm<double>(n, f, a, w, eps, nmax, strategy, depth, rstrategy, prec).h;
        },
        "n"_a,
        "f"_a,
        "a"_a,
        "w"_a,
        "eps"_a=0.01,
        "nmax"_a=4,
        "strategy"_a=pm::init_t::UNIFORM,
        "depth"_a=0u,
        "rstrategy"_a=pm::init_t::UNIFORM,
        "prec"_a=165ul,
        "Parks-McClellan routine for implementing type I and II FIR filters.");
}
