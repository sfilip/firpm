#include "firpm.h"
#include "firpm/pmmath.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <atomic>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>

namespace py = pybind11;
using pm::firpm;
using namespace pybind11::literals;

static double to_double(double x) { return x; }
static double to_double(long double x) { return static_cast<double>(x); }
#if HAVE_MPFR
static double to_double(const mpfr::mpreal& x) { return x.toDouble(); }
#endif

static const char* status_name(pm::status_t s)
{
    switch(s) {
    case pm::status_t::STATUS_SUCCESS:                    return "STATUS_SUCCESS";
    case pm::status_t::STATUS_FREQUENCY_INVALID_INTERVAL: return "STATUS_FREQUENCY_INVALID_INTERVAL";
    case pm::status_t::STATUS_AMPLITUDE_VECTOR_MISMATCH:  return "STATUS_AMPLITUDE_VECTOR_MISMATCH";
    case pm::status_t::STATUS_AMPLITUDE_DISCONTINUITY:    return "STATUS_AMPLITUDE_DISCONTINUITY";
    case pm::status_t::STATUS_WEIGHT_NEGATIVE:            return "STATUS_WEIGHT_NEGATIVE";
    case pm::status_t::STATUS_WEIGHT_VECTOR_MISMATCH:     return "STATUS_WEIGHT_VECTOR_MISMATCH";
    case pm::status_t::STATUS_WEIGHT_DISCONTINUITY:       return "STATUS_WEIGHT_DISCONTINUITY";
    case pm::status_t::STATUS_SCALING_INVALID:            return "STATUS_SCALING_INVALID";
    case pm::status_t::STATUS_AFP_INVALID:                return "STATUS_AFP_INVALID";
    case pm::status_t::STATUS_COEFFICIENT_SET_INVALID:    return "STATUS_COEFFICIENT_SET_INVALID";
    case pm::status_t::STATUS_EXCHANGE_FAILURE:           return "STATUS_EXCHANGE_FAILURE";
    case pm::status_t::STATUS_CONVERGENCE_WARNING:        return "STATUS_CONVERGENCE_WARNING";
    case pm::status_t::STATUS_UNKNOWN_FAILURE:            return "STATUS_UNKNOWN_FAILURE";
    }
    return "STATUS_UNKNOWN";
}

/* Band callbacks are evaluated from inside the OpenMP parallel regions of
 * exchange().  An exception escaping an OpenMP structured block terminates
 * the process, so a Python (or C++) exception raised by a callback cannot be
 * propagated directly.  Instead the first one is captured here, every later
 * callback short-circuits with NaN (which makes the exchange loop exit at the
 * next iteration), and the exception is re-raised on the calling thread once
 * pm::firpm has returned. */
struct callback_state {
    std::atomic<bool> failed{false};
    std::unique_ptr<py::error_already_set> pyerr;   /* touched only with the GIL held */
    std::string cpperr;
};

template <typename T>
static T call_python(py::function const& fn, callback_state& state,
        pm::space_t s, T x)
{
    if(state.failed.load(std::memory_order_acquire))
        return T(std::numeric_limits<double>::quiet_NaN());
    if(s == pm::space_t::CHEBY)
        x = pm::pmmath::acos(x);

    py::gil_scoped_acquire gil;
    try {
        return T(fn(to_double(x)).template cast<double>());
    } catch(py::error_already_set& e) {
        if(!state.failed.exchange(true, std::memory_order_acq_rel))
            state.pyerr.reset(new py::error_already_set(std::move(e)));
    } catch(std::exception& e) {
        if(!state.failed.exchange(true, std::memory_order_acq_rel))
            state.cpperr = e.what();
    }
    return T(std::numeric_limits<double>::quiet_NaN());
}

/* Custom-band Parks-McClellan: per-band amplitude and weight are arbitrary
 * Python callables evaluated at omega in [0, pi] (FREQ space).  Supports
 * type I (even n) and type II (odd n).  Delegates to pm::firpm(n, fbands)
 * which validates the bands and handles the cos(omega/2) basis change and
 * tap recovery internally. */
template <typename T>
static std::vector<double> firpm_bands_impl(std::size_t n,
        std::vector<std::tuple<double, double, py::function, py::function>> bands,
        double eps, std::size_t nmax, unsigned long prec)
{
    auto state = std::make_shared<callback_state>();

    std::vector<pm::band_t<T>> fbands(bands.size());
    for(std::size_t i{0u}; i < bands.size(); ++i) {
        double f0 = std::get<0>(bands[i]), f1 = std::get<1>(bands[i]);
        /* pm::firpm copies band_t (and so these closures) freely while the
         * GIL is released.  Holding the callables through shared_ptr keeps
         * those copies off the Python reference counts; the originals in
         * `bands` outlive the call and are released with the GIL held. */
        auto amp = std::make_shared<py::function>(std::get<2>(bands[i]));
        auto wt  = std::make_shared<py::function>(std::get<3>(bands[i]));
        fbands[i].start = pm::pmmath::const_pi<T>() * T(f0);
        fbands[i].stop  = pm::pmmath::const_pi<T>() * T(f1);
        fbands[i].part  = {fbands[i].start, fbands[i].stop};
        fbands[i].space = pm::space_t::FREQ;
        fbands[i].amplitude = [amp, state](pm::space_t s, T x) -> T {
            return call_python<T>(*amp, *state, s, x);
        };
        fbands[i].weight = [wt, state](pm::space_t s, T x) -> T {
            return call_python<T>(*wt, *state, s, x);
        };
    }

    pm::pmoutput_t<T> output = [&]() {
        py::gil_scoped_release nogil;
        return pm::firpm<T>(n, fbands, eps, nmax, pm::init_t::UNIFORM, 0u,
                pm::init_t::UNIFORM, prec);
    }();

    if(state->pyerr)
        throw py::error_already_set(*state->pyerr);
    if(!state->cpperr.empty())
        throw std::runtime_error("firpm_bands: callback failed: " + state->cpperr);

    switch(output.status) {
    case pm::status_t::STATUS_SUCCESS:
        break;
    case pm::status_t::STATUS_CONVERGENCE_WARNING:
        if(PyErr_WarnEx(PyExc_RuntimeWarning,
                "firpm_bands: the exchange algorithm did not converge", 1) < 0)
            throw py::error_already_set();
        break;
    default:
        throw std::runtime_error(std::string("firpm_bands failed: ") +
                status_name(output.status));
    }

    std::vector<double> hd(output.h.size());
    for(std::size_t i{0u}; i < output.h.size(); ++i)
        hd[i] = to_double(output.h[i]);
    return hd;
}

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

    m.def("firpm_bands", &firpm_bands_impl<double>,
        "n"_a, "bands"_a, "eps"_a=0.01, "nmax"_a=4, "prec"_a=165ul,
        "Parks-McClellan with per-band callable amplitude/weight (type I/II, "
        "double precision).  bands is a list of (f0, f1, amplitude, weight) "
        "tuples with f in [0, 1] and callables evaluated at omega in [0, pi].  "
        "Returns tap vector h; raises RuntimeError on failure and re-raises "
        "exceptions thrown by the callbacks.");

#if HAVE_MPFR
    m.def("firpm_bands_mp", [](std::size_t n,
              std::vector<std::tuple<double, double, py::function, py::function>> bands,
              double eps, std::size_t nmax, unsigned long prec) {
            mpfr::mpreal::set_default_prec(prec);
            return firpm_bands_impl<mpfr::mpreal>(n, bands, eps, nmax, prec);
        },
        "n"_a, "bands"_a, "eps"_a=0.01, "nmax"_a=4, "prec"_a=165ul,
        "As firpm_bands (type I/II), computed in MPFR arbitrary precision (prec bits); "
        "taps extracted before rounding to double.  Returns tap vector h.");
#endif
}
