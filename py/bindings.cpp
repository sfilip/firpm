#include "firpm.h"
#include "firpm/pmmath.h"

#include <pybind11/numpy.h>
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

/* Turns a design's status into a ValueError (as scipy.signal.remez raises
 * for invalid specifications and failed designs alike) or RuntimeWarning,
 * and returns its taps.  Every entry point goes through here, so that no
 * failure comes back as a (possibly empty) tap vector. */
template <typename T>
static py::array_t<double> taps(pm::pmoutput_t<T> const& output)
{
    const char* failure = nullptr;
    switch(output.status) {
    case pm::status_t::STATUS_SUCCESS:
        break;
    case pm::status_t::STATUS_CONVERGENCE_WARNING:
        if(PyErr_WarnEx(PyExc_RuntimeWarning,
                "the exchange algorithm did not converge", 1) < 0)
            throw py::error_already_set();
        break;
    case pm::status_t::STATUS_FREQUENCY_INVALID_INTERVAL:
        failure = "invalid band edges: they must come in pairs, be "
            "nondecreasing and lie between 0 and Nyquist, and not every band "
            "can be a single point";
        break;
    case pm::status_t::STATUS_AMPLITUDE_VECTOR_MISMATCH:
        failure = "there must be one amplitude per band edge";
        break;
    case pm::status_t::STATUS_AMPLITUDE_DISCONTINUITY:
        failure = "adjacent bands that share an edge must have the same "
            "amplitude there";
        break;
    case pm::status_t::STATUS_WEIGHT_NEGATIVE:
        failure = "band weights must be positive";
        break;
    case pm::status_t::STATUS_WEIGHT_VECTOR_MISMATCH:
        failure = "there must be one weight per band";
        break;
    case pm::status_t::STATUS_WEIGHT_DISCONTINUITY:
        failure = "incompatible weights within a band";
        break;
    case pm::status_t::STATUS_SCALING_INVALID:
        failure = "reference scaling failed to produce a starting reference";
        break;
    case pm::status_t::STATUS_AFP_INVALID:
        failure = "AFP initialization failed to produce a starting reference";
        break;
    case pm::status_t::STATUS_COEFFICIENT_SET_INVALID:
        failure = "the design produced an invalid set of coefficients";
        break;
    case pm::status_t::STATUS_EXCHANGE_FAILURE:
        failure = "the exchange algorithm failed to produce a valid reference";
        break;
    case pm::status_t::STATUS_UNKNOWN_FAILURE:
        failure = "the design failed";
        break;
    }
    if(failure)
        throw py::value_error(failure);

    py::array_t<double> h(output.h.size());
    auto hd = h.mutable_unchecked<1>();
    for(std::size_t i{0u}; i < output.h.size(); ++i)
        hd(i) = to_double(output.h[i]);
    return h;
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

/* The band callbacks see omega in [0, pi] (FREQ) or cos(omega) (CHEBY);
 * Python sees frequency in the units of its band edges, in [0, nyquist]. */
template <typename T>
static T frequency(pm::space_t s, T x, T nyquist)
{
    if(s == pm::space_t::CHEBY)
        x = pm::pmmath::acos(x);
    return x / pm::pmmath::const_pi<T>() * nyquist;
}

template <typename T>
static T call_python(py::function const& fn, callback_state& state, T f)
{
    if(state.failed.load(std::memory_order_acquire))
        return T(std::numeric_limits<double>::quiet_NaN());

    py::gil_scoped_acquire gil;
    try {
        return T(fn(to_double(f)).template cast<double>());
    } catch(py::error_already_set& e) {
        if(!state.failed.exchange(true, std::memory_order_acq_rel))
            state.pyerr.reset(new py::error_already_set(std::move(e)));
    } catch(std::exception& e) {
        if(!state.failed.exchange(true, std::memory_order_acq_rel))
            state.cpperr = e.what();
    }
    return T(std::numeric_limits<double>::quiet_NaN());
}

/* A band's amplitude and weight are each a constant or a Python callable.
 * The firpm package sends purely numeric specifications to firpm_impl, so
 * these only arrive when at least one band has a callable. */
using band_spec_t = std::tuple<double, double, py::object, py::object>;

template <typename T>
static std::function<T(pm::space_t, T)> band_function(py::object const& v,
        T nyquist, std::shared_ptr<callback_state> const& state)
{
    if(!PyCallable_Check(v.ptr())) {
        T value = v.cast<double>();
        return [value](pm::space_t, T) -> T { return value; };
    }
    /* pm::firpm copies band_t (and so these closures) freely while the
     * GIL is released.  Holding the callable through shared_ptr keeps
     * those copies off the Python reference counts; the original in the
     * band specification outlives the call and is released with the GIL
     * held. */
    auto fp = std::make_shared<py::function>(
            py::reinterpret_borrow<py::function>(v));
    return [fp, nyquist, state](pm::space_t s, T x) -> T {
        return call_python<T>(*fp, *state, frequency(s, x, nyquist));
    };
}

/* Parks-McClellan with per-band amplitude and weight functions, for type I
 * (even n) and type II (odd n).  Band edges run from 0 to nyquist, in the
 * caller's units.  Delegates to pm::firpm(n, fbands), which validates the
 * bands and handles the cos(omega/2) basis change and tap recovery
 * internally. */
template <typename T>
static py::array_t<double> firpm_bands_impl(std::size_t n,
        std::vector<band_spec_t> const& bands, double nyquist, double eps,
        std::size_t nmax, pm::init_t strategy, std::size_t depth,
        pm::init_t rstrategy, unsigned long prec)
{
    auto state = std::make_shared<callback_state>();
    T ny = nyquist;

    std::vector<pm::band_t<T>> fbands(bands.size());
    for(std::size_t i{0u}; i < bands.size(); ++i) {
        py::object const& amplitude = std::get<2>(bands[i]);
        py::object const& weight = std::get<3>(bands[i]);
        /* As the flat specification's parseSpecification() requires. */
        if(!PyCallable_Check(weight.ptr()) && !(weight.cast<double>() > 0.0))
            throw py::value_error("band weights must be positive");

        fbands[i].start = pm::pmmath::const_pi<T>() * T(std::get<0>(bands[i])) / ny;
        fbands[i].stop  = pm::pmmath::const_pi<T>() * T(std::get<1>(bands[i])) / ny;
        fbands[i].part  = {fbands[i].start, fbands[i].stop};
        fbands[i].space = pm::space_t::FREQ;
        fbands[i].amplitude = band_function<T>(amplitude, ny, state);
        fbands[i].weight = band_function<T>(weight, ny, state);
    }

    pm::pmoutput_t<T> output = [&]() {
        py::gil_scoped_release nogil;
        return pm::firpm<T>(n, fbands, eps, nmax, strategy, depth, rstrategy,
                prec);
    }();

    if(state->pyerr)
        throw py::error_already_set(*state->pyerr);
    if(!state->cpperr.empty())
        throw std::runtime_error("callback failed: " + state->cpperr);

    return taps(output);
}

/* Parks-McClellan with MATLAB-style band edges, amplitudes and weights.
 * type is "bandpass" (types I and II), or "differentiator" or "hilbert"
 * (types III and IV). */
template <typename T>
static py::array_t<double> firpm_impl(std::size_t n,
        std::vector<double> const& f, std::vector<double> const& a,
        std::vector<double> const& w, std::string const& type, double eps,
        std::size_t nmax, pm::init_t strategy, std::size_t depth,
        pm::init_t rstrategy, unsigned long prec)
{
    bool bandpass = type == "bandpass";
    pm::filter_t ftype = pm::filter_t::FIR_DIFFERENTIATOR;
    if(type == "hilbert")
        ftype = pm::filter_t::FIR_HILBERT;
    else if(!bandpass && type != "differentiator")
        throw py::value_error("unknown filter type: " + type);

    std::vector<T> ft(f.begin(), f.end()), at(a.begin(), a.end()),
        wt(w.begin(), w.end());

    pm::pmoutput_t<T> output = [&]() {
        py::gil_scoped_release nogil;
        if(bandpass)
            return firpm<T>(n, ft, at, wt, eps, nmax, strategy, depth,
                    rstrategy, prec);
        return firpm<T>(n, ft, at, wt, ftype, eps, nmax, strategy, depth,
                rstrategy, prec);
    }();
    return taps(output);
}

/* The private half of the firpm package: firpm.firpm() checks and
 * dispatches its arguments to these. */
PYBIND11_MODULE(_firpm, m){

    py::enum_<pm::init_t>(m, "Strategy")
        .value("UNIFORM", pm::init_t::UNIFORM)
        .value("SCALING", pm::init_t::SCALING)
        .value("AFP", pm::init_t::AFP);

    m.attr("__version__") = FIRPM_VERSION;

    m.def("firpm", &firpm_impl<double>,
        "n"_a, "f"_a, "a"_a, "w"_a, "type"_a, "eps"_a, "nmax"_a, "strategy"_a,
        "depth"_a, "rstrategy"_a, "prec"_a);
    m.def("firpm_bands", &firpm_bands_impl<double>,
        "n"_a, "bands"_a, "nyquist"_a, "eps"_a, "nmax"_a, "strategy"_a,
        "depth"_a, "rstrategy"_a, "prec"_a);

#if HAVE_MPFR
    /* Values of mpreal take the default precision when constructed, so it
     * must be set before the specification is converted. */
    m.def("firpm_mp", [](std::size_t n, std::vector<double> const& f,
              std::vector<double> const& a, std::vector<double> const& w,
              std::string const& type, double eps, std::size_t nmax,
              pm::init_t strategy, std::size_t depth, pm::init_t rstrategy,
              unsigned long prec) {
            mpfr::mpreal::set_default_prec(prec);
            return firpm_impl<mpfr::mpreal>(n, f, a, w, type, eps, nmax,
                    strategy, depth, rstrategy, prec);
        },
        "n"_a, "f"_a, "a"_a, "w"_a, "type"_a, "eps"_a, "nmax"_a, "strategy"_a,
        "depth"_a, "rstrategy"_a, "prec"_a);
    m.def("firpm_bands_mp", [](std::size_t n,
              std::vector<band_spec_t> const& bands, double nyquist,
              double eps, std::size_t nmax, pm::init_t strategy,
              std::size_t depth, pm::init_t rstrategy, unsigned long prec) {
            mpfr::mpreal::set_default_prec(prec);
            return firpm_bands_impl<mpfr::mpreal>(n, bands, nyquist, eps,
                    nmax, strategy, depth, rstrategy, prec);
        },
        "n"_a, "bands"_a, "nyquist"_a, "eps"_a, "nmax"_a, "strategy"_a,
        "depth"_a, "rstrategy"_a, "prec"_a);
#endif
}
