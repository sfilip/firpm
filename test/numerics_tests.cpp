// Unit tests for the numerical kernels: pmmath wrappers and the barycentric
// routines, checked against MPFR where it is available.

#include <vector>
#include <limits>
#include <type_traits>
#include "firpm.h"
#include "gtest/gtest.h"

#ifdef HAVE_MPFR
    using types = testing::Types<double, long double, mpfr::mpreal>;
    template<typename T> constexpr bool is_mpreal = std::is_same_v<T, mpfr::mpreal>;
#else
    using types = testing::Types<double, long double>;
    template<typename T> constexpr bool is_mpreal = false;
#endif

template<typename _T>
struct pmmath_test : public testing::Test { using T = _T; };
TYPED_TEST_SUITE(pmmath_test, types);

TYPED_TEST(pmmath_test, frexp_ldexp_roundtrip)
{
    using T = typename TestFixture::T;
#ifdef HAVE_MPFR
    mpfr::mpreal::set_default_prec(165ul);
#endif
    std::vector<double> samples = {1.0, -1.0, 0.75, 3.0, 1e-300, -2.5e200,
                                   123456.789, 6.02214076e23, 0.0};
    for (double d : samples) {
        T x = T(d);
        int e = 0;
        T m = pm::pmmath::frexp(x, &e);
        if (d == 0.0) {
            ASSERT_EQ(m, T(0));
            continue;
        }
        ASSERT_GE(pm::pmmath::fabs(m), T(0.5));
        ASSERT_LT(pm::pmmath::fabs(m), T(1));
        ASSERT_EQ(pm::pmmath::ldexp(m, e), x);
        ASSERT_EQ(pm::pmmath::ldexp(x, 3), x * 8);
        ASSERT_EQ(pm::pmmath::ldexp(x, -3), x / 8);
    }
}

#ifdef HAVE_MPFR
// The uniform starting reference of the partition1 case from the scaling
// tests: six nodes clustered within 2e-3 of x = 1 and 197 spread over the
// stopband, giving barycentric weights spanning 36 orders of magnitude.
static std::vector<double> clustered_reference()
{
    std::vector<pm::band_t<double>> fbands(2);
    double pi = pm::pmmath::const_pi<double>();
    fbands[0].start = 0;         fbands[0].stop = 0.021 * pi;
    fbands[1].start = 0.112333 * pi; fbands[1].stop = pi;
    for (auto& b : fbands) {
        b.space = pm::space_t::FREQ;
        b.amplitude = [](pm::space_t, double) { return 1.0; };
        b.weight = [](pm::space_t, double) { return 1.0; };
    }
    std::vector<double> omega, x;
    pm::uniform(omega, fbands, 203u);
    pm::cos(x, omega);
    std::ranges::sort(x);
    return x;
}

TEST(barycentric_test, weights_match_mpfr_on_clustered_nodes)
{
    mpfr::mpreal::set_default_prec(256ul);
    std::vector<double> x = clustered_reference();
    std::vector<mpfr::mpreal> xm(x.begin(), x.end());   // exact conversion
    std::vector<double> w(x.size());
    std::vector<mpfr::mpreal> wm(x.size());
    pm::baryweights(w, x);
    pm::baryweights(wm, xm);
    double wmin = 1e300, wmax = 0;
    for (std::size_t i = 0; i < x.size(); ++i) {
        double ref = wm[i].toDouble();
        ASSERT_TRUE(std::isfinite(w[i]));
        ASSERT_LT(std::fabs(w[i] - ref), 1e-13 * std::fabs(ref)) << "i=" << i;
        wmin = std::min(wmin, std::fabs(ref));
        wmax = std::max(wmax, std::fabs(ref));
    }
    ASSERT_GT(wmax / wmin, 1e30);   // the case is only interesting if ill-scaled
}
#endif

#ifdef HAVE_MPFR
// partition1's bands, with the piecewise-linear passband amplitude and the
// band partition points, in the space the exchange works in.
template<typename T>
static std::vector<pm::band_t<T>> partition1_bands()
{
    T pi = pm::pmmath::const_pi<T>();
    std::vector<T> f = {0.0, 0.0126667, 0.018, 0.021};
    std::vector<T> a = {1.0, 1.00055, 1.00111, 1.00151};
    std::vector<pm::band_t<T>> fbands(2);
    fbands[0].start = 0;
    fbands[0].stop = f[3] * pi;
    for (auto fi : f)
        fbands[0].part.push_back(fi * pi);
    fbands[0].amplitude = [f, a, pi](pm::space_t s, T x) -> T {
        if (s == pm::space_t::CHEBY)
            x = pm::pmmath::acos(x);
        x /= pi;
        std::size_t j = (x <= f[1]) ? 0 : (x <= f[2]) ? 1 : 2;
        return a[j] + (a[j + 1] - a[j]) * (x - f[j]) / (f[j + 1] - f[j]);
    };
    fbands[1].start = T(0.112333) * pi;
    fbands[1].stop = pi;
    fbands[1].amplitude = [](pm::space_t, T) -> T { return 0; };
    for (auto& b : fbands) {
        b.space = pm::space_t::FREQ;
        b.weight = [](pm::space_t, T) -> T { return 1; };
    }
    std::vector<pm::band_t<T>> cbands;
    pm::bandconv(cbands, fbands, pm::convdir_t::FROMFREQ);
    return cbands;
}

// Error evaluation at the clustered reference.  The leveled error there is
// around 1e-40, far below double resolution relative to the response, and
// the barycentric weights are badly scaled.  The double result must agree
// with 256-bit MPFR to within the rounding bound of the form comperror()
// uses at each probe: the quotient of sums where the denominator sum is
// well above its rounding floor, the first barycentric form
//
//   P(x) - D(x) = sum_i l_i(x) (C_i - D(x)),  l_i the Lagrange basis,
//
// otherwise, whose error is a small multiple of eps times
// sum_i |l_i(x)| |C_i - D(x)|.  That bound is tight (below 1e-9 relative)
// across most of the stopband, which is where the exchange needs the error;
// with the quotient of sums throughout, the stopband values were rounding
// noise there.  Near the six clustered passband nodes the Lebesgue
// function reaches 1e36, the numerator is below its own rounding floor, and
// no double evaluation can be accurate: there the value only has to be of
// the order of the data.
TEST(barycentric_test, error_matches_mpfr_on_clustered_nodes)
{
    mpfr::mpreal::set_default_prec(256ul);
    using M = mpfr::mpreal;
    std::vector<double> x = clustered_reference();
    std::vector<M> xm(x.begin(), x.end());
    auto cb = partition1_bands<double>();
    auto cbm = partition1_bands<M>();
    // Band edges computed at 256 bits differ from the double ones by an
    // ulp, which can leave a double node just outside its band: give the
    // MPFR bands the double edges (exactly representable) instead.
    for (std::size_t i = 0; i < cb.size(); ++i) {
        cbm[i].start = cb[i].start;
        cbm[i].stop = cb[i].stop;
        cbm[i].part.assign(cb[i].part.begin(), cb[i].part.end());
    }

    std::vector<double> w(x.size()), D(x.size()), S(x.size());
    std::vector<M> wm(x.size()), Dm(x.size()), Sm(x.size());
    double delta; M deltam;
    pm::baryweights(w, x);   pm::compdelta(delta, w, x, cb);   pm::compref(D, S, x, cb);
    pm::baryweights(wm, xm); pm::compdelta(deltam, wm, xm, cbm); pm::compref(Dm, Sm, xm, cbm);
    ASSERT_LT(std::fabs(delta - deltam.toDouble()), 1e-9 * std::fabs(deltam.toDouble()));
    ASSERT_LT(std::fabs(delta), 1e-30);

    const double eps = std::numeric_limits<double>::epsilon();
    std::size_t tight = 0, probes = 0, informative = 0;
    for (auto const& b : cb) {
        for (int k = 1; k < 40; ++k) {
            double xv = b.start + (b.stop - b.start) * k / 40.0;
            double err; M errm;
            pm::comperror(err, xv, delta, x, D, S, w, cb);
            pm::comperror(errm, M(xv), deltam, xm, Dm, Sm, wm, cbm);
            ++probes;

            M Dx, Wx;
            pm::idealvals(Dx, Wx, M(xv), cbm);
            M num = 0, sabs = 0, sden = 0, dabs = 0, cmax = 0;
            for (std::size_t i = 0; i < x.size(); ++i) {
                M bi = wm[i] / (M(xv) - xm[i]);
                M ci = (Dm[i] - Dx) + deltam * Sm[i];
                num += bi * ci;
                sabs += pm::pmmath::fabs(bi * ci);
                sden += bi;
                dabs += pm::pmmath::fabs(bi);
                cmax = pm::pmmath::fabs(ci) > cmax ? pm::pmmath::fabs(ci) : cmax;
            }
            if (!(pm::pmmath::fabs(num) > sabs * (eps * 64))) {
                // below the numerator's rounding floor: only boundedness
                if (std::isfinite(err))
                    ASSERT_LT(std::fabs(err), 1e3 * (Wx * cmax).toDouble()) << "x=" << xv;
                continue;
            }
            ++informative;
            ASSERT_TRUE(std::isfinite(err)) << "x=" << xv;
            double bound;
            if (pm::pmmath::fabs(sden) > dabs * (eps * 1e6))
                // quotient of two sums: |d num| / |sden| + |num| |d sden| / sden^2
                bound = 64 * eps * (Wx * (sabs + pm::pmmath::fabs(num) * dabs
                            / pm::pmmath::fabs(sden)) / pm::pmmath::fabs(sden)).toDouble();
            else
                // first form: eps * W * sum |l_i (C_i - D)|
                bound = 64 * eps * (Wx * sabs / pm::pmmath::fabs(sden)).toDouble();
            ASSERT_LT(std::fabs(err - errm.toDouble()), bound)
                << "x=" << xv << " double=" << err << " mpfr=" << errm.toDouble();
            if (bound < 1e-9 * std::fabs(errm.toDouble()))
                ++tight;
        }
    }
    ASSERT_GE(informative, 35u) << "of " << probes << " probes";
    ASSERT_GE(tight, 25u) << "of " << informative << " informative probes";

    for (std::size_t i = 0; i < x.size(); ++i) {
        double err;
        pm::comperror(err, x[i], delta, x, D, S, w, cb);
        ASSERT_EQ(err, (i % 2 == 0) ? delta : -delta);
    }
}
#endif
