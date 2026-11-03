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
