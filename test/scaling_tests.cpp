#include <vector>
#include <fstream>
#include <chrono>
#include <type_traits>
#include "firpm.h"
#include "gtest/gtest.h"

using testing::Types;
using pm::firpm;
using pm::firpmRS;
using pm::firpmAFP;
using pm::pmoutput_t;
using pm::filter_t;

#ifdef HAVE_MPFR
    #include <unsupported/Eigen/MPRealSupport>
    using types = testing::Types<double, long double, mpfr::mpreal>;
#else
    using types = testing::Types<double, long double>;
#endif

template<typename _T>
struct firpm_scaling_test : public testing::Test { using T = _T; };
TYPED_TEST_SUITE(firpm_scaling_test, types);

template<typename _T>
struct firpm_lebesgue_test : public testing::Test { using T = _T; };
TYPED_TEST_SUITE(firpm_lebesgue_test, types);

template<typename _T>
struct firpm_cic_test : public testing::Test { using T = _T; };
TYPED_TEST_SUITE(firpm_cic_test, types);

template<typename _T>
struct firpm_issues_test : public testing::Test { using T = _T; };
TYPED_TEST_SUITE(firpm_issues_test, types);

// example of how to use the exchange method directly
TYPED_TEST(firpm_scaling_test, lowpass50a)
{
#ifdef HAVE_MPFR
    mpfr::mpreal::set_default_prec(165ul);
#endif

    using T = typename TestFixture::T;

    std::vector<pm::band_t<T>> freqBands(2);

    freqBands[0].start = 0;
    freqBands[0].stop = pm::pmmath::const_pi<T>() * 0.4;
    freqBands[0].weight = [] (pm::space_t, T) -> T {return 1; };
    freqBands[0].space = pm::space_t::FREQ;
    freqBands[0].amplitude = [](pm::space_t, T) -> T { return 1; };

    freqBands[1].start = pm::pmmath::const_pi<T>() * 0.5;
    freqBands[1].stop = pm::pmmath::const_pi<T>();
    freqBands[1].weight = [] (pm::space_t, T) -> T {return 1; };
    freqBands[1].space = pm::space_t::FREQ;
    freqBands[1].amplitude = [](pm::space_t, T) -> T { return 0; };



    auto start = std::chrono::steady_clock::now();

    std::size_t degree = 50;
    std::vector<T> a;
    std::vector<pm::band_t<T>> chebyBands;
    std::vector<T> omega(degree + 2u);
    std::vector<T> x(degree + 2u);
    pm::uniform(omega, freqBands, degree + 2u);
    pm::cos(x, omega);
    pm::bandconv(chebyBands, freqBands, pm::convdir_t::FROMFREQ);

    auto output = pm::exchange(x, chebyBands);
    ASSERT_LT(output.q, 1e-2);

    for(std::size_t counter = 0; counter < 1; ++counter) {
        std::vector<T> newX;
        pm::refscaling(output.status, newX, chebyBands, freqBands, 2 * degree + 2,
                output.x, chebyBands, freqBands);
        degree = 2 * degree;

        output = pm::exchange(newX, chebyBands);
        ASSERT_LT(output.q, 1e-2);
    }
        auto stop  = std::chrono::steady_clock::now();
        double elapsedTime = std::chrono::duration_cast<
            std::chrono::duration<double>>(stop - start).count();
        std::cout << "Elapsed time = " << elapsedTime << std::endl;

}

// Filters appearing inside Section 1 & Section 4 (minus the degree 53348 one)
TYPED_TEST(firpm_scaling_test, lowpass50)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {1.0, 1.0};

    std::size_t degree = 50;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto start = std::chrono::steady_clock::now();
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    auto stop  = std::chrono::steady_clock::now();
    double elapsedTime = std::chrono::duration_cast<
        std::chrono::duration<double>>(stop - start).count();
    std::cout << "Elapsed time = " << elapsedTime << std::endl;

    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    start = std::chrono::steady_clock::now();
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    stop  = std::chrono::steady_clock::now();
    elapsedTime = std::chrono::duration_cast<
        std::chrono::duration<double>>(stop - start).count();
    std::cout << "Elapsed time = " << elapsedTime << std::endl;

    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    start = std::chrono::steady_clock::now();
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);
    stop  = std::chrono::steady_clock::now();
    elapsedTime = std::chrono::duration_cast<
        std::chrono::duration<double>>(stop - start).count();
    std::cout << "Elapsed time = " << elapsedTime << std::endl;

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TYPED_TEST(firpm_scaling_test, lowpass80)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {1.0, 1.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TYPED_TEST(firpm_scaling_test, lowpass100)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {1.0, 1.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}



TYPED_TEST(firpm_scaling_test, bandstop50)
{
    using T = typename TestFixture::T;
    std::vector<T> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<T> w = {1.0, 1.0, 1.0};

    std::size_t degree = 50;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TYPED_TEST(firpm_scaling_test, bandstop80)
{
    using T = typename TestFixture::T;
    std::vector<T> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<T> w = {1.0, 1.0, 1.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TYPED_TEST(firpm_scaling_test, bandstop100)
{
    using T = typename TestFixture::T;
    std::vector<T> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<T> w = {1.0, 1.0, 1.0};

    std::size_t degree = 100;
    pm::pmoutput_t<T> output1;
    if(!std::is_same<T, mpfr::mpreal>::value) {
        std::cout << "START Parks-McClellan with uniform initialization\n";
        output1 = firpm<T>(degree * 2u, f, a, w);
        std::cout << "Final Delta     = " << output1.delta << std::endl;
        std::cout << "Iteration count = " << output1.iter  << std::endl;
        std::cout << "FINISH Parks-McClellan with uniform initialization\n";
        ASSERT_LT(output1.q, 1e-2);
    }

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output2.delta-output3.delta)/output2.delta), 2e-2);

    if(!std::is_same<T, mpfr::mpreal>::value) {
        std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
        std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;
        ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);
    }

}


TYPED_TEST(firpm_scaling_test, combfir)
{
    using T = typename TestFixture::T;
    std::vector<T> f = {0.0, 0.99, 1.0, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {1.0, 1.0};

    std::size_t degree = 520;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

// might take dozens of minutes to compute, depending on your system
// prefer to comment it out
/*TYPED_TEST(firpm_scaling_test, hugefilter)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 1.0 / 8192, 3.0 / 8192, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {1.0, 1.0};

    if(std::is_same<T, long double>::value) {
        std::size_t degree = 53248u;
        std::cout << "START Parks-McClellan with reference scaling\n";
        auto output2 = firpmRS<T>(degree * 2u + 1u, f, a, w, 0.01, 4u, 2u);
        std::cout << "Final Delta     = " << output2.delta << std::endl;
        std::cout << "Iteration count = " << output2.iter  << std::endl;
        std::cout << "FINISH Parks-McClellan with reference scaling\n";
        ASSERT_LT(output2.q, 1e-2);
    } else {
        std::cout << "This test is not run for the double and mpfr::mpreal versions of the code.\n"
            << "REASON:\n"
            << " * double:       failure to converge\n"
            << " * mpfr::mpreal: takes too long to compute\n";
    }
}*/

TYPED_TEST(firpm_lebesgue_test, lowpass500)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.49, 0.5, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {1.0, 10.0};

    std::size_t degree = 500;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}


TYPED_TEST(firpm_lebesgue_test, lowpass1000)
{
    using T = typename TestFixture::T;
    std::vector<T> f = {0.0, 0.49, 0.5, 1.0};
    std::vector<T> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {1.0, 10.0};

    std::size_t degree = 1000;

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);

    pm::pmoutput_t<T> output3;
    // takes way too long with the mpfr version of the code to generate
    // the AFP initialization (AFP is meant to be used in conjunction with
    // double/long double versions of the code for filter lengths that
    // are not too high (e.g. < 5000))
    if(!std::is_same<T, mpfr::mpreal>::value) {
        std::cout << "START Parks-McClellan with AFP\n";
        output3 = firpmAFP<T>(degree * 2u, f, a, w);
        std::cout << "Final Delta     = " << output3.delta << std::endl;
        std::cout << "Iteration count = " << output3.iter  << std::endl;
        std::cout << "FINISH Parks-McClellan with AFP\n";
        ASSERT_LT(output3.q, 1e-2);
        ASSERT_LE(pm::pmmath::fabs((output2.delta-output3.delta)/output2.delta), 2e-2);
    }
}

TYPED_TEST(firpm_lebesgue_test, bandpass60)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<T> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {10.0, 1.0, 5.0};

    std::size_t degree = 60;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}



TYPED_TEST(firpm_lebesgue_test, bandpass70)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<T> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {10.0, 1.0, 5.0};

    std::size_t degree = 70;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}


TYPED_TEST(firpm_lebesgue_test, bandpass80)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<T> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {10.0, 1.0, 5.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}


TYPED_TEST(firpm_lebesgue_test, bandpass100)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<T> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {10.0, 1.0, 5.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TYPED_TEST(firpm_lebesgue_test, multiband100)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<T> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TYPED_TEST(firpm_lebesgue_test, multiband200)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<T> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 200;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TYPED_TEST(firpm_lebesgue_test, multiband300)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<T> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 300;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TYPED_TEST(firpm_lebesgue_test, multiband600)
{
    using T = typename TestFixture::T;

    std::vector<T> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<T> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<T> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 600;
    // convergence with uniform initialization is not guaranteed

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, f, a, w, 1e-5);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-5);

    pm::pmoutput_t<T> output3;
    // idem to lowpass1000
    if(!std::is_same<T, mpfr::mpreal>::value) {
        std::cout << "START Parks-McClellan with AFP\n";
        output3 = firpmAFP<T>(degree * 2u, f, a, w, 1e-5);
        std::cout << "Final Delta     = " << output3.delta << std::endl;
        std::cout << "Iteration count = " << output3.iter  << std::endl;
        std::cout << "FINISH Parks-McClellan with AFP\n";
        ASSERT_LT(output3.q, 1e-5);
        ASSERT_LE(pm::pmmath::fabs((output2.delta-output3.delta)/output2.delta), 2e-2);
    }
}

TYPED_TEST(firpm_cic_test, cic119) {
    using T = typename TestFixture::T;
    std::vector<T>cf_freq = {
    0.000000000000000e+00, 1.666666666666667e-02, 3.333333333333333e-02,
    5.000000000000000e-02, 6.666666666666667e-02, 8.333333333333333e-02,
    1.000000000000000e-01, 1.166666666666667e-01, 1.333333333333333e-01,
    1.500000000000000e-01, 1.666666666666667e-01, 1.833333333333333e-01,
    2.000000000000000e-01, 1.000000000000000e+00
    };

    std::vector<T>cf_mag = {
    1.000000000000000e+00, 1.000237504182675e+00, 1.000950371083591e+00,
    1.002139664683764e+00, 1.003807161362865e+00, 1.005955354528478e+00,
    1.008587461121836e+00, 1.011707430024841e+00, 1.015319952400556e+00,
    1.019430474006924e+00, 1.024045209531033e+00, 1.029171158999276e+00,
    1.034816126326940e+00, 0.000000000000000e+00
    };

    std::vector<T>cf_weight(cf_mag.size() / 2, 2000);

    cf_weight[cf_weight.size() - 1] = 1;

    std::size_t deg = 59;

    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(2u*deg+1u, cf_freq, cf_mag, cf_weight);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2u*deg+1u, cf_freq, cf_mag, cf_weight);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2u*deg+1u, cf_freq, cf_mag, cf_weight);
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;
}

TYPED_TEST(firpm_issues_test, smallfir1) {

    using T = typename TestFixture::T;
    std::size_t degree = 10u;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, {0.0, 0.4, 0.6, 0.64, 0.69, 0.74, 0.79, 0.83, 0.88, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 1.0, 1.0, 1.0, 1.0});
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, {0.0, 0.4, 0.6, 0.64, 0.69, 0.74, 0.79, 0.83, 0.88, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 1.0, 1.0, 1.0, 1.0});
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, {0.0, 0.4, 0.6, 0.64, 0.69, 0.74, 0.79, 0.83, 0.88, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 1.0, 1.0, 1.0, 1.0});
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;
}

TYPED_TEST(firpm_issues_test, smallfir2) {
    using T = typename TestFixture::T;
    std::size_t degree = 6u;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u, {0.0, 0.4, 0.6, 0.6089285714285714, 0.6178571428571429, 0.6267857142857143, 0.6357142857142857, 0.6446428571428571, 0.6535714285714286, 0.6625, 0.6714285714285714, 0.6803571428571428, 0.6892857142857143, 0.6982142857142857, 0.7071428571428571, 0.7160714285714285, 0.725, 0.7339285714285714, 0.7428571428571429, 0.7517857142857143, 0.7607142857142857, 0.7696428571428571, 0.7785714285714285, 0.7874999999999999, 0.7964285714285714, 0.8053571428571428, 0.8142857142857143, 0.8232142857142857, 0.8321428571428571, 0.8410714285714285, 0.85, 0.8589285714285714, 0.8678571428571429, 0.8767857142857143, 0.8857142857142857, 0.8946428571428571, 0.9035714285714285, 0.9124999999999999, 0.9214285714285714, 0.9303571428571428, 0.9392857142857143, 0.9482142857142857, 0.9571428571428571, 0.9660714285714285, 0.975, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {87.35985578898647, 1000.0, 1029.7619047619048, 1059.5238095238094, 1089.2857142857144, 1119.047619047619, 1148.8095238095239, 1178.5714285714284, 1208.3333333333333, 1238.0952380952383, 1267.857142857143, 1297.6190476190475, 1327.3809523809523, 1357.142857142857, 1386.904761904762, 1416.6666666666667, 1446.4285714285716, 1476.1904761904761, 1505.952380952381, 1535.7142857142858, 1565.4761904761906, 1595.2380952380952, 1625.0});
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u, {0.0, 0.4, 0.6, 0.6089285714285714, 0.6178571428571429, 0.6267857142857143, 0.6357142857142857, 0.6446428571428571, 0.6535714285714286, 0.6625, 0.6714285714285714, 0.6803571428571428, 0.6892857142857143, 0.6982142857142857, 0.7071428571428571, 0.7160714285714285, 0.725, 0.7339285714285714, 0.7428571428571429, 0.7517857142857143, 0.7607142857142857, 0.7696428571428571, 0.7785714285714285, 0.7874999999999999, 0.7964285714285714, 0.8053571428571428, 0.8142857142857143, 0.8232142857142857, 0.8321428571428571, 0.8410714285714285, 0.85, 0.8589285714285714, 0.8678571428571429, 0.8767857142857143, 0.8857142857142857, 0.8946428571428571, 0.9035714285714285, 0.9124999999999999, 0.9214285714285714, 0.9303571428571428, 0.9392857142857143, 0.9482142857142857, 0.9571428571428571, 0.9660714285714285, 0.975, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {87.35985578898647, 1000.0, 1029.7619047619048, 1059.5238095238094, 1089.2857142857144, 1119.047619047619, 1148.8095238095239, 1178.5714285714284, 1208.3333333333333, 1238.0952380952383, 1267.857142857143, 1297.6190476190475, 1327.3809523809523, 1357.142857142857, 1386.904761904762, 1416.6666666666667, 1446.4285714285716, 1476.1904761904761, 1505.952380952381, 1535.7142857142858, 1565.4761904761906, 1595.2380952380952, 1625.0});
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u, {0.0, 0.4, 0.6, 0.6089285714285714, 0.6178571428571429, 0.6267857142857143, 0.6357142857142857, 0.6446428571428571, 0.6535714285714286, 0.6625, 0.6714285714285714, 0.6803571428571428, 0.6892857142857143, 0.6982142857142857, 0.7071428571428571, 0.7160714285714285, 0.725, 0.7339285714285714, 0.7428571428571429, 0.7517857142857143, 0.7607142857142857, 0.7696428571428571, 0.7785714285714285, 0.7874999999999999, 0.7964285714285714, 0.8053571428571428, 0.8142857142857143, 0.8232142857142857, 0.8321428571428571, 0.8410714285714285, 0.85, 0.8589285714285714, 0.8678571428571429, 0.8767857142857143, 0.8857142857142857, 0.8946428571428571, 0.9035714285714285, 0.9124999999999999, 0.9214285714285714, 0.9303571428571428, 0.9392857142857143, 0.9482142857142857, 0.9571428571428571, 0.9660714285714285, 0.975, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}, {87.35985578898647, 1000.0, 1029.7619047619048, 1059.5238095238094, 1089.2857142857144, 1119.047619047619, 1148.8095238095239, 1178.5714285714284, 1208.3333333333333, 1238.0952380952383, 1267.857142857143, 1297.6190476190475, 1327.3809523809523, 1357.142857142857, 1386.904761904762, 1416.6666666666667, 1446.4285714285716, 1476.1904761904761, 1505.952380952381, 1535.7142857142858, 1565.4761904761906, 1595.2380952380952, 1625.0});
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;
}

TYPED_TEST(firpm_issues_test, partition1) {

    using T = typename TestFixture::T;
    std::size_t degree = 201u;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u + 1u,
                            {0.0, 0.0126667, 0.0126667, 0.018, 0.018, 0.021, 0.112333, 1.0},
                            {1.0, 1.00055, 1.00055, 1.00111, 1.00111, 1.00151, 0.0, 0.0},
                            {1, 1, 1, 1});
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpm<T>(degree * 2u + 1u,
                            {0.0, 0.0126667, 0.0126667, 0.018, 0.018, 0.021, 0.112333, 1.0},
                            {1.0, 1.00055, 1.00055, 1.00111, 1.00111, 1.00151, 0.0, 0.0},
                            {1, 1, 1, 1});
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpm<T>(degree * 2u + 1u,
                            {0.0, 0.0126667, 0.0126667, 0.018, 0.018, 0.021, 0.112333, 1.0},
                            {1.0, 1.00055, 1.00055, 1.00111, 1.00111, 1.00151, 0.0, 0.0},
                            {1, 1, 1, 1});
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;
}

TYPED_TEST(firpm_issues_test, partition2) {

    using T = typename TestFixture::T;
    std::size_t degree = 50u;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(degree * 2u,
                            {0.0, 0.00566667, 0.00566667, 0.008, 0.008, 0.00833333, 0.112333, 1.0},
                            {1.0, 1.00011, 1.00011, 1.00022, 1.00022, 1.00024, 0.0, 0.0},
                            {1, 1, 1, 750});
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.q, 1e-2);

    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(degree * 2u,
                            {0.0, 0.00566667, 0.00566667, 0.008, 0.008, 0.00833333, 0.112333, 1.0},
                            {1.0, 1.00011, 1.00011, 1.00022, 1.00022, 1.00024, 0.0, 0.0},
                            {1, 1, 1, 750});
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output2.delta)/output1.delta), 2e-2);

    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(degree * 2u,
                            {0.0, 0.00566667, 0.00566667, 0.008, 0.008, 0.00833333, 0.112333, 1.0},
                            {1.0, 1.00011, 1.00011, 1.00022, 1.00022, 1.00024, 0.0, 0.0},
                            {1, 1, 1, 750});
    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);
    ASSERT_LE(pm::pmmath::fabs((output1.delta-output3.delta)/output1.delta), 2e-2);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;
}