#include <vector>
#include <fstream>
#include <chrono>
#include "firpm/util.h"
#include "firpm/barycentric.h"
#include "firpm/pm.h"
#include "firpm/cheby.h"
#include "firpm/band.h"
#include "firpm/eigenvalue.h"
#include "gtest/gtest.h"

// it takes quite a while (several hours) to generate the next filter, so run at your own risk
/*TEST(firpm_scaling_test, hugefilter)
{
    std::vector<mpfr::mpreal> f = {0.0, 1.0 / 8192, 3.0 / 8192, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {1.0, 1.0};

    std::size_t degree = 53248;
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u + 1u, f, a, w, 0.01, 2u, 4);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);
}*/




// example of how to use the exchange method directly
TEST(firpm_scaling_test, lowpass100a)
{
    using mpfr::mpreal;
    std::size_t prec = 165ul;
    mpreal::set_default_prec(prec);
    mpfr::mpreal pi = mpfr::const_pi(prec);

    std::vector<Band> freqBands(2);

    freqBands[0].start = 0;
    freqBands[0].stop = pi * 0.4;
    freqBands[0].weight = [] (BandSpace, mpfr::mpreal) -> mpfr::mpreal {return mpfr::mpreal(1); };
    freqBands[0].space = BandSpace::FREQ;
    freqBands[0].amplitude = [](BandSpace, mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(1); };

    freqBands[1].start = pi * 0.5;
    freqBands[1].stop = pi;
    freqBands[1].weight = [] (BandSpace, mpfr::mpreal) -> mpfr::mpreal {return mpfr::mpreal(1); };
    freqBands[1].space = BandSpace::FREQ;
    freqBands[1].amplitude = [](BandSpace, mpfr::mpreal) -> mpfr::mpreal { return mpfr::mpreal(0); };



    auto start = std::chrono::steady_clock::now();

    std::size_t degree = 50;
    std::vector<mpfr::mpreal> a;
    mpfr::mpreal finalDelta;
    std::vector<Band> chebyBands;
    std::vector<mpfr::mpreal> omega(degree + 2u);
    std::vector<mpfr::mpreal> x(degree + 2u);
    initUniformExtremas(omega, freqBands, prec);
    applyCos(x, omega);
    bandConversion(chebyBands, freqBands, ConversionDirection::FROMFREQ);

    PMOutput output = exchange(x, chebyBands);
    ASSERT_LT(output.Q.toDouble(), 0.1e-1);

    for(std::size_t counter = 0; counter < 1; ++counter) {
        std::vector<mpfr::mpreal> newX;
        referenceScaling(newX, chebyBands, freqBands, 2 * degree + 2, output.x, chebyBands,
                freqBands, prec);
        degree = 2 * degree;

        output = exchange(newX, chebyBands);
        ASSERT_LT(output.Q.toDouble(), 0.1e-1);
    }
        auto stop  = std::chrono::steady_clock::now();
        mpfr::mpreal elapsedTime = std::chrono::duration_cast<
            std::chrono::duration<mpfr::mpreal>>(stop - start).count();
        std::cout << "Elapsed time = " << elapsedTime << std::endl;

}

// Filters appearing inside Section 1 & Section 4
TEST(firpm_scaling_test, lowpass50)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {1.0, 1.0};

    std::size_t degree = 50;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}

TEST(firpm_scaling_test, lowpass80)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {1.0, 1.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;
}

TEST(firpm_scaling_test, lowpass100)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {1.0, 1.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}



TEST(firpm_scaling_test, bandstop50)
{
    std::vector<mpfr::mpreal> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<mpfr::mpreal> w = {1.0, 1.0, 1.0};

    std::size_t degree = 50;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}

TEST(firpm_scaling_test, bandstop80)
{
    std::vector<mpfr::mpreal> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<mpfr::mpreal> w = {1.0, 1.0, 1.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;


}

TEST(firpm_scaling_test, bandstop100)
{
    std::vector<mpfr::mpreal> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<mpfr::mpreal> w = {1.0, 1.0, 1.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}


TEST(firpm_scaling_test, combfir)
{
    std::vector<mpfr::mpreal> f = {0.0, 0.99, 1.0, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {1.0, 1.0};

    std::size_t degree = 520;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}


TEST(firpm_lebesgue_test, lowpass500)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.49, 0.5, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {1.0, 10.0};

    std::size_t degree = 500;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}


TEST(firpm_lebesgue_test, lowpass1000)
{
    std::vector<mpfr::mpreal> f = {0.0, 0.49, 0.5, 1.0};
    std::vector<mpfr::mpreal> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {1.0, 10.0};

    std::size_t degree = 1000;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, bandpass60)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<mpfr::mpreal> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {10.0, 1.0, 5.0};

    std::size_t degree = 60;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}



TEST(firpm_lebesgue_test, bandpass70)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<mpfr::mpreal> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {10.0, 1.0, 5.0};

    std::size_t degree = 70;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}


TEST(firpm_lebesgue_test, bandpass80)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<mpfr::mpreal> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {10.0, 1.0, 5.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}


TEST(firpm_lebesgue_test, bandpass100)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<mpfr::mpreal> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {10.0, 1.0, 5.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, multiband50)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<mpfr::mpreal> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 50;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, multiband100)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<mpfr::mpreal> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, multiband200)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<mpfr::mpreal> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 200;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, multiband300)
{

    std::vector<mpfr::mpreal> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<mpfr::mpreal> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<mpfr::mpreal> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 300;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q.toDouble(), 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q.toDouble(), 0.1e-1);

    std::cout << "Iteration count reduction for final filter: " << 1.0 - (mpfr::mpreal)output2.iter / output1.iter << std::endl;

}
