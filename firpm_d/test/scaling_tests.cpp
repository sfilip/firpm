#include <vector>
#include <fstream>
#include <chrono>
#include "firpm/barycentric.h"
#include "firpm/pm.h"
#include "firpm/cheby.h"
#include "gtest/gtest.h"


// example of how to use the exchange method directly
TEST(firpm_scaling_test, lowpass100a)
{

    std::vector<Band> freqBands(2);

    freqBands[0].start = 0;
    freqBands[0].stop = M_PI * 0.4;
    freqBands[0].weight = [] (BandSpace, double) -> double {return 1; };
    freqBands[0].space = BandSpace::FREQ;
    freqBands[0].amplitude = [](BandSpace, double) -> double { return 1; };

    freqBands[1].start = M_PI * 0.5;
    freqBands[1].stop = M_PI;
    freqBands[1].weight = [] (BandSpace, double) -> double {return 1; };
    freqBands[1].space = BandSpace::FREQ;
    freqBands[1].amplitude = [](BandSpace, double) -> double { return 0; };



    auto start = std::chrono::steady_clock::now();

    std::size_t degree = 50;
    std::vector<double> a;
    double finalDelta;
    std::vector<Band> chebyBands;
    std::vector<double> omega(degree + 2u);
    std::vector<double> x(degree + 2u);
    initUniformExtremas(omega, freqBands);
    applyCos(x, omega);
    bandConversion(chebyBands, freqBands, ConversionDirection::FROMFREQ);

    PMOutput output = exchange(x, chebyBands);
    ASSERT_LT(output.Q, 0.1e-1);

    for(std::size_t counter = 0; counter < 1; ++counter) {
        std::vector<double> newX;
        referenceScaling(newX, chebyBands, freqBands, 2 * degree + 2, output.x, chebyBands,
                freqBands);
        degree = 2 * degree;

        output = exchange(newX, chebyBands);
        ASSERT_LT(output.Q, 0.1e-1);
    }
        auto stop  = std::chrono::steady_clock::now();
        double elapsedTime = std::chrono::duration_cast<
            std::chrono::duration<double>>(stop - start).count();
        std::cout << "Elapsed time = " << elapsedTime << std::endl;

}

// Filters appearing inside Section 1 & Section 4 (minus the degree 53348 one)
TEST(firpm_scaling_test, lowpass50)
{

    std::vector<double> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {1.0, 1.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto start = std::chrono::steady_clock::now();
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    auto stop  = std::chrono::steady_clock::now();
    double elapsedTime = std::chrono::duration_cast<
        std::chrono::duration<double>>(stop - start).count();
    std::cout << "Elapsed time = " << elapsedTime << std::endl;

    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    start = std::chrono::steady_clock::now();
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    stop  = std::chrono::steady_clock::now();
    elapsedTime = std::chrono::duration_cast<
        std::chrono::duration<double>>(stop - start).count();
    std::cout << "Elapsed time = " << elapsedTime << std::endl;

    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    start = std::chrono::steady_clock::now();
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);
    stop  = std::chrono::steady_clock::now();
    elapsedTime = std::chrono::duration_cast<
        std::chrono::duration<double>>(stop - start).count();
    std::cout << "Elapsed time = " << elapsedTime << std::endl;

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);


    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_scaling_test, lowpass80)
{

    std::vector<double> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {1.0, 1.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_scaling_test, lowpass100)
{

    std::vector<double> f = {0.0, 0.4, 0.5, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {1.0, 1.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}



TEST(firpm_scaling_test, bandstop50)
{
    std::vector<double> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<double> w = {1.0, 1.0, 1.0};

    std::size_t degree = 50;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_scaling_test, bandstop80)
{
    std::vector<double> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<double> w = {1.0, 1.0, 1.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_scaling_test, bandstop100)
{
    std::vector<double> f = {0.0, 0.2, 0.3, 0.5, 0.6, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0, 1.0, 1.0};
    std::vector<double> w = {1.0, 1.0, 1.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}


TEST(firpm_scaling_test, combfir)
{
    std::vector<double> f = {0.0, 0.99, 1.0, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {1.0, 1.0};

    std::size_t degree = 520;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}


TEST(firpm_lebesgue_test, lowpass500)
{

    std::vector<double> f = {0.0, 0.49, 0.5, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {1.0, 10.0};

    std::size_t degree = 500;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}


TEST(firpm_lebesgue_test, lowpass1000)
{
    std::vector<double> f = {0.0, 0.49, 0.5, 1.0};
    std::vector<double> a = {1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {1.0, 10.0};

    std::size_t degree = 1000;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, bandpass60)
{

    std::vector<double> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<double> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {10.0, 1.0, 5.0};

    std::size_t degree = 60;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}



TEST(firpm_lebesgue_test, bandpass70)
{

    std::vector<double> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<double> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {10.0, 1.0, 5.0};

    std::size_t degree = 70;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}


TEST(firpm_lebesgue_test, bandpass80)
{

    std::vector<double> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<double> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {10.0, 1.0, 5.0};

    std::size_t degree = 80;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}


TEST(firpm_lebesgue_test, bandpass100)
{

    std::vector<double> f = {0.0, 0.15, 0.25, 0.6, 0.7, 1.0};
    std::vector<double> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {10.0, 1.0, 5.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, multiband100)
{

    std::vector<double> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<double> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 100;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, multiband200)
{

    std::vector<double> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<double> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 200;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, multiband300)
{

    std::vector<double> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<double> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 300;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}

TEST(firpm_lebesgue_test, multiband600)
{

    std::vector<double> f = {0.0, 0.18, 0.2, 0.4, 0.42, 0.55, 0.57, 0.65, 0.67, 0.75, 0.77, 0.85, 0.87, 1.0};
    std::vector<double> a = {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0};
    std::vector<double> w = {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0};

    std::size_t degree = 600;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(degree * 2u, f, a, w, 0.00001);
    std::cout << "Final Delta     = " << output1.delta << std::endl;
    std::cout << "Iteration count = " << output1.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    //ASSERT_LT(output1.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(degree * 2u, f, a, w, 0.00001);
    std::cout << "Final Delta     = " << output2.delta << std::endl;
    std::cout << "Iteration count = " << output2.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);

    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(degree * 2u, f, a, w, 0.00001);

    std::cout << "Final Delta     = " << output3.delta << std::endl;
    std::cout << "Iteration count = " << output3.iter  << std::endl;
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    std::cout << "Iteration count reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    std::cout << "Iteration count reduction for final filter AFP: " << 1.0 - (double)output3.iter / output1.iter << std::endl;

}
