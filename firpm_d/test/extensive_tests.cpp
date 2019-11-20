#include <vector>
#include <fstream>
#include <chrono>
#include "firpm/barycentric.h"
#include "firpm/pm.h"
#include "firpm/cheby.h"
#include "gtest/gtest.h"

void printInfo(pmoutput_t<double>& output, double eps)
{
	if(output.q < eps)
	{
		std::cout << "Final delta     = " << output.delta << std::endl;
		std::cout << "Iteration count = " << output.iter << std::endl;
	}
	else
	{
		std::cout << "Iteration count = NC\n";
	}
}


void compareInfoRS(pmoutput_t<double>& output1, pmoutput_t<double>& output2, double eps)
{
    if(output1.q < eps)
    {
        std::cout << "Iteration reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    }
}

void compareInfoAFP(pmoutput_t<double>& output1, pmoutput_t<double>& output2, double eps)
{
    if(output1.q < eps)
    {
        std::cout << "Iteration reduction for final filter AFP: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    }
}

// Specifications: 17
// Filters: 51

// Type 1 filters
// Specification 1
TEST(firpm_extensive_test, extensive1)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TEST(firpm_extensive_test, extensive2)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TEST(firpm_extensive_test, extensive3)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TEST(firpm_extensive_test, extensive4)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "Delta = " << output2.delta << std::endl;
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}


// Specification 2


TEST(firpm_extensive_test, extensive5)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "Delta = " << output2.delta << std::endl;
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TEST(firpm_extensive_test, extensive6)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive7)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive8)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive9)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}





// Specification 3
TEST(firpm_extensive_test, extensive10)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive11)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TEST(firpm_extensive_test, extensive12)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}


// Specification 4

TEST(firpm_extensive_test, extensive13)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TEST(firpm_extensive_test, extensive14)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TEST(firpm_extensive_test, extensive15)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0}, 1e-2, 4u, 2u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    //compareInfoRS(output1, output2, 1e-2);
    //compareInfoAFP(output1, output3, 1e-2);

}


// Specification 5

TEST(firpm_extensive_test, extensive16)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.0001);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive17)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TEST(firpm_extensive_test, extensive18)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 6
TEST(firpm_extensive_test, extensive19)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive20)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive21)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0}, 1e-2, 4u, 2u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0}, 1e-2);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 7
TEST(firpm_extensive_test, extensive22)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive23)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 8
TEST(firpm_extensive_test, extensive24)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive25)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 9
TEST(firpm_extensive_test, extensive26)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


TEST(firpm_extensive_test, extensive27)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive28)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 10
TEST(firpm_extensive_test, extensive29)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive30)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 8u, 1u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 8);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


TEST(firpm_extensive_test, extensive31)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 11
TEST(firpm_extensive_test, extensive32)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive33)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u, 2u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


TEST(firpm_extensive_test, extensive34)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u, 1u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 12
TEST(firpm_extensive_test, extensive35)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive36)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive37)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive38)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 13

TEST(firpm_extensive_test, extensive39)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<double>(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive40)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 16u);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 8u, 1u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 8u);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive41)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Type 2 filters
// Specification 14
TEST(firpm_extensive_test, extensive42)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}
TEST(firpm_extensive_test, extensive43)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);


    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive44)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 15
TEST(firpm_extensive_test, extensive45)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive46)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive47)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive48)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 16
TEST(firpm_extensive_test, extensive49)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
}

TEST(firpm_extensive_test, extensive50)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
}

// Type 3 & 4 Hilbert transformers
// Specification 17
TEST(firpm_extensive_test, extensive51)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);


    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive52)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);


    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 18
TEST(firpm_extensive_test, extensive53)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive54)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);


    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Type 3 & 4 differentiators
// Specification 19
TEST(firpm_extensive_test, extensive55)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive56)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive57)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive58)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 20
TEST(firpm_extensive_test, extensive59)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TEST(firpm_extensive_test, extensive60)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<double>(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<double>(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<double>(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}
