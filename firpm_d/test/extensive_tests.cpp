#include <vector>
#include <fstream>
#include <chrono>
#include "firpm/barycentric.h"
#include "firpm/pm.h"
#include "firpm/cheby.h"
#include "gtest/gtest.h"

void printInfo(PMOutput& output, double eps)
{
	if(output.Q < eps)
	{
		std::cout << "Final delta     = " << output.delta << std::endl;
		std::cout << "Iteration count = " << output.iter << std::endl;
	}
	else
	{
		std::cout << "Iteration count = NC\n";
	}
}


void compareInfoRS(PMOutput& output1, PMOutput& output2, double eps)
{
    if(output1.Q < eps)
    {
        std::cout << "Iteration reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    }
}

void compareInfoAFP(PMOutput& output1, PMOutput& output2, double eps)
{
    if(output1.Q < eps)
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
    PMOutput output1 = firpm(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}

TEST(firpm_extensive_test, extensive2)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}

TEST(firpm_extensive_test, extensive3)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}

TEST(firpm_extensive_test, extensive4)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "Delta = " << output2.delta << std::endl;
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}


// Specification 2


TEST(firpm_extensive_test, extensive5)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "Delta = " << output2.delta << std::endl;
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}

TEST(firpm_extensive_test, extensive6)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive7)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive8)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive9)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}





// Specification 3
TEST(firpm_extensive_test, extensive10)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive11)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}

TEST(firpm_extensive_test, extensive12)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}


// Specification 4

TEST(firpm_extensive_test, extensive13)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}

TEST(firpm_extensive_test, extensive14)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}

TEST(firpm_extensive_test, extensive15)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0}, 0.01, 2u, 4);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    //compareInfoRS(output1, output2, 0.01);
    //compareInfoAFP(output1, output3, 0.01);

}


// Specification 5

TEST(firpm_extensive_test, extensive16)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.0001);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive17)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.01);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);

}

TEST(firpm_extensive_test, extensive18)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.01);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.01);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

// Specification 6
TEST(firpm_extensive_test, extensive19)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive20)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive21)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0}, 0.01, 2, 4u);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0}, 0.01);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


// Specification 7
TEST(firpm_extensive_test, extensive22)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive23)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

// Specification 8
TEST(firpm_extensive_test, extensive24)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive25)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


// Specification 9
TEST(firpm_extensive_test, extensive26)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


TEST(firpm_extensive_test, extensive27)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive28)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


// Specification 10
TEST(firpm_extensive_test, extensive29)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive30)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0}, 0.01, 1, 8);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0}, 0.01, 8);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


TEST(firpm_extensive_test, extensive31)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

// Specification 11
TEST(firpm_extensive_test, extensive32)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive33)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 0.01, 8);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 0.01, 2, 8);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 0.01, 8);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


TEST(firpm_extensive_test, extensive34)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 0.01, 8);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 0.01, 1, 8);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 0.01, 8);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

// Specification 12
TEST(firpm_extensive_test, extensive35)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive36)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive37)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive38)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


// Specification 13

TEST(firpm_extensive_test, extensive39)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";

    PMOutput output1 = firpm(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive40)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 0.01, 16);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 0.01, 1, 8);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 0.01, 8);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive41)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


// Type 2 filters
// Specification 14
TEST(firpm_extensive_test, extensive42)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}
TEST(firpm_extensive_test, extensive43)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);


    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive44)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}


// Specification 15
TEST(firpm_extensive_test, extensive45)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive46)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive47)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive48)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

// Specification 16
TEST(firpm_extensive_test, extensive49)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
}

TEST(firpm_extensive_test, extensive50)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 0.01);
}

// Type 3 & 4 Hilbert transformers
// Specification 17
TEST(firpm_extensive_test, extensive51)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);


    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive52)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);


    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

// Specification 18
TEST(firpm_extensive_test, extensive53)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive54)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, ftype::FIR_HILBERT);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);


    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

// Type 3 & 4 differentiators
// Specification 19
TEST(firpm_extensive_test, extensive55)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive56)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive57)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive58)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

// Specification 20
TEST(firpm_extensive_test, extensive59)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}

TEST(firpm_extensive_test, extensive60)
{
    std::cout << "START Parks-McClellan with uniform initialization\n";
    PMOutput output1 = firpm(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output1, 0.01);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    PMOutput output2 = firpmRS(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output2, 0.01);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.Q, 0.1e-1);
    std::cout << "START Parks-McClellan with AFP\n";
    PMOutput output3 = firpmAFP(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, ftype::FIR_DIFFERENTIATOR);
    printInfo(output3, 0.01);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.Q, 0.1e-1);

    compareInfoRS(output1, output2, 0.01);
    compareInfoAFP(output1, output3, 0.01);
}
