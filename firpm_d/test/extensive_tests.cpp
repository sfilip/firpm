#include <vector>
#include <fstream>
#include <chrono>
#include "firpm/barycentric.h"
#include "firpm/pm.h"
#include "firpm/cheby.h"
#include "gtest/gtest.h"

using testing::Types;

using types = testing::Types<double, long double>;

template<typename _T>
struct firpm_extensive_test : public testing::Test { using T = _T; };
TYPED_TEST_CASE(firpm_extensive_test, types);

template<typename T>
void printInfo(pmoutput_t<T>& output, double eps)
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

template<typename T>
void compareInfoRS(pmoutput_t<T>& output1, pmoutput_t<T>& output2, double eps)
{
    if(output1.q < eps)
    {
        std::cout << "Iteration reduction for final filter  RS: " << 1.0 - (double)output2.iter / output1.iter << std::endl;
    }
}

template<typename T>
void compareInfoAFP(pmoutput_t<T>& output1, pmoutput_t<T>& output2, double eps)
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
TYPED_TEST(firpm_extensive_test, extensive1)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TYPED_TEST(firpm_extensive_test, extensive2)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TYPED_TEST(firpm_extensive_test, extensive3)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(440, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TYPED_TEST(firpm_extensive_test, extensive4)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "Delta = " << output2.delta << std::endl;
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}


// Specification 2


TYPED_TEST(firpm_extensive_test, extensive5)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "Delta = " << output2.delta << std::endl;
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(400, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TYPED_TEST(firpm_extensive_test, extensive6)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(401, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive7)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(402, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive8)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(441, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive9)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(442, {0.0, 0.38, 0.45, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}





// Specification 3
TYPED_TEST(firpm_extensive_test, extensive10)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(150, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive11)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(160, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TYPED_TEST(firpm_extensive_test, extensive12)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(161, {0.0, 0.2, 0.4, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}


// Specification 4

TYPED_TEST(firpm_extensive_test, extensive13)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1000, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TYPED_TEST(firpm_extensive_test, extensive14)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2002, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TYPED_TEST(firpm_extensive_test, extensive15)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0}, 1e-2, 4u, 2u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2483, {0.0, 0.8, 0.81, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    //compareInfoRS(output1, output2, 1e-2);
    //compareInfoAFP(output1, output3, 1e-2);

}


// Specification 5

TYPED_TEST(firpm_extensive_test, extensive16)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 0.0001);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2002, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive17)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(4422, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);

}

TYPED_TEST(firpm_extensive_test, extensive18)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(4560, {0.0, 0.1, 0.105, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 6
TYPED_TEST(firpm_extensive_test, extensive19)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1400, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive20)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(3002, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive21)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0}, 1e-2, 4u, 2u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(4200, {0.0, 0.1, 0.105, 0.6, 0.605, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0}, 1e-2);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 7
TYPED_TEST(firpm_extensive_test, extensive22)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(3002, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive23)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2256, {0.0, 0.2, 0.205, 0.7, 0.705, 1.0}, {1.0, 1.0, 0.0, 0.0, 0.7, 0.7}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 8
TYPED_TEST(firpm_extensive_test, extensive24)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(600, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive25)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(362, {0.0, 0.3, 0.35, 0.7, 0.75, 1.0}, {1.0, 1.0, 0.0, 0.0, 1.0, 1.0}, {1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 9
TYPED_TEST(firpm_extensive_test, extensive26)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(200, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


TYPED_TEST(firpm_extensive_test, extensive27)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(202, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive28)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(250, {0.0, 0.4, 0.5, 0.7, 0.8, 1.0}, {0.0, 0.0, 1.0, 1.0, 0.0, 0.0}, {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 10
TYPED_TEST(firpm_extensive_test, extensive29)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(140, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive30)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 8u, 1u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(290, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 8);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


TYPED_TEST(firpm_extensive_test, extensive31)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(422, {0.0, 0.2, 0.25, 0.4, 0.45, 0.5, 0.55, 0.7, 0.75, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 11
TYPED_TEST(firpm_extensive_test, extensive32)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1600, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive33)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u, 2u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2300, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


TYPED_TEST(firpm_extensive_test, extensive34)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u, 1u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2414, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 1.0},
            {1.0, 1.0, 0.0, 0.0, 0.7, 0.7, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0}, 1e-2, 8u);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 12
TYPED_TEST(firpm_extensive_test, extensive35)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2000, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive36)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2800, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive37)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(3042, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive38)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1200, {0.0, 0.21, 0.215, 0.42, 0.425, 0.61, 0.615, 0.74, 0.745, 0.8, 0.805, 0.9, 0.905, 1.0},
            {1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0},
            {1.0, 10.0, 1.0, 10.0, 1.0, 10.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 13

TYPED_TEST(firpm_extensive_test, extensive39)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";

    auto output1 = firpm<T>(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1200, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive40)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 16u);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 8u, 1u);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1800, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0}, 1e-2, 8u);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive41)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
            {0.0, 0.0, 0.6, 0.6, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0, 1.0, 10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(774, {0.0, 0.24, 0.25, 0.39, 0.4, 0.57, 0.58, 0.7, 0.71, 0.8, 0.81, 0.92, 0.93, 1.0},
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
TYPED_TEST(firpm_extensive_test, extensive42)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(161, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}
TYPED_TEST(firpm_extensive_test, extensive43)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(201, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);


    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive44)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(223, {0.0, 0.4, 0.5, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}


// Specification 15
TYPED_TEST(firpm_extensive_test, extensive45)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(401, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive46)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(801, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive47)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1601, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive48)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1847, {0.0, 0.7, 0.71, 1.0}, {1.0, 1.0, 0.0, 0.0}, {1.0, 1.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 16
TYPED_TEST(firpm_extensive_test, extensive49)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(1801, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive50)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(2203, {0.0, 0.29, 0.3, 0.8, 0.81, 1.0},
            {0.0, 0.0, 1.0, 1.0, 0.0, 0.0},
            {10.0, 1.0, 10.0});
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";

    compareInfoRS(output1, output2, 1e-2);
}

// Type 3 & 4 Hilbert transformers
// Specification 17
TYPED_TEST(firpm_extensive_test, extensive51)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(4000, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);


    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive52)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(4001, {0.001, 0.999}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);


    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 18
TYPED_TEST(firpm_extensive_test, extensive53)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(100, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive54)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(101, {0.1, 0.9}, {1.0, 1.0}, {1.0}, filter_t::FIR_HILBERT);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);


    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Type 3 & 4 differentiators
// Specification 19
TYPED_TEST(firpm_extensive_test, extensive55)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(100, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive56)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(101, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive57)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(200, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive58)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(201, {0, 0.5, 0.55, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

// Specification 20
TYPED_TEST(firpm_extensive_test, extensive59)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(600, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}

TYPED_TEST(firpm_extensive_test, extensive60)
{
    using T = typename TestFixture::T;
    std::cout << "START Parks-McClellan with uniform initialization\n";
    auto output1 = firpm<T>(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output1, 1e-2);
    std::cout << "FINISH Parks-McClellan with uniform initialization\n";
    std::cout << "START Parks-McClellan with reference scaling\n";
    auto output2 = firpmRS<T>(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output2, 1e-2);
    std::cout << "FINISH Parks-McClellan with reference scaling\n";
    ASSERT_LT(output2.q, 1e-2);
    std::cout << "START Parks-McClellan with AFP\n";
    auto output3 = firpmAFP<T>(601, {0, 0.7, 0.71, 1.0}, {0.0, 1.0, 0.0, 0.0}, {1.0, 1.0}, filter_t::FIR_DIFFERENTIATOR);
    printInfo(output3, 1e-2);
    std::cout << "FINISH Parks-McClellan with AFP\n";
    ASSERT_LT(output3.q, 1e-2);

    compareInfoRS(output1, output2, 1e-2);
    compareInfoAFP(output1, output3, 1e-2);
}
