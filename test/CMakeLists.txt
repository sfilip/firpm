include(FetchContent)

FetchContent_Declare(
  googletest
  GIT_REPOSITORY https://github.com/google/googletest.git
  GIT_TAG        v1.14.0
)
FetchContent_MakeAvailable(googletest)
add_library(GTest::GTest INTERFACE IMPORTED)
target_link_libraries(GTest::GTest INTERFACE gtest_main)

enable_testing()
set(PROJECT_TEST_SCALING ${PROJECT_NAME_STR}_scaling_test)
set(PROJECT_TEST_EXTENSIVE ${PROJECT_NAME_STR}_extensive_test)

set(CMAKE_CTEST_COMMAND ctest -V)
include_directories(${GTEST_INCLUDE_DIR} ${COMMON_INCLUDES} ${EIGEN3_INCLUDE_DIRS})

set(TEST_SRC_SCALING scaling_tests.cpp)
set(TEST_SRC_EXTENSIVE ${PROJECT_SOURCE_DIR}/test/extensive_tests.cpp)

add_executable(${PROJECT_TEST_SCALING} ${TEST_SRC_SCALING})
add_executable(${PROJECT_TEST_EXTENSIVE} ${TEST_SRC_EXTENSIVE})


if( MPFR_FOUND AND GMP_FOUND )
    target_link_libraries(${PROJECT_TEST_SCALING}
        GTest::GTest pthread
        OpenMP::OpenMP_CXX
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES} firpm
    )

    target_link_libraries(${PROJECT_TEST_EXTENSIVE}
        GTest::GTest pthread
        OpenMP::OpenMP_CXX
        ${GMP_LIBRARIES}
        ${MPFR_LIBRARIES} firpm
    )

else()
    target_link_libraries(${PROJECT_TEST_SCALING}
        GTest::GTest pthread
        OpenMP::OpenMP_CXX
        firpm
    )

    target_link_libraries(${PROJECT_TEST_EXTENSIVE}
        GTest::GTest pthread
        OpenMP::OpenMP_CXX
        firpm
    )

endif()



add_test(ScalingTests ${PROJECT_TEST_SCALING})
add_test(ExtensiveTests ${PROJECT_TEST_EXTENSIVE})