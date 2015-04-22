cmake_minimum_required (VERSION 2.8.8)
set(PROJECT_NAME_STR firpmlib)
project(${PROJECT_NAME_STR})
# The project version number.
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   1   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)
set(doxy_main_page ${PROJECT_SOURCE_DIR}/README.md)

#----------------------------------------
# dependencies testing
#----------------------------------------
list (INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)

find_package(OpenMP REQUIRED)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()
find_package(Eigen3 REQUIRED)
find_package(GTest)
find_package(Doxygen)


#----------------------------------------
# build configuration
#----------------------------------------
set(CMAKE_BUILD_TYPE Release)
add_definitions(-std=c++11)

#----------------------------------------
# common include directories
#----------------------------------------
set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
if(DOXYGEN_FOUND)
    add_subdirectory(${PROJECT_SOURCE_DIR}/doc)
endif()

#----------------------------------------
# tests
#----------------------------------------
if(GTEST_FOUND)
    enable_testing()
    set(PROJECT_TEST_SCALING ${PROJECT_NAME_STR}_scaling_test)
    set(PROJECT_TEST_EXTENSIVE ${PROJECT_NAME_STR}_extensive_test)
    set(CMAKE_CTEST_COMMAND ctest -V)
    include_directories(${GTEST_INCLUDE_DIR} ${COMMON_INCLUDES})


    set(TEST_SRC_SCALING ${PROJECT_SOURCE_DIR}/test/scaling_tests.cpp)
    set(TEST_SRC_EXTENSIVE ${PROJECT_SOURCE_DIR}/test/extensive_tests.cpp)


    add_executable(${PROJECT_TEST_SCALING} ${TEST_SRC_SCALING})
    add_executable(${PROJECT_TEST_EXTENSIVE} ${TEST_SRC_EXTENSIVE})


    target_link_libraries(${PROJECT_TEST_SCALING}
        gtest
        gtest_main
        pthread
        firpm
        gomp
    )

    target_link_libraries(${PROJECT_TEST_EXTENSIVE}
        gtest
        gtest_main
        pthread
        firpm
        gomp
    )


    add_test(ScalingTests ${PROJECT_TEST_SCALING})
    add_test(ExtensiveTests ${PROJECT_TEST_EXTENSIVE})
endif()
