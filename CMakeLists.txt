cmake_minimum_required (VERSION 3.16)
set(PROJECT_NAME_STR firpmlib)
project(${PROJECT_NAME_STR})
# The project version number.
set(VERSION_MAJOR   0   CACHE STRING "Project major version number.")
set(VERSION_MINOR   2   CACHE STRING "Project minor version number.")
set(VERSION_PATCH   1   CACHE STRING "Project patch version number.")
mark_as_advanced(VERSION_MAJOR VERSION_MINOR VERSION_PATCH)
set(doxy_main_page ${PROJECT_SOURCE_DIR}/README.md)

#-----------------------------------------------------------------------
# dependencies testing
#-----------------------------------------------------------------------
list (INSERT CMAKE_MODULE_PATH 0 ${PROJECT_SOURCE_DIR}/cmake)

find_package(OpenMP REQUIRED)
find_package(Eigen3 3.3 REQUIRED)
find_package(GMP 6.0.0)
find_package(MPFR 4.0.0)
# find_package(Doxygen)
find_package(Python COMPONENTS Interpreter Development)
find_package(pybind11 CONFIG)

if( MPFR_FOUND AND GMP_FOUND)
    add_definitions(-DHAVE_MPFR)
endif( MPFR_FOUND AND GMP_FOUND)

#-----------------------------------------------------------------------
# build configuration
#-----------------------------------------------------------------------
set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

#-----------------------------------------------------------------------
# common include directories
#-----------------------------------------------------------------------
if( MPFR_FOUND AND GMP_FOUND )
    set(MPREAL_HEADER mpreal.h)
    find_path(MPREAL_INCLUDE_DIR NAMES ${MPREAL_HEADER})
    if (NOT MPREAL_INCLUDE_DIR)
        message(FATAL_ERROR "\t ${MPREAL_HEADER} header not found.")
    endif()
    message("${MPREAL_HEADER} found in ${MPREAL_INCLUDE_DIR}")
    set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include ${MPFR_INCLUDES} ${MPREAL_INCLUDE_DIR})
else()
    set(COMMON_INCLUDES ${PROJECT_SOURCE_DIR}/include)
endif( MPFR_FOUND AND GMP_FOUND )

add_subdirectory(${PROJECT_SOURCE_DIR}/src)
if(DOXYGEN_FOUND)
    add_subdirectory(${PROJECT_SOURCE_DIR}/doc)
endif()

add_subdirectory(${PROJECT_SOURCE_DIR}/test)

#-----------------------------------------------------------------------
# Python bindings
#-----------------------------------------------------------------------

if( pybind11_FOUND )
	pybind11_add_module(pyfirpm py/bindings.cpp)
    include_directories(${COMMON_INCLUDES} ${EIGEN3_INCLUDE_DIRS})
	target_link_libraries(pyfirpm PRIVATE firpm)
endif()


if( pybind11_FOUND AND Python_FOUND )
	add_test(NAME PythonTests COMMAND ${Python3_EXECUTABLE} -m pytest test/firpm_tests.py)
endif()
