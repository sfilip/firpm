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
# Eigen's config only accepts requests for its own major version.
find_package(Eigen3 3.3 QUIET)
if(NOT Eigen3_FOUND)
    find_package(Eigen3 5 REQUIRED)
endif()
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
set(CMAKE_CXX_STANDARD 20)
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
	# The extension is the private half of the firpm package. Both halves
	# are assembled under python/ in the build tree, so it can be imported
	# from there without installing.
	set(FIRPM_PY_DIR ${CMAKE_CURRENT_BINARY_DIR}/python/firpm)
	pybind11_add_module(_firpm py/bindings.cpp)
    include_directories(${COMMON_INCLUDES})
	target_link_libraries(_firpm PRIVATE firpm)
	# The generator expression keeps multi-config generators from adding a
	# per-configuration subdirectory.
	set_target_properties(_firpm PROPERTIES
		LIBRARY_OUTPUT_DIRECTORY $<1:${FIRPM_PY_DIR}>)
	target_compile_definitions(_firpm PRIVATE
		FIRPM_VERSION="${VERSION_MAJOR}.${VERSION_MINOR}.${VERSION_PATCH}")
	set(FIRPM_PY_FILES __init__.py _firpm.pyi py.typed)
	foreach(file ${FIRPM_PY_FILES})
		configure_file(py/firpm/${file} ${FIRPM_PY_DIR}/${file} COPYONLY)
	endforeach()
endif()


if( pybind11_FOUND AND Python_FOUND )
	# No cache, so that running the tests leaves nothing in the source tree.
	add_test(NAME PythonTests
		COMMAND ${Python_EXECUTABLE} -m pytest -p no:cacheprovider
			${PROJECT_SOURCE_DIR}/test/firpm_tests.py)
	set_tests_properties(PythonTests PROPERTIES
		ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_BINARY_DIR}/python")
endif()
