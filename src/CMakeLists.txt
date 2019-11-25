cmake_minimum_required(VERSION 3.12)
set(PROJECT_NAME_STR firpmlib-build)
project(${PROJECT_NAME_STR})

set(CMAKE_BUILD_TYPE Release)
add_definitions("-std=c++11")

file(GLOB PROJECT_SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../include)
add_definitions     (${EIGEN3_DEFINITIONS})
include_directories (${PROJECT_INCLUDE_DIR} ${EIGEN3_INCLUDE_DIRS})

add_library(firpm SHARED ${PROJECT_SRC_FILES})

if( MPFR_FOUND AND GMP_FOUND )
    target_link_libraries(firpm PRIVATE gmp mpfr OpenMP::OpenMP_CXX)
else()
    target_link_libraries(firpm PRIVATE OpenMP::OpenMP_CXX)
endif()

install(TARGETS firpm DESTINATION lib)
install(DIRECTORY ${PROJECT_INCLUDE_DIR}/firpm DESTINATION include)