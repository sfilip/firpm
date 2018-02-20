cmake_minimum_required(VERSION 3.0)
set(PROJECT_NAME_STR firpmlib-build)
project(${PROJECT_NAME_STR})

set(CMAKE_BUILD_TYPE Release)
# use C++11 standard (might switch to a newer standard in the future)
add_definitions("-std=c++11")
# Enforce export of no symbols by default (if this is supported
# by the compiler). The goal is to make the visibility behavior
# of GCC consistent to that of MSVC
set(CMAKE_C_VISIBILITY_PRESET hidden)
set(CMAKE_CXX_VISIBILITY_PRESET hidden)

set(CMAKE_INCLUDE_CURRENT_DIR ON)

file(GLOB PROJECT_SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../include)
include_directories(${PROJECT_INCLUDE_DIR})

# create the shared library
add_library(firpm SHARED ${PROJECT_SRC_FILES})

# generate export symbol header automatically
include (GenerateExportHeader)
generate_export_header(firpm
    BASE_NAME firpm
    EXPORT_MACRO_NAME FIRPM_EXPORT
    EXPORT_FILE_NAME ${CMAKE_BINARY_DIR}/exports/firpm_EXPORT.h
    STATIC_DEFINE FIRPM_EXPORT_BUILT_AS_STATIC)
# add the directory with the generated headers to the include path
target_include_directories(firpm PUBLIC ${CMAKE_BINARY_DIR}/exports)


target_link_libraries(firpm gomp)


install(TARGETS firpm DESTINATION lib)
install(DIRECTORY ${PROJECT_INCLUDE_DIR}/firpm DESTINATION include)
