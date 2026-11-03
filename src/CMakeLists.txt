cmake_minimum_required(VERSION 3.16)
set(PROJECT_NAME_STR firpmlib-build)
project(${PROJECT_NAME_STR})

set(CMAKE_BUILD_TYPE Release)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

file(GLOB PROJECT_SRC_FILES ${PROJECT_SOURCE_DIR}/*.cpp)
set(PROJECT_INCLUDE_DIR ${PROJECT_SOURCE_DIR}/../include/)
if( MPFR_FOUND AND GMP_FOUND )
    include_directories (${PROJECT_INCLUDE_DIR} ${MPFR_INCLUDES} ${MPREAL_INCLUDE_DIR})
else()
    include_directories (${PROJECT_INCLUDE_DIR})
endif( MPFR_FOUND AND GMP_FOUND )

add_library(firpm ${PROJECT_SRC_FILES})
target_link_libraries(firpm PUBLIC Eigen3::Eigen)

if( MPFR_FOUND AND GMP_FOUND )
    target_link_libraries(firpm PRIVATE ${GMP_LIBRARIES} ${MPFR_LIBRARIES} OpenMP::OpenMP_CXX)
    target_include_directories(firpm PUBLIC ${PROJECT_INCLUDE_DIR} ${MPFR_INCLUDES} ${MPREAL_INCLUDE_DIR})
else()
    target_link_libraries(firpm PRIVATE OpenMP::OpenMP_CXX)
    target_include_directories(firpm PUBLIC ${PROJECT_INCLUDE_DIR})
endif()

install(TARGETS firpm DESTINATION lib)
install(DIRECTORY ${PROJECT_INCLUDE_DIR} DESTINATION include)
