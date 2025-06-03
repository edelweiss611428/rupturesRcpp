cmake_minimum_required(VERSION 3.30)
project(cpdL2)

set(CMAKE_CXX_STANDARD 20)

# Find and link BLAS and LAPACK
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

# Add the path to the Eigen library
include_directories(/Users/charles/Documents/code/armadillo-14.4.1/include)

add_executable(var main.cpp
        utils.cpp
        utils.h
        costL2.cpp
        costL2.h
        peltL2.cpp
        peltL2.h
)
target_link_libraries(var PRIVATE ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})