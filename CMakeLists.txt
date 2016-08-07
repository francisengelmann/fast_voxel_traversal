cmake_minimum_required(VERSION 2.6)
project(FastVoxelTraversal)

add_definitions(-Wall -std=c++0x)

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

add_executable(fast_voxel_traversal main.cpp)
target_link_libraries(fast_voxel_traversal)

install(TARGETS fast_voxel_traversal RUNTIME DESTINATION bin)
