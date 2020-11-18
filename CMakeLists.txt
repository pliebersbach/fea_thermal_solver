cmake_minimum_required(VERSION 3.0.0)
project(myTest VERSION 0.1.0)

set(CMAKE_BUILD_TYPE debug)
#set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

add_executable(matrixTest matrix_test.cpp math_utilities.h math_utilities.cpp)
add_executable(myTest main.cpp elements.cpp mesh.cpp stiffnessMatrix.cpp math_utilities.cpp elements.h mesh.h stiffnessMatrix.h math_utilities.h)
add_executable(stiff stiffness_test.cpp elements.cpp mesh.cpp stiffnessMatrix.cpp math_utilities.cpp mesh.h elements.h stiffnessMatrix.h math_utilities.h)
add_executable(plate plate.cpp elements.cpp mesh.cpp stiffnessMatrix.cpp math_utilities.cpp mesh.h elements.h stiffnessMatrix.h math_utilities.h)
add_executable(slab slab.cpp elements.cpp mesh.cpp stiffnessMatrix.cpp math_utilities.cpp mesh.h elements.h stiffnessMatrix.h math_utilities.h)