cmake_minimum_required(VERSION 3.0.0)
project(myTest VERSION 0.1.0)

set(CMAKE_BUILD_TYPE debug)
set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

#add_subdirectory(${..}/src)
#set(HEADER_FILES ${../include}/elements.h ${../include}/mesh.h ${../include}/math_utilities.h ${../include}/stiffnessMatrix.h)

#add_library(mylib ${../src}elements.cpp ${HEADER_FILES})
#target_include_directories(mylib PRIVATE ${../src})
#add_executable(slab slab.cpp ${HEADER_FILES})
#target_include_directories(slab PRIVATE ${../src})

#add_executable(matrixTest matrix_test.cpp math_utilities.h math_utilities.cpp)
#add_executable(myTest main.cpp elements.cpp mesh.cpp stiffnessMatrix.cpp math_utilities.cpp elements.h mesh.h stiffnessMatrix.h math_utilities.h)
#add_executable(stiff stiffness_test.cpp elements.cpp mesh.cpp stiffnessMatrix.cpp math_utilities.cpp mesh.h elements.h stiffnessMatrix.h math_utilities.h)
#add_executable(plate plate.cpp elements.cpp mesh.cpp stiffnessMatrix.cpp math_utilities.cpp mesh.h elements.h stiffnessMatrix.h math_utilities.h)
add_executable(slab slab.cpp ../src/elements.cpp ../src/mesh.cpp ../src/stiffnessMatrix.cpp ../src/math_utilities.cpp ../include/mesh.h ../include/elements.h ../include/stiffnessMatrix.h ../include/math_utilities.h)
add_executable(slab2 slab2.cpp ../src/elements.cpp ../src/mesh.cpp ../src/stiffnessMatrix.cpp ../src/math_utilities.cpp ../include/mesh.h ../include/elements.h ../include/stiffnessMatrix.h ../include/math_utilities.h)
add_executable(slab3 slab3.cpp ../src/elements.cpp ../src/mesh.cpp ../src/stiffnessMatrix.cpp ../src/math_utilities.cpp ../include/mesh.h ../include/elements.h ../include/stiffnessMatrix.h ../include/math_utilities.h)
add_executable(plate plate.cpp ../src/elements.cpp ../src/mesh.cpp ../src/stiffnessMatrix.cpp ../src/math_utilities.cpp ../include/mesh.h ../include/elements.h ../include/stiffnessMatrix.h ../include/math_utilities.h)
#add_executable(stiff stiffness_test.cpp ../src/elements.cpp ../src/mesh.cpp ../src/stiffnessMatrix.cpp ../src/math_utilities.cpp ../include/mesh.h ../include/elements.h ../include/stiffnessMatrix.h ../include/math_utilities.h)
#add_executable(matrixTest matrix_test.cpp ../src/elements.cpp ../src/mesh.cpp ../src/stiffnessMatrix.cpp ../src/math_utilities.cpp ../include/mesh.h ../include/elements.h ../include/stiffnessMatrix.h ../include/math_utilities.h)
