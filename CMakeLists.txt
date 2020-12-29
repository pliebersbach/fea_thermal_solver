cmake_minimum_required(VERSION 3.0.0)
project(myTest VERSION 0.1.0)

set(CMAKE_BUILD_TYPE debug)
#set(CMAKE_CXX_COMPILER "/usr/bin/clang++")

ADD_SUBDIRECTORY(test)

INCLUDE_DIRECTORIES(include)

ADD_LIBRARY(mySrc 
src/elements.cpp
src/mesh.cpp
src/stiffnessMatrix.cpp
src/math_utilities.cpp)

ADD_EXECUTABLE(meshreader mesh_reader.cpp)

TARGET_LINK_LIBRARIES(meshreader mySrc)
