# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.10

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:


#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:


# Remove some rules from gmake that .SUFFIXES does not remove.
SUFFIXES =

.SUFFIXES: .hpux_make_needs_suffix_list


# Suppress display of executed commands.
$(VERBOSE).SILENT:


# A target that is always out of date.
cmake_force:

.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /usr/bin/cmake

# The command to remove a file.
RM = /usr/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /mnt/d/grad_school/fea_solver/solver/test

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /mnt/d/grad_school/fea_solver/solver/test

# Include any dependencies generated for this target.
include CMakeFiles/matrixTest.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/matrixTest.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/matrixTest.dir/flags.make

CMakeFiles/matrixTest.dir/matrix_test.cpp.o: CMakeFiles/matrixTest.dir/flags.make
CMakeFiles/matrixTest.dir/matrix_test.cpp.o: matrix_test.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/grad_school/fea_solver/solver/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/matrixTest.dir/matrix_test.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrixTest.dir/matrix_test.cpp.o -c /mnt/d/grad_school/fea_solver/solver/test/matrix_test.cpp

CMakeFiles/matrixTest.dir/matrix_test.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrixTest.dir/matrix_test.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/grad_school/fea_solver/solver/test/matrix_test.cpp > CMakeFiles/matrixTest.dir/matrix_test.cpp.i

CMakeFiles/matrixTest.dir/matrix_test.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrixTest.dir/matrix_test.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/grad_school/fea_solver/solver/test/matrix_test.cpp -o CMakeFiles/matrixTest.dir/matrix_test.cpp.s

CMakeFiles/matrixTest.dir/matrix_test.cpp.o.requires:

.PHONY : CMakeFiles/matrixTest.dir/matrix_test.cpp.o.requires

CMakeFiles/matrixTest.dir/matrix_test.cpp.o.provides: CMakeFiles/matrixTest.dir/matrix_test.cpp.o.requires
	$(MAKE) -f CMakeFiles/matrixTest.dir/build.make CMakeFiles/matrixTest.dir/matrix_test.cpp.o.provides.build
.PHONY : CMakeFiles/matrixTest.dir/matrix_test.cpp.o.provides

CMakeFiles/matrixTest.dir/matrix_test.cpp.o.provides.build: CMakeFiles/matrixTest.dir/matrix_test.cpp.o


CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o: CMakeFiles/matrixTest.dir/flags.make
CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o: /mnt/d/grad_school/fea_solver/solver/src/elements.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/grad_school/fea_solver/solver/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o -c /mnt/d/grad_school/fea_solver/solver/src/elements.cpp

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/grad_school/fea_solver/solver/src/elements.cpp > CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.i

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/grad_school/fea_solver/solver/src/elements.cpp -o CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.s

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o.requires:

.PHONY : CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o.requires

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o.provides: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o.requires
	$(MAKE) -f CMakeFiles/matrixTest.dir/build.make CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o.provides.build
.PHONY : CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o.provides

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o.provides.build: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o


CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o: CMakeFiles/matrixTest.dir/flags.make
CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o: /mnt/d/grad_school/fea_solver/solver/src/mesh.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/grad_school/fea_solver/solver/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o -c /mnt/d/grad_school/fea_solver/solver/src/mesh.cpp

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/grad_school/fea_solver/solver/src/mesh.cpp > CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.i

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/grad_school/fea_solver/solver/src/mesh.cpp -o CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.s

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o.requires:

.PHONY : CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o.requires

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o.provides: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o.requires
	$(MAKE) -f CMakeFiles/matrixTest.dir/build.make CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o.provides.build
.PHONY : CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o.provides

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o.provides.build: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o


CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o: CMakeFiles/matrixTest.dir/flags.make
CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o: /mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/grad_school/fea_solver/solver/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o -c /mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp > CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.i

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp -o CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.s

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o.requires:

.PHONY : CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o.requires

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o.provides: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o.requires
	$(MAKE) -f CMakeFiles/matrixTest.dir/build.make CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o.provides.build
.PHONY : CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o.provides

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o.provides.build: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o


CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o: CMakeFiles/matrixTest.dir/flags.make
CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o: /mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/mnt/d/grad_school/fea_solver/solver/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o"
	/usr/bin/clang++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o -c /mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.i"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp > CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.i

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.s"
	/usr/bin/clang++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp -o CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.s

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o.requires:

.PHONY : CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o.requires

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o.provides: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o.requires
	$(MAKE) -f CMakeFiles/matrixTest.dir/build.make CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o.provides.build
.PHONY : CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o.provides

CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o.provides.build: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o


# Object files for target matrixTest
matrixTest_OBJECTS = \
"CMakeFiles/matrixTest.dir/matrix_test.cpp.o" \
"CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o" \
"CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o" \
"CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o" \
"CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o"

# External object files for target matrixTest
matrixTest_EXTERNAL_OBJECTS =

matrixTest: CMakeFiles/matrixTest.dir/matrix_test.cpp.o
matrixTest: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o
matrixTest: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o
matrixTest: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o
matrixTest: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o
matrixTest: CMakeFiles/matrixTest.dir/build.make
matrixTest: CMakeFiles/matrixTest.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/mnt/d/grad_school/fea_solver/solver/test/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Linking CXX executable matrixTest"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/matrixTest.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/matrixTest.dir/build: matrixTest

.PHONY : CMakeFiles/matrixTest.dir/build

CMakeFiles/matrixTest.dir/requires: CMakeFiles/matrixTest.dir/matrix_test.cpp.o.requires
CMakeFiles/matrixTest.dir/requires: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/elements.cpp.o.requires
CMakeFiles/matrixTest.dir/requires: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/mesh.cpp.o.requires
CMakeFiles/matrixTest.dir/requires: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/stiffnessMatrix.cpp.o.requires
CMakeFiles/matrixTest.dir/requires: CMakeFiles/matrixTest.dir/mnt/d/grad_school/fea_solver/solver/src/math_utilities.cpp.o.requires

.PHONY : CMakeFiles/matrixTest.dir/requires

CMakeFiles/matrixTest.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/matrixTest.dir/cmake_clean.cmake
.PHONY : CMakeFiles/matrixTest.dir/clean

CMakeFiles/matrixTest.dir/depend:
	cd /mnt/d/grad_school/fea_solver/solver/test && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /mnt/d/grad_school/fea_solver/solver/test /mnt/d/grad_school/fea_solver/solver/test /mnt/d/grad_school/fea_solver/solver/test /mnt/d/grad_school/fea_solver/solver/test /mnt/d/grad_school/fea_solver/solver/test/CMakeFiles/matrixTest.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/matrixTest.dir/depend

