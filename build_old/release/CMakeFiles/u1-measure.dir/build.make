# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.16

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
CMAKE_SOURCE_DIR = /home/simone/Documents/bonn_work/U1_project/su2

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/simone/Documents/bonn_work/U1_project/su2/build/release

# Include any dependencies generated for this target.
include CMakeFiles/u1-measure.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/u1-measure.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/u1-measure.dir/flags.make

CMakeFiles/u1-measure.dir/measure-u1.cc.o: CMakeFiles/u1-measure.dir/flags.make
CMakeFiles/u1-measure.dir/measure-u1.cc.o: ../../measure-u1.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/simone/Documents/bonn_work/U1_project/su2/build/release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/u1-measure.dir/measure-u1.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/u1-measure.dir/measure-u1.cc.o -c /home/simone/Documents/bonn_work/U1_project/su2/measure-u1.cc

CMakeFiles/u1-measure.dir/measure-u1.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/u1-measure.dir/measure-u1.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/simone/Documents/bonn_work/U1_project/su2/measure-u1.cc > CMakeFiles/u1-measure.dir/measure-u1.cc.i

CMakeFiles/u1-measure.dir/measure-u1.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/u1-measure.dir/measure-u1.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/simone/Documents/bonn_work/U1_project/su2/measure-u1.cc -o CMakeFiles/u1-measure.dir/measure-u1.cc.s

# Object files for target u1-measure
u1__measure_OBJECTS = \
"CMakeFiles/u1-measure.dir/measure-u1.cc.o"

# External object files for target u1-measure
u1__measure_EXTERNAL_OBJECTS =

u1-measure: CMakeFiles/u1-measure.dir/measure-u1.cc.o
u1-measure: CMakeFiles/u1-measure.dir/build.make
u1-measure: libsu2.a
u1-measure: libsu2.a
u1-measure: /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.71.0
u1-measure: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.71.0
u1-measure: /usr/lib/x86_64-linux-gnu/libboost_system.so.1.71.0
u1-measure: CMakeFiles/u1-measure.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/simone/Documents/bonn_work/U1_project/su2/build/release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable u1-measure"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/u1-measure.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/u1-measure.dir/build: u1-measure

.PHONY : CMakeFiles/u1-measure.dir/build

CMakeFiles/u1-measure.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/u1-measure.dir/cmake_clean.cmake
.PHONY : CMakeFiles/u1-measure.dir/clean

CMakeFiles/u1-measure.dir/depend:
	cd /home/simone/Documents/bonn_work/U1_project/su2/build/release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/simone/Documents/bonn_work/U1_project/su2 /home/simone/Documents/bonn_work/U1_project/su2 /home/simone/Documents/bonn_work/U1_project/su2/build/release /home/simone/Documents/bonn_work/U1_project/su2/build/release /home/simone/Documents/bonn_work/U1_project/su2/build/release/CMakeFiles/u1-measure.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/u1-measure.dir/depend

