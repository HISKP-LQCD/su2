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
CMAKE_BINARY_DIR = /home/simone/Documents/bonn_work/U1_project/su2/build/debug

# Include any dependencies generated for this target.
include CMakeFiles/scaling.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/scaling.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/scaling.dir/flags.make

CMakeFiles/scaling.dir/scaling.cc.o: CMakeFiles/scaling.dir/flags.make
CMakeFiles/scaling.dir/scaling.cc.o: ../../scaling.cc
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/simone/Documents/bonn_work/U1_project/su2/build/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/scaling.dir/scaling.cc.o"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/scaling.dir/scaling.cc.o -c /home/simone/Documents/bonn_work/U1_project/su2/scaling.cc

CMakeFiles/scaling.dir/scaling.cc.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/scaling.dir/scaling.cc.i"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/simone/Documents/bonn_work/U1_project/su2/scaling.cc > CMakeFiles/scaling.dir/scaling.cc.i

CMakeFiles/scaling.dir/scaling.cc.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/scaling.dir/scaling.cc.s"
	/usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/simone/Documents/bonn_work/U1_project/su2/scaling.cc -o CMakeFiles/scaling.dir/scaling.cc.s

# Object files for target scaling
scaling_OBJECTS = \
"CMakeFiles/scaling.dir/scaling.cc.o"

# External object files for target scaling
scaling_EXTERNAL_OBJECTS =

scaling: CMakeFiles/scaling.dir/scaling.cc.o
scaling: CMakeFiles/scaling.dir/build.make
scaling: libsu2.a
scaling: libsu2.a
scaling: /usr/lib/x86_64-linux-gnu/libboost_program_options.so.1.71.0
scaling: /usr/lib/x86_64-linux-gnu/libboost_filesystem.so.1.71.0
scaling: /usr/lib/x86_64-linux-gnu/libboost_system.so.1.71.0
scaling: CMakeFiles/scaling.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/simone/Documents/bonn_work/U1_project/su2/build/debug/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable scaling"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/scaling.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/scaling.dir/build: scaling

.PHONY : CMakeFiles/scaling.dir/build

CMakeFiles/scaling.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/scaling.dir/cmake_clean.cmake
.PHONY : CMakeFiles/scaling.dir/clean

CMakeFiles/scaling.dir/depend:
	cd /home/simone/Documents/bonn_work/U1_project/su2/build/debug && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/simone/Documents/bonn_work/U1_project/su2 /home/simone/Documents/bonn_work/U1_project/su2 /home/simone/Documents/bonn_work/U1_project/su2/build/debug /home/simone/Documents/bonn_work/U1_project/su2/build/debug /home/simone/Documents/bonn_work/U1_project/su2/build/debug/CMakeFiles/scaling.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/scaling.dir/depend

