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
CMAKE_SOURCE_DIR = /home/melvin/cernbox/development/projects/BEMMM

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/melvin/cernbox/development/projects/BEMMM/build

# Include any dependencies generated for this target.
include examples/CMakeFiles/example_uq_curved_domain.out.dir/depend.make

# Include the progress variables for this target.
include examples/CMakeFiles/example_uq_curved_domain.out.dir/progress.make

# Include the compile flags for this target's objects.
include examples/CMakeFiles/example_uq_curved_domain.out.dir/flags.make

examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o: examples/CMakeFiles/example_uq_curved_domain.out.dir/flags.make
examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o: ../examples/example_uq_curved_domain.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/melvin/cernbox/development/projects/BEMMM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o"
	cd /home/melvin/cernbox/development/projects/BEMMM/build/examples && /usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o -c /home/melvin/cernbox/development/projects/BEMMM/examples/example_uq_curved_domain.cpp

examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.i"
	cd /home/melvin/cernbox/development/projects/BEMMM/build/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/melvin/cernbox/development/projects/BEMMM/examples/example_uq_curved_domain.cpp > CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.i

examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.s"
	cd /home/melvin/cernbox/development/projects/BEMMM/build/examples && /usr/bin/c++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/melvin/cernbox/development/projects/BEMMM/examples/example_uq_curved_domain.cpp -o CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.s

examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o.requires:

.PHONY : examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o.requires

examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o.provides: examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o.requires
	$(MAKE) -f examples/CMakeFiles/example_uq_curved_domain.out.dir/build.make examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o.provides.build
.PHONY : examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o.provides

examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o.provides.build: examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o


# Object files for target example_uq_curved_domain.out
example_uq_curved_domain_out_OBJECTS = \
"CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o"

# External object files for target example_uq_curved_domain.out
example_uq_curved_domain_out_EXTERNAL_OBJECTS =

examples/example_uq_curved_domain.out: examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o
examples/example_uq_curved_domain.out: examples/CMakeFiles/example_uq_curved_domain.out.dir/build.make
examples/example_uq_curved_domain.out: examples/CMakeFiles/example_uq_curved_domain.out.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/melvin/cernbox/development/projects/BEMMM/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable example_uq_curved_domain.out"
	cd /home/melvin/cernbox/development/projects/BEMMM/build/examples && $(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/example_uq_curved_domain.out.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
examples/CMakeFiles/example_uq_curved_domain.out.dir/build: examples/example_uq_curved_domain.out

.PHONY : examples/CMakeFiles/example_uq_curved_domain.out.dir/build

examples/CMakeFiles/example_uq_curved_domain.out.dir/requires: examples/CMakeFiles/example_uq_curved_domain.out.dir/example_uq_curved_domain.cpp.o.requires

.PHONY : examples/CMakeFiles/example_uq_curved_domain.out.dir/requires

examples/CMakeFiles/example_uq_curved_domain.out.dir/clean:
	cd /home/melvin/cernbox/development/projects/BEMMM/build/examples && $(CMAKE_COMMAND) -P CMakeFiles/example_uq_curved_domain.out.dir/cmake_clean.cmake
.PHONY : examples/CMakeFiles/example_uq_curved_domain.out.dir/clean

examples/CMakeFiles/example_uq_curved_domain.out.dir/depend:
	cd /home/melvin/cernbox/development/projects/BEMMM/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/melvin/cernbox/development/projects/BEMMM /home/melvin/cernbox/development/projects/BEMMM/examples /home/melvin/cernbox/development/projects/BEMMM/build /home/melvin/cernbox/development/projects/BEMMM/build/examples /home/melvin/cernbox/development/projects/BEMMM/build/examples/CMakeFiles/example_uq_curved_domain.out.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : examples/CMakeFiles/example_uq_curved_domain.out.dir/depend

