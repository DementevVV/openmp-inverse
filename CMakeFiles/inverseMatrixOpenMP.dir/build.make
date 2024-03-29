# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 2.8

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
CMAKE_SOURCE_DIR = /home/i/Documents/parallel/openmp-inverse

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/i/Documents/parallel/openmp-inverse

# Include any dependencies generated for this target.
include CMakeFiles/inverseMatrixOpenMP.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/inverseMatrixOpenMP.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/inverseMatrixOpenMP.dir/flags.make

CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o: CMakeFiles/inverseMatrixOpenMP.dir/flags.make
CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o: src/main.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/i/Documents/parallel/openmp-inverse/CMakeFiles $(CMAKE_PROGRESS_1)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o   -c /home/i/Documents/parallel/openmp-inverse/src/main.c

CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/i/Documents/parallel/openmp-inverse/src/main.c > CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.i

CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/i/Documents/parallel/openmp-inverse/src/main.c -o CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.s

CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o.requires:
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o.requires

CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o.provides: CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o.requires
	$(MAKE) -f CMakeFiles/inverseMatrixOpenMP.dir/build.make CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o.provides.build
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o.provides

CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o.provides.build: CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o

CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o: CMakeFiles/inverseMatrixOpenMP.dir/flags.make
CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o: src/matrix.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/i/Documents/parallel/openmp-inverse/CMakeFiles $(CMAKE_PROGRESS_2)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o   -c /home/i/Documents/parallel/openmp-inverse/src/matrix.c

CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/i/Documents/parallel/openmp-inverse/src/matrix.c > CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.i

CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/i/Documents/parallel/openmp-inverse/src/matrix.c -o CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.s

CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o.requires:
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o.requires

CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o.provides: CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o.requires
	$(MAKE) -f CMakeFiles/inverseMatrixOpenMP.dir/build.make CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o.provides.build
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o.provides

CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o.provides.build: CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o

CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o: CMakeFiles/inverseMatrixOpenMP.dir/flags.make
CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o: src/lu.c
	$(CMAKE_COMMAND) -E cmake_progress_report /home/i/Documents/parallel/openmp-inverse/CMakeFiles $(CMAKE_PROGRESS_3)
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Building C object CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -o CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o   -c /home/i/Documents/parallel/openmp-inverse/src/lu.c

CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -E /home/i/Documents/parallel/openmp-inverse/src/lu.c > CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.i

CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_FLAGS) -S /home/i/Documents/parallel/openmp-inverse/src/lu.c -o CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.s

CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o.requires:
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o.requires

CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o.provides: CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o.requires
	$(MAKE) -f CMakeFiles/inverseMatrixOpenMP.dir/build.make CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o.provides.build
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o.provides

CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o.provides.build: CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o

# Object files for target inverseMatrixOpenMP
inverseMatrixOpenMP_OBJECTS = \
"CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o" \
"CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o" \
"CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o"

# External object files for target inverseMatrixOpenMP
inverseMatrixOpenMP_EXTERNAL_OBJECTS =

inverseMatrixOpenMP: CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o
inverseMatrixOpenMP: CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o
inverseMatrixOpenMP: CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o
inverseMatrixOpenMP: CMakeFiles/inverseMatrixOpenMP.dir/build.make
inverseMatrixOpenMP: CMakeFiles/inverseMatrixOpenMP.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --red --bold "Linking C executable inverseMatrixOpenMP"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/inverseMatrixOpenMP.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/inverseMatrixOpenMP.dir/build: inverseMatrixOpenMP
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/build

CMakeFiles/inverseMatrixOpenMP.dir/requires: CMakeFiles/inverseMatrixOpenMP.dir/src/main.c.o.requires
CMakeFiles/inverseMatrixOpenMP.dir/requires: CMakeFiles/inverseMatrixOpenMP.dir/src/matrix.c.o.requires
CMakeFiles/inverseMatrixOpenMP.dir/requires: CMakeFiles/inverseMatrixOpenMP.dir/src/lu.c.o.requires
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/requires

CMakeFiles/inverseMatrixOpenMP.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/inverseMatrixOpenMP.dir/cmake_clean.cmake
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/clean

CMakeFiles/inverseMatrixOpenMP.dir/depend:
	cd /home/i/Documents/parallel/openmp-inverse && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/i/Documents/parallel/openmp-inverse /home/i/Documents/parallel/openmp-inverse /home/i/Documents/parallel/openmp-inverse /home/i/Documents/parallel/openmp-inverse /home/i/Documents/parallel/openmp-inverse/CMakeFiles/inverseMatrixOpenMP.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/inverseMatrixOpenMP.dir/depend

