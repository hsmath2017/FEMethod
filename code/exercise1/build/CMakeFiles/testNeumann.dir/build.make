# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.22

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
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
RM = /usr/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/shuang/FEMethod/code/exercise1

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/shuang/FEMethod/code/exercise1/build

# Include any dependencies generated for this target.
include CMakeFiles/testNeumann.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/testNeumann.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/testNeumann.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/testNeumann.dir/flags.make

CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o: CMakeFiles/testNeumann.dir/flags.make
CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o: ../src/Neumanntest.cpp
CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o: CMakeFiles/testNeumann.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/shuang/FEMethod/code/exercise1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -MD -MT CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o -MF CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o.d -o CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o -c /home/shuang/FEMethod/code/exercise1/src/Neumanntest.cpp

CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.i"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/shuang/FEMethod/code/exercise1/src/Neumanntest.cpp > CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.i

CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.s"
	/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/shuang/FEMethod/code/exercise1/src/Neumanntest.cpp -o CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.s

# Object files for target testNeumann
testNeumann_OBJECTS = \
"CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o"

# External object files for target testNeumann
testNeumann_EXTERNAL_OBJECTS =

../bin/testNeumann: CMakeFiles/testNeumann.dir/src/Neumanntest.cpp.o
../bin/testNeumann: CMakeFiles/testNeumann.dir/build.make
../bin/testNeumann: CMakeFiles/testNeumann.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/shuang/FEMethod/code/exercise1/build/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CXX executable ../bin/testNeumann"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/testNeumann.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/testNeumann.dir/build: ../bin/testNeumann
.PHONY : CMakeFiles/testNeumann.dir/build

CMakeFiles/testNeumann.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/testNeumann.dir/cmake_clean.cmake
.PHONY : CMakeFiles/testNeumann.dir/clean

CMakeFiles/testNeumann.dir/depend:
	cd /home/shuang/FEMethod/code/exercise1/build && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/shuang/FEMethod/code/exercise1 /home/shuang/FEMethod/code/exercise1 /home/shuang/FEMethod/code/exercise1/build /home/shuang/FEMethod/code/exercise1/build /home/shuang/FEMethod/code/exercise1/build/CMakeFiles/testNeumann.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/testNeumann.dir/depend

