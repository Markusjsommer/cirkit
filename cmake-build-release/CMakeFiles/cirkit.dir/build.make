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


# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

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
CMAKE_COMMAND = /home/markus/cmake/cmake-3.16.6-Linux-x86_64/bin/cmake

# The command to remove a file.
RM = /home/markus/cmake/cmake-3.16.6-Linux-x86_64/bin/cmake -E remove -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /tmp/tmp.fbBXt0ZQ6K

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /tmp/tmp.fbBXt0ZQ6K/cmake-build-release

# Include any dependencies generated for this target.
include CMakeFiles/cirkit.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/cirkit.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/cirkit.dir/flags.make

CMakeFiles/cirkit.dir/main.cpp.o: CMakeFiles/cirkit.dir/flags.make
CMakeFiles/cirkit.dir/main.cpp.o: ../main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/tmp.fbBXt0ZQ6K/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CXX object CMakeFiles/cirkit.dir/main.cpp.o"
	/opt/rh/devtoolset-8/root/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cirkit.dir/main.cpp.o -c /tmp/tmp.fbBXt0ZQ6K/main.cpp

CMakeFiles/cirkit.dir/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cirkit.dir/main.cpp.i"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /tmp/tmp.fbBXt0ZQ6K/main.cpp > CMakeFiles/cirkit.dir/main.cpp.i

CMakeFiles/cirkit.dir/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cirkit.dir/main.cpp.s"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /tmp/tmp.fbBXt0ZQ6K/main.cpp -o CMakeFiles/cirkit.dir/main.cpp.s

CMakeFiles/cirkit.dir/include/FastaReader.cpp.o: CMakeFiles/cirkit.dir/flags.make
CMakeFiles/cirkit.dir/include/FastaReader.cpp.o: ../include/FastaReader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/tmp.fbBXt0ZQ6K/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Building CXX object CMakeFiles/cirkit.dir/include/FastaReader.cpp.o"
	/opt/rh/devtoolset-8/root/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cirkit.dir/include/FastaReader.cpp.o -c /tmp/tmp.fbBXt0ZQ6K/include/FastaReader.cpp

CMakeFiles/cirkit.dir/include/FastaReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cirkit.dir/include/FastaReader.cpp.i"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /tmp/tmp.fbBXt0ZQ6K/include/FastaReader.cpp > CMakeFiles/cirkit.dir/include/FastaReader.cpp.i

CMakeFiles/cirkit.dir/include/FastaReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cirkit.dir/include/FastaReader.cpp.s"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /tmp/tmp.fbBXt0ZQ6K/include/FastaReader.cpp -o CMakeFiles/cirkit.dir/include/FastaReader.cpp.s

CMakeFiles/cirkit.dir/include/build_index.cpp.o: CMakeFiles/cirkit.dir/flags.make
CMakeFiles/cirkit.dir/include/build_index.cpp.o: ../include/build_index.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/tmp.fbBXt0ZQ6K/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Building CXX object CMakeFiles/cirkit.dir/include/build_index.cpp.o"
	/opt/rh/devtoolset-8/root/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cirkit.dir/include/build_index.cpp.o -c /tmp/tmp.fbBXt0ZQ6K/include/build_index.cpp

CMakeFiles/cirkit.dir/include/build_index.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cirkit.dir/include/build_index.cpp.i"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /tmp/tmp.fbBXt0ZQ6K/include/build_index.cpp > CMakeFiles/cirkit.dir/include/build_index.cpp.i

CMakeFiles/cirkit.dir/include/build_index.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cirkit.dir/include/build_index.cpp.s"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /tmp/tmp.fbBXt0ZQ6K/include/build_index.cpp -o CMakeFiles/cirkit.dir/include/build_index.cpp.s

CMakeFiles/cirkit.dir/include/FileReader.cpp.o: CMakeFiles/cirkit.dir/flags.make
CMakeFiles/cirkit.dir/include/FileReader.cpp.o: ../include/FileReader.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/tmp.fbBXt0ZQ6K/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/cirkit.dir/include/FileReader.cpp.o"
	/opt/rh/devtoolset-8/root/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cirkit.dir/include/FileReader.cpp.o -c /tmp/tmp.fbBXt0ZQ6K/include/FileReader.cpp

CMakeFiles/cirkit.dir/include/FileReader.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cirkit.dir/include/FileReader.cpp.i"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /tmp/tmp.fbBXt0ZQ6K/include/FileReader.cpp > CMakeFiles/cirkit.dir/include/FileReader.cpp.i

CMakeFiles/cirkit.dir/include/FileReader.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cirkit.dir/include/FileReader.cpp.s"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /tmp/tmp.fbBXt0ZQ6K/include/FileReader.cpp -o CMakeFiles/cirkit.dir/include/FileReader.cpp.s

CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.o: CMakeFiles/cirkit.dir/flags.make
CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.o: ../include/kmer_flathash.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/tmp.fbBXt0ZQ6K/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.o"
	/opt/rh/devtoolset-8/root/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.o -c /tmp/tmp.fbBXt0ZQ6K/include/kmer_flathash.cpp

CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.i"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /tmp/tmp.fbBXt0ZQ6K/include/kmer_flathash.cpp > CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.i

CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.s"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /tmp/tmp.fbBXt0ZQ6K/include/kmer_flathash.cpp -o CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.s

CMakeFiles/cirkit.dir/include/BitsetManager.cpp.o: CMakeFiles/cirkit.dir/flags.make
CMakeFiles/cirkit.dir/include/BitsetManager.cpp.o: ../include/BitsetManager.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/tmp/tmp.fbBXt0ZQ6K/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/cirkit.dir/include/BitsetManager.cpp.o"
	/opt/rh/devtoolset-8/root/usr/bin/g++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/cirkit.dir/include/BitsetManager.cpp.o -c /tmp/tmp.fbBXt0ZQ6K/include/BitsetManager.cpp

CMakeFiles/cirkit.dir/include/BitsetManager.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/cirkit.dir/include/BitsetManager.cpp.i"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /tmp/tmp.fbBXt0ZQ6K/include/BitsetManager.cpp > CMakeFiles/cirkit.dir/include/BitsetManager.cpp.i

CMakeFiles/cirkit.dir/include/BitsetManager.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/cirkit.dir/include/BitsetManager.cpp.s"
	/opt/rh/devtoolset-8/root/usr/bin/g++ $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /tmp/tmp.fbBXt0ZQ6K/include/BitsetManager.cpp -o CMakeFiles/cirkit.dir/include/BitsetManager.cpp.s

# Object files for target cirkit
cirkit_OBJECTS = \
"CMakeFiles/cirkit.dir/main.cpp.o" \
"CMakeFiles/cirkit.dir/include/FastaReader.cpp.o" \
"CMakeFiles/cirkit.dir/include/build_index.cpp.o" \
"CMakeFiles/cirkit.dir/include/FileReader.cpp.o" \
"CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.o" \
"CMakeFiles/cirkit.dir/include/BitsetManager.cpp.o"

# External object files for target cirkit
cirkit_EXTERNAL_OBJECTS =

cirkit: CMakeFiles/cirkit.dir/main.cpp.o
cirkit: CMakeFiles/cirkit.dir/include/FastaReader.cpp.o
cirkit: CMakeFiles/cirkit.dir/include/build_index.cpp.o
cirkit: CMakeFiles/cirkit.dir/include/FileReader.cpp.o
cirkit: CMakeFiles/cirkit.dir/include/kmer_flathash.cpp.o
cirkit: CMakeFiles/cirkit.dir/include/BitsetManager.cpp.o
cirkit: CMakeFiles/cirkit.dir/build.make
cirkit: /usr/lib64/libz.so
cirkit: CMakeFiles/cirkit.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/tmp/tmp.fbBXt0ZQ6K/cmake-build-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Linking CXX executable cirkit"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/cirkit.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/cirkit.dir/build: cirkit

.PHONY : CMakeFiles/cirkit.dir/build

CMakeFiles/cirkit.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/cirkit.dir/cmake_clean.cmake
.PHONY : CMakeFiles/cirkit.dir/clean

CMakeFiles/cirkit.dir/depend:
	cd /tmp/tmp.fbBXt0ZQ6K/cmake-build-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /tmp/tmp.fbBXt0ZQ6K /tmp/tmp.fbBXt0ZQ6K /tmp/tmp.fbBXt0ZQ6K/cmake-build-release /tmp/tmp.fbBXt0ZQ6K/cmake-build-release /tmp/tmp.fbBXt0ZQ6K/cmake-build-release/CMakeFiles/cirkit.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/cirkit.dir/depend

