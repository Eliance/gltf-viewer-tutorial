# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.7

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
CMAKE_SOURCE_DIR = /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer

# Include any dependencies generated for this target.
include CMakeFiles/gltf-viewer.dir/depend.make

# Include the progress variables for this target.
include CMakeFiles/gltf-viewer.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/gltf-viewer.dir/flags.make

bin/shaders/gltf-viewer/forward.vs.glsl: apps/gltf-viewer/shaders/forward.vs.glsl
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Generating bin/shaders/gltf-viewer/forward.vs.glsl"
	/usr/bin/cmake -E copy /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/shaders/forward.vs.glsl /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/bin/shaders/gltf-viewer/forward.vs.glsl

bin/shaders/gltf-viewer/magenta.fs.glsl: apps/gltf-viewer/shaders/magenta.fs.glsl
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Generating bin/shaders/gltf-viewer/magenta.fs.glsl"
	/usr/bin/cmake -E copy /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/shaders/magenta.fs.glsl /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/bin/shaders/gltf-viewer/magenta.fs.glsl

bin/shaders/gltf-viewer/normals.fs.glsl: apps/gltf-viewer/shaders/normals.fs.glsl
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --blue --bold --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Generating bin/shaders/gltf-viewer/normals.fs.glsl"
	/usr/bin/cmake -E copy /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/shaders/normals.fs.glsl /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/bin/shaders/gltf-viewer/normals.fs.glsl

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o: apps/gltf-viewer/ViewerApplication.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_4) "Building CXX object CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/ViewerApplication.cpp

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/ViewerApplication.cpp > CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.i

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/ViewerApplication.cpp -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.s

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o.requires

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o.provides: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o.provides

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o


CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o: apps/gltf-viewer/main.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_5) "Building CXX object CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/main.cpp

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/main.cpp > CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.i

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/main.cpp -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.s

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o.requires

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o.provides: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o.provides

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o


CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o: apps/gltf-viewer/tiny_gltf_impl.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_6) "Building CXX object CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/tiny_gltf_impl.cpp

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/tiny_gltf_impl.cpp > CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.i

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/tiny_gltf_impl.cpp -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.s

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o.requires

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o.provides: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o.provides

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o


CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o: apps/gltf-viewer/utils/cameras.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_7) "Building CXX object CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/cameras.cpp

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/cameras.cpp > CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.i

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/cameras.cpp -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.s

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o.requires

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o.provides: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o.provides

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o


CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o: apps/gltf-viewer/utils/gl_debug_output.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_8) "Building CXX object CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/gl_debug_output.cpp

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/gl_debug_output.cpp > CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.i

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/gl_debug_output.cpp -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.s

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o.requires

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o.provides: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o.provides

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o


CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o: apps/gltf-viewer/utils/gltf.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_9) "Building CXX object CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/gltf.cpp

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/gltf.cpp > CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.i

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/apps/gltf-viewer/utils/gltf.cpp -o CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.s

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o.requires

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o.provides: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o.provides

CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o


CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o: third-party/imgui-1.74/imgui.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_10) "Building CXX object CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui.cpp

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui.cpp > CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.i

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui.cpp -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.s

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o.requires

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o.provides: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o.provides

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o


CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o: third-party/imgui-1.74/imgui_demo.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_11) "Building CXX object CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_demo.cpp

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_demo.cpp > CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.i

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_demo.cpp -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.s

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o.requires

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o.provides: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o.provides

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o


CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o: third-party/imgui-1.74/imgui_draw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_12) "Building CXX object CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_draw.cpp

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_draw.cpp > CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.i

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_draw.cpp -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.s

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o.requires

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o.provides: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o.provides

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o


CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o: third-party/imgui-1.74/imgui_widgets.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_13) "Building CXX object CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_widgets.cpp

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_widgets.cpp > CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.i

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/imgui_widgets.cpp -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.s

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o.requires

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o.provides: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o.provides

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o


CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o: third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_14) "Building CXX object CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp > CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.i

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.s

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o.requires

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o.provides: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o.provides

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o


CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o: third-party/imgui-1.74/examples/imgui_impl_glfw.cpp
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_15) "Building CXX object CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o"
	/usr/bin/c++   $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CXX source to CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.i"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp > CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.i

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CXX source to assembly CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.s"
	/usr/bin/c++  $(CXX_DEFINES) $(CXX_INCLUDES) $(CXX_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp -o CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.s

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o.requires

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o.provides: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o.provides

CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o.provides.build: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o


CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o: CMakeFiles/gltf-viewer.dir/flags.make
CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o: third-party/glad/src/glad.c
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_16) "Building C object CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -o CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o   -c /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/glad/src/glad.c

CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing C source to CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.i"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -E /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/glad/src/glad.c > CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.i

CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling C source to assembly CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.s"
	/usr/bin/cc  $(C_DEFINES) $(C_INCLUDES) $(C_FLAGS) -S /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/third-party/glad/src/glad.c -o CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.s

CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o.requires:

.PHONY : CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o.requires

CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o.provides: CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o.requires
	$(MAKE) -f CMakeFiles/gltf-viewer.dir/build.make CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o.provides.build
.PHONY : CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o.provides

CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o.provides.build: CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o


# Object files for target gltf-viewer
gltf__viewer_OBJECTS = \
"CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o" \
"CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o" \
"CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o" \
"CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o" \
"CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o" \
"CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o" \
"CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o" \
"CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o" \
"CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o" \
"CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o" \
"CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o" \
"CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o" \
"CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o"

# External object files for target gltf-viewer
gltf__viewer_EXTERNAL_OBJECTS =

bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/build.make
bin/gltf-viewer: /usr/lib/x86_64-linux-gnu/libGLU.so
bin/gltf-viewer: /usr/lib/x86_64-linux-gnu/libGL.so
bin/gltf-viewer: third-party/glfw-3.3.1/src/libglfw3.a
bin/gltf-viewer: /usr/lib/x86_64-linux-gnu/librt.so
bin/gltf-viewer: /usr/lib/x86_64-linux-gnu/libm.so
bin/gltf-viewer: /usr/lib/x86_64-linux-gnu/libX11.so
bin/gltf-viewer: CMakeFiles/gltf-viewer.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles --progress-num=$(CMAKE_PROGRESS_17) "Linking CXX executable bin/gltf-viewer"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/gltf-viewer.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/gltf-viewer.dir/build: bin/gltf-viewer

.PHONY : CMakeFiles/gltf-viewer.dir/build

CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/ViewerApplication.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/main.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/tiny_gltf_impl.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/cameras.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gl_debug_output.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/apps/gltf-viewer/utils/gltf.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_demo.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_draw.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/imgui_widgets.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_opengl3.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/third-party/imgui-1.74/examples/imgui_impl_glfw.cpp.o.requires
CMakeFiles/gltf-viewer.dir/requires: CMakeFiles/gltf-viewer.dir/third-party/glad/src/glad.c.o.requires

.PHONY : CMakeFiles/gltf-viewer.dir/requires

CMakeFiles/gltf-viewer.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/gltf-viewer.dir/cmake_clean.cmake
.PHONY : CMakeFiles/gltf-viewer.dir/clean

CMakeFiles/gltf-viewer.dir/depend: bin/shaders/gltf-viewer/forward.vs.glsl
CMakeFiles/gltf-viewer.dir/depend: bin/shaders/gltf-viewer/magenta.fs.glsl
CMakeFiles/gltf-viewer.dir/depend: bin/shaders/gltf-viewer/normals.fs.glsl
	cd /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer /home/6im2/ejarcet/OpenGL/glTF-VIEWER-ROOT/gltf-viewer/CMakeFiles/gltf-viewer.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/gltf-viewer.dir/depend

