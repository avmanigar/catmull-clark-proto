cmake_minimum_required(VERSION 2.8)

project(subsurf)

find_package(OpenGL REQUIRED)

if (APPLE)
  set(CMAKE_MACOSX_RPATH 1)
endif()

if (UNIX)
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --std=gnu++11")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
  set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wno-unused-variable")
  # recommended but not set by default
  # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
elseif(MSVC)
  # recommended but not set by default
  # set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -WX")
endif()

set (SUBSURF_LIBS ${OPENGL_gl_LIBRARY})

# GLFW
set(GLFW_INSTALL OFF)
set(GLFW_BUILD_DOCS OFF)
set(GLFW_BUILD_TESTS OFF)
set(GLFW_BUILD_EXAMPLES OFF)
add_subdirectory(glfw)
list(APPEND SUBSURF_LIBS glfw)
list(APPEND SUBSURF_INCLUDES glfw/include)

# GLEW - not needed on OS X
# we add glew source/header directly to the build, no glew library build.
if (NOT APPLE)
  add_definitions(-DGLEW_STATIC)
  list(APPEND SUBSURF_INCLUDES glew/include)
  list(APPEND SUBSURF_SRC glew/src/glew.c)
  SOURCE_GROUP(GLEW FILES glew/src/glew.c)
endif()

# vecmath include directory
include_directories(vecmath/include)
add_subdirectory(vecmath)
list (APPEND SUBSURF_LIBS vecmath)
list (APPEND SUBSURF_INCLUDES vecmath/include)
list (APPEND SUBSURF_SRC
  src/camera.cpp
  src/main.cpp
  src/starter0_util.cpp
  src/vertexrecorder.cpp
)
list (APPEND SUBSURF_HEADER
  src/camera.h
  src/starter0_util.h
  src/teapot.h
  src/tuple.h
  src/gl.h
  src/vertexrecorder.h
)

add_executable(subsurf ${SUBSURF_SRC} ${SUBSURF_HEADER})
target_include_directories(subsurf PUBLIC ${SUBSURF_INCLUDES})
target_link_libraries(subsurf ${SUBSURF_LIBS})
