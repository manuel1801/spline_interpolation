#
# This is based on orocos sample cmakelists, modified slightly
# The find_package macro for Orocos-RTT works best with
# cmake >= 2.8.3
#
cmake_minimum_required(VERSION 2.8.3)
 
#
# This creates a standard cmake project. You may extend this file with
# any cmake macro you see fit.
#
project(spline_interpolation)

if(CMAKE_COMPILER_IS_GNUCXX OR ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang"))
    set(CMAKE_CXX_FLAGS "-std=c++11 ${CMAKE_CXX_FLAGS}")
endif()

# find installed eigen library
find_package(Eigen3 REQUIRED)
include_directories(SYSTEM ${EIGEN3_INCLUDE_DIRS})

# current source and include:
include_directories(${PROJECT_SOURCE_DIR}/include)
 
# Set the CMAKE_PREFIX_PATH in case you're not using Orocos through ROS
# for helping these find commands find RTT.
# find_package(OROCOS-RTT REQUIRED ${RTT_HINTS})
# if (NOT OROCOS-RTT_FOUND)
#   message (FATAL_ERROR "\nCould not find Orocos. Please use the shell command\n 'source orocos_toolchain/env.sh' and then run cmake again.")
# endif()

# find_package(RST-RT REQUIRED)
# IF (RST-RT_FOUND)
# message(STATUS "RST-RT found. Version: ${RST-RT_VERSION}")
# include_directories(BEFORE SYSTEM ${RST-RT_INCLUDE_DIRS})
# add_definitions(${RST-RT_CFLAGS})
# link_directories(${RST-RT_LIBRARY_DIRS})
# ENDIF()



find_package(orocos_kdl REQUIRED)
message(STATUS "orocos_kdl_INCLUDE_DIRS: ${orocos_kdl_INCLUDE_DIRS}")
include_directories(${orocos_kdl_INCLUDE_DIRS})

# find_package(kdl_parser REQUIRED)
# message(STATUS "kdl_parser_INCLUDE_DIRS: ${kdl_parser_INCLUDE_DIRS}")
# include_directories(${kdl_parser_INCLUDE_DIRS})

#find_package(urdf REQUIRED)

# find_library(PORT PORT)



# Defines the orocos_* cmake macros. See that file for additional
# documentation.
# include(${OROCOS-RTT_USE_FILE_PATH}/UseOROCOS-RTT.cmake)
 
#
# Components, types and plugins.
#
# The CMake 'target' names are identical to the first argument of the
# macros below, except for orocos_typegen_headers, where the target is fully
# controlled by generated code of 'typegen'.
#

 
# Creates a component library libexample-<target>.so
# and installs in the directory lib/orocos/example/
#

# add_subdirectory(urdf-to-kdl)
# include_directories(${KDL_Parser_INCLUDE_DIRS})

add_executable(spline_interpolation src/spline_interpolation.cpp)
target_link_libraries(spline_interpolation ${orocos_kdl_LIBRARIES} ${Eigen3})
# target_link_libraries(trajectory ${RST-RT_LIBRARIES} ${RTT-CORE-EXTENSIONS_LIBRARIES} ${orocos_kdl_LIBRARIES} ${kdl_parser_LIBRARIES} urdf ${PORT})

