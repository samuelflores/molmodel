# FindSimbody.cmake
# 
# Simbios National Center for Physics Based Simulation of Biological Structures
# Stanford University
# This cmake file created 2012 by Michael Sherman and is in the public 
# domain; Simbody itself is open source under the Apache 2.0 license.
#
# This is a CMake "find" module that will try to find the Simbody multibody 
# dynamics (physics) package installed somewhere on your computer. Simbody
# is part of the SimTK biosimulation toolkit. For more information, see
# https://simtk.org/home/simbody.
#
# To use this in a find_package() command from your own CMakeLists.txt file, 
# make sure this file is in a directory that is in the CMAKE_MODULE_PATH. 
# You can add a directory to that path with a line like this:
#     list(APPEND CMAKE_MODULE_PATH "myModuleDir")
#
# Then, to use Simbody include these lines:
#
#     find_package(Simbody REQUIRED)
#     include_directories(${Simbody_INCLUDE_DIR})
#     link_directories(${Simbody_LIB_DIR})
#     add_executable(myexe ${my_source_files} ${my_header_files})
#     target_link_libraries(myexe ${Simbody_LIBRARIES})
#
# If you don't want to make it REQUIRED, you can check Simbody_FOUND after
# find_package() returns.
# TODO: no version selection is implemented here yet; if you provide it
# to find_package() it will be ignored.
#
# This includes several libraries:
#     SimTKsimbody
#     SimTKmath
#     SimTKcommon
#     Windows only: liblapack libblas pthreadVC2[_x64]
# The above will be in Simbody_ROOT_DIR/lib.
#
# On Mac and Linux we don't provide our own lapack but expect it to be 
# available.
#     Mac/Linux only: lapack blas
#
# Once done this will define:
#
#   Simbody_FOUND - Whether search for Simbody libraries and headers succeeded.
#   Simbody_ROOT_DIR - the installation directory; all the pieces must be
#                      found together
#   Simbody_INCLUDE_DIR - location of Simbody.h
#   Simbody_LIB_DIR     - location of libSimTKsimbody.{a,so,dylib} or SimTKsimbody.lib
#   Simbody_BIN_DIR     - location of VisualizerGUI and .dll's on Windows
#   Simbody_LIBRARIES   - suitable for target_link_libraries(); includes
#                           both optimized and debug libraries if both are
#                           available
#   Simbody_STATIC_LIBRARIES - suitable for target_link_libraries(); includes
#                              both optimized and debug static libraries if
#                              both are available
#
# The following environment variables are used if set, in order of decreasing
# preference:
#   SIMBODY_HOME
#   SimTK_INSTALL_DIR
#
# Otherwise we look in platform-specific standard places, in this order:
#   <standardPlaces>/Simbody, /simbody
#   <standardPlaces>/SimTK, /simtk
#
# This module looks for certain CMake variables on input and behaves 
# accordingly if they are present:
#
#   SimTK_INSTALL_DIR
#       This is commonly set by other SimTK software and overrides the
#       environment variables if present. Note that this has the same name
#       as one of the environment variables but is distinct.

cmake_minimum_required(VERSION 3.10)

# Get values of relevant environment variables for convenient testing.
set(ENV_SIMBODY_HOME_VALUE $ENV{SIMBODY_HOME})
set(ENV_SimTK_INSTALL_DIR_VALUE $ENV{SimTK_INSTALL_DIR})

if(ENV_SIMBODY_HOME_VALUE)
    set(HINT_DIR ${ENV_SIMBODY_HOME_VALUE})
else()
    set(HINT_DIR ${ENV_SimTK_INSTALL_DIR_VALUE})
endif()

if (HINT_DIR)
    find_package(Simbody CONFIG REQUIRED HINTS "${HINT_DIR}")
else()
    find_package(Simbody CONFIG REQUIRED)
endif()
