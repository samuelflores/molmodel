#---------------------------------------------------
# Molmodel
#
# Creates SimTK library, base name=SimTKmolmodel.
# Default libraries are shared & optimized. Variants
# are created for static (_static) and debug (_d) and
# provision is made for an optional "namespace" (ns)
# and version number (vn).
#
# Windows:
#   [ns_]SimTKmolmodel[_vn][_d].dll
#   [ns_]SimTKmolmodel[_vn][_d].lib
#   [ns_]SimTKmolmodel[_vn]_static[_d].lib
# Unix:
#   lib[ns_]SimTKmolmodel[_vn][_d].so
#   lib[ns_]SimTKmolmodel[_vn]_static[_d].a
#
# All libraries are installed in
#   %ProgramFiles%\SimTK\lib    (Windows)
#   /usr/local/SimTK/lib[64]    (Linux, Mac)
#
# Also creates an OpenMM plugin DLL that is used at
# runtime to determine whether OpenMM is available.
# That DLL is named
#   OpenMMPlugin[_d].dll
#   libOpenMMPlugin[_d].so
#   libOpenMMPlugin[_d].dylib
# And there is no static version.
#----------------------------------------------------

CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

PROJECT(Molmodel)

SET(MOLMODEL_MAJOR_VERSION 3)
SET(MOLMODEL_MINOR_VERSION 0)
SET(MOLMODEL_PATCH_VERSION 0)

SET(MOLMODEL_COPYRIGHT_YEARS "2006-12")

# underbar separated list of dotted authors, no spaces or commas
SET(MOLMODEL_AUTHORS "Christopher.Bruns_Michael.Sherman")
SET(CMAKE_CXX_STANDARD 14)
SET(CMAKE_CXX_EXTENSIONS OFF)

# Report the version number to the CMake UI. Don't include the
# build version if it is zero.
SET(PATCH_VERSION_STRING)
IF(MOLMODEL_PATCH_VERSION)
    SET(PATCH_VERSION_STRING ".${MOLMODEL_PATCH_VERSION}")
ENDIF()

SET(MOLMODEL_VERSION
    "${MOLMODEL_MAJOR_VERSION}.${MOLMODEL_MINOR_VERSION}${PATCH_VERSION_STRING}"
    CACHE STRING
    "This is the version that will be built (can't be changed here)."
    FORCE)

# This is the suffix if we're building and depending on versioned libraries.
SET(VN "_${MOLMODEL_VERSION}")

SET(BUILD_BINARY_DIR ${CMAKE_BINARY_DIR}
    CACHE PATH
    "The Molmodel build (not the install) puts all the libraries and executables together here (with /Release, etc. appended on some platforms).")

# Make everything go in the same binary directory. (These are CMake-defined
# variables.)
SET(EXECUTABLE_OUTPUT_PATH ${BUILD_BINARY_DIR})
SET(LIBRARY_OUTPUT_PATH ${BUILD_BINARY_DIR})

SET(BUILD_SHARED_LIBRARIES TRUE CACHE BOOL
    "Build dynamically linked (shared) version of the libraries")

# Static libraries, tests, and examples won't be built unless this
# is SET.
SET(BUILD_STATIC_LIBRARIES FALSE CACHE BOOL
    "Build '_static' versions of the libraries libraries?")
IF((NOT ${BUILD_SHARED_LIBRARIES}) AND (NOT ${BUILD_STATIC_LIBRARIES}))
    MESSAGE(FATAL_ERROR "Neither shared nor static build has been enabled. Aborting")
ENDIF()

# Use this to generate a private SET of libraries whose names
# won't conflict with installed versions.
SET(BUILD_USING_NAMESPACE "" CACHE STRING
	"All library names will be prefixed with 'xxx_' if this is SET to xxx.")

SET(BUILD_UNVERSIONED_LIBRARIES TRUE CACHE BOOL
    "Build library names, and assume dependency names, with no version numbers?")

SET(BUILD_VERSIONED_LIBRARIES FALSE CACHE BOOL
    "Build library names, and assume dependency names, with version numbers?")

SET(GEMMI_PATH "" CACHE PATH "Where is the Gemmi library installed?")

IF(WIN32)
    SET(zlib_INSTALL_DIR "" CACHE PATH
	    "Where to find zlib files")
ENDIF()

SET(NS)
IF(BUILD_USING_NAMESPACE)
    SET(NS "${BUILD_USING_NAMESPACE}_")
ENDIF()

SET(MOLMODEL_LIBRARY_NAME ${NS}SimTKmolmodel CACHE STRING
    "Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
    FORCE)


# Permit use of custom FindOpenMM and FindSimbody modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")
SET(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Control volume of build output" )

SET(OpenMM_ON OFF)
FIND_PACKAGE(OpenMM)
IF(OpenMM_FOUND)
    SET(OpenMM_ON ON)
ENDIF()
SET(MOLMODEL_USE_OpenMM ${OpenMM_ON} CACHE BOOL
    "Control whether Molmodel builds the OpenMM Plugin (requires that OpenMM is installed on the build machine)." )

# Caution: this variable is automatically created by the CMake
# ENABLE_TESTING() command, but we'll take it over here for
# our own purposes too.
SET(BUILD_TESTING ON CACHE BOOL
	"Control building of Molmodel test programs.")

SET(BUILD_EXAMPLES ON CACHE BOOL
	"Control building of Molmodel example programs.")
# Turning this off reduces the build time (and space) substantially,
# but you may miss the occasional odd bug. Also currently on Windows it
# is easier to debug the static tests than the DLL-liked ones.
SET(BUILD_TESTING_STATIC ON CACHE BOOL
    "If building static libraries, build static test and example programs too?")

SET(BUILD_TESTING_SHARED ON CACHE BOOL
    "If building test or example programs, include dynamically-linked ones?" )
# Create a platform name useful for some platform-specific stuff.

IF(WIN32)
    SET(NATIVE_COPY_CMD copy)
ELSEIF(APPLE)
    SET(NATIVE_COPY_CMD cp)
ELSE()
    SET(NATIVE_COPY_CMD cp)
ENDIF()

# In addition to the platform name we need to know the Application Binary
# Interface (ABI) we're building for. Currently that is either x86, meaning
# 32 bit Intel instruction SET, or x64 for 64 bit Intel instruction SET.
IF(${CMAKE_SIZEOF_VOID_P} EQUAL 8)
    SET(PLATFORM_ABI x64)
ELSE()
    SET(PLATFORM_ABI x86)
ENDIF()

IF(WIN32)
    ADD_DEFINITIONS(-DNOMINMAX -DWIN32_LEAN_AND_MEAN -D_USE_MATH_DEFINES -D_WIN32_WINNT=0x0501)

	IF(MSVC)
	    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
	ENDIF()
ENDIF()

SET(BUILD_PLATFORM "${CMAKE_HOST_SYSTEM_NAME}:${PLATFORM_ABI}" CACHE STRING
    "This is the platform and ABI we're building for. Not changeable here; use a different CMake generator instead."
    FORCE)

# Find Simbody using the local FindSimbody.cmake if present.
find_package(Simbody REQUIRED)
include_directories(${Simbody_INCLUDE_DIR})
link_directories(${Simbody_LIB_DIR})

IF(WIN32)
    IF (zlib_INSTALL_DIR STREQUAL "")
        MESSAGE(FATAL_ERROR "Path to zlib is not set")
	ENDIF()

	INCLUDE_DIRECTORIES(${zlib_INSTALL_DIR}/include)
	LINK_DIRECTORIES(${zlib_INSTALL_DIR}/lib)
ENDIF()

# If CMAKE_INSTALL_PREFIX is /usr/local, then the LIBDIR should necessarily be
# lib/. Sometimes (on Linux), LIBDIR is something like x86_64-linux-gnu. The
# linker will search /usr/lib/x86_64-linux-gnu (this path is in
# HOWEVER, it WILL search /usr/local/lib. So that Linux users needn't modify
# their LD_LIBRARY_PATH if installing to /usr/local, we force the LIBDIR to be
# lib/.
# Note: CMake 3.0 fixes this issue. When we move to CMake 3.0, we can
# remove this if-statement. See issue #151.
IF("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local" OR
   "${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local/")
    # Overwrite both of these variables; we use both of them.
    SET(CMAKE_INSTALL_LIBDIR "lib")
    SET(CMAKE_INSTALL_FULL_LIBDIR
        "${CMAKE_INSTALL_PREFIX}/${CMAKE_INSTALL_LIBDIR}")
ENDIF()


IF(NOT MSVC AND NOT XCODE AND NOT CMAKE_BUILD_TYPE)
    SET(CMAKE_BUILD_TYPE RelWithDebInfo CACHE STRING
        "Debug, RelWithDebInfo (recommended), or Release build"
        FORCE)
ENDIF()

## Choose the maximum level of x86 instruction SET that the compiler is
## allowed to use.
## Was using sse2 but changed to let the compilers choose. Most will
## probably use sse2 or later by default.
## On 64 bit MSVC 2013, the default is sse2 and the argument
## isn't recognized so don't specify it.
IF(CMAKE_CL_64)
    SET(default_build_inst_SET)
ELSE()
    SET(default_build_inst_SET)
ENDIF()

## This can be SET to a different value by the person running CMake.
SET(BUILD_INST_SET ""
    CACHE STRING
    "CPU instruction level compiler is permitted to use (default: let compiler decide).")
MARK_AS_ADVANCED(BUILD_INST_SET)

IF(BUILD_INST_SET)
    SET(inst_SET_to_use ${BUILD_INST_SET})
ELSE()
    SET(inst_SET_to_use ${default_build_inst_SET})
ENDIF()


## When building in any of the Release modes, tell VC++ cl compiler to use
## intrinsics (i.e. sqrt instruction rather than sqrt subroutine) with
## flag /Oi.
## Caution: can't use CMAKE_CXX_COMPILER_ID MATCHES MSVC here because
## "MSVC" is a predefined CMAKE variable and will get expanded to 1 or 0
IF(MSVC)
    IF(inst_SET_to_use)
        STRING(TOUPPER ${inst_SET_to_use} CL_INST_SET)
        SET(CL_INST_SET "/arch:${CL_INST_SET}")
    ELSE()
        SET(CL_INST_SET)
    ENDIF()

    SET(BUILD_LIMIT_PARALLEL_COMPILES "" CACHE STRING
        "Set a maximum number of simultaneous compilations.")
    MARK_AS_ADVANCED(BUILD_LIMIT_PARALLEL_COMPILES)
    SET(mxcpu ${BUILD_LIMIT_PARALLEL_COMPILES}) # abbreviation

    ## C++
    SET(BUILD_CXX_FLAGS_DEBUG
        "/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- ${CL_INST_SET}")
    SET(BUILD_CXX_FLAGS_RELEASE
        "/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- ${CL_INST_SET}")
    SET(BUILD_CXX_FLAGS_RELWITHDEBINFO
        "/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- ${CL_INST_SET}")
    SET(BUILD_CXX_FLAGS_MINSIZEREL
        "/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- ${CL_INST_SET}")

    SET(CMAKE_CXX_FLAGS_DEBUG "/MP${mxcpu} ${BUILD_CXX_FLAGS_DEBUG}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_RELEASE "/MP${mxcpu} ${BUILD_CXX_FLAGS_RELEASE}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MP${mxcpu} ${BUILD_CXX_FLAGS_RELWITHDEBINFO}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
    SET(CMAKE_CXX_FLAGS_MINSIZEREL "/MP${mxcpu} ${BUILD_CXX_FLAGS_MINSIZEREL}"
        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)
ENDIF()

# The source is organized into subdirectories, but we handle them all from
# this CMakeLists file rather than letting CMake visit them as SUBDIRS.
SET(MOLMODEL_SOURCE_SUBDIRS .)


# Collect up information about the version of the molmodel library we're building
# and make it available to the code so it can be built into the binaries.

# Get the subversion revision number if we can
# It's possible that WIN32 installs use svnversion through cygwin
# so we'll try for both svnversion.exe and svnversion. Note that
# this will result in warnings if all you have is Tortoise without
# Cygwin, and your "about" string will say "unknown" rather than
# providing the SVN version of the source.

# If there is an actual reason for the SVNVERSION vs. SVNVERSION_EXE
# nonsense that was here before, come talk to me.  --Chris Bruns
# 1) CMake already understands that Windows executables have ".exe" at the end.
# 2) Providing a PATHS argument to find_program is *NOT* hard coding a single path.
# Please read the CMake documentation for more details.
FIND_PROGRAM(SVNVERSION_EXECUTABLE svnversion PATHS "C:/cygwin/bin")

IF(SVNVERSION_EXECUTABLE)
    EXEC_PROGRAM(${SVNVERSION_EXECUTABLE}
              # Works better on Windows to SET working directory rather than
              # passing argument to svnversion
                 ${CMAKE_SOURCE_DIR} # cwd for run
                 OUTPUT_VARIABLE OUT)
    SET(MOLMODEL_SVN_REVISION "${OUT}" CACHE STRING "Molmodel svn revision number" FORCE)
ELSE(SVNVERSION_EXECUTABLE)
    MESSAGE(STATUS
            "Could not find 'svnversion' executable; 'about' will be wrong. (Cygwin provides one on Windows.)")
    SET(MOLMODEL_SVN_REVISION "unknown" CACHE STRING "Molmodel svn revision number")
ENDIF(SVNVERSION_EXECUTABLE)
MARK_AS_ADVANCED(MOLMODEL_SVN_REVISION)

# Remove colon from build version, for easier placement in directory names
STRING(REPLACE ":" "_" MOLMODEL_SVN_REVISION ${MOLMODEL_SVN_REVISION})

ADD_DEFINITIONS(-DSimTK_MOLMODEL_LIBRARY_NAME=${MOLMODEL_LIBRARY_NAME}
                -DSimTK_MOLMODEL_MAJOR_VERSION=${MOLMODEL_MAJOR_VERSION}
                -DSimTK_MOLMODEL_MINOR_VERSION=${MOLMODEL_MINOR_VERSION}
                -DSimTK_MOLMODEL_PATCH_VERSION=${MOLMODEL_PATCH_VERSION})

# CMake quotes automatically when building Visual Studio projects but we need
# to add them ourselves for Linux or Cygwin. Two cases to avoid duplicate quotes
# in Visual Studio which end up in the binary.

IF(${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(NEED_QUOTES FALSE)
ELSE(${CMAKE_GENERATOR} MATCHES "Visual Studio")
   SET(NEED_QUOTES TRUE)
ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio")

##TODO: doesn't work without quotes in nightly build
SET(NEED_QUOTES TRUE)

IF(NEED_QUOTES)
   ADD_DEFINITIONS(-DSimTK_MOLMODEL_SVN_REVISION="${MOLMODEL_SVN_REVISION}"
                   -DSimTK_MOLMODEL_COPYRIGHT_YEARS="${MOLMODEL_COPYRIGHT_YEARS}"
                   -DSimTK_MOLMODEL_AUTHORS="${MOLMODEL_AUTHORS}")
ELSE(NEED_QUOTES)
   ADD_DEFINITIONS(-DSimTK_MOLMODEL_SVN_REVISION=${MOLMODEL_SVN_REVISION}
                   -DSimTK_MOLMODEL_COPYRIGHT_YEARS=${MOLMODEL_COPYRIGHT_YEARS}
                   -DSimTK_MOLMODEL_AUTHORS=${MOLMODEL_AUTHORS})
ENDIF(NEED_QUOTES)

# -DSimTK_MOLMODEL_TYPE has to be defined in the target subdirectories.
# -Dmolmodel_EXPORTS defined automatically when Windows DLL build is being done.

SET(SHARED_TARGET ${MOLMODEL_LIBRARY_NAME})
SET(STATIC_TARGET ${MOLMODEL_LIBRARY_NAME}_static)
SET(SHARED_TARGET_VN ${MOLMODEL_LIBRARY_NAME}${VN})
SET(STATIC_TARGET_VN ${MOLMODEL_LIBRARY_NAME}${VN}_static)

# Ensure that debug libraries have "_d" appended to their names.
# CMake gets this right on Windows automatically with this definition.
IF(${CMAKE_GENERATOR} MATCHES "Visual Studio")
    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)
ENDIF(${CMAKE_GENERATOR} MATCHES "Visual Studio")

# But on Unix or Cygwin we have to add the suffix manually
IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    SET(SHARED_TARGET ${SHARED_TARGET}_d)
    SET(STATIC_TARGET ${STATIC_TARGET}_d)
    SET(SHARED_TARGET_VN ${SHARED_TARGET_VN}_d)
    SET(STATIC_TARGET_VN ${STATIC_TARGET_VN}_d)
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

## Test against the unversioned libraries if they are being built;
## otherwise against the versioned libraries.
IF(BUILD_UNVERSIONED_LIBRARIES)
    SET(TEST_SHARED_TARGET ${SHARED_TARGET} ${Simbody_LIBRARIES})
    SET(TEST_STATIC_TARGET ${STATIC_TARGET} ${Simbody_STATIC_LIBRARIES})
ELSE(BUILD_UNVERSIONED_LIBRARIES)
    SET(TEST_SHARED_TARGET ${SHARED_TARGET_VN} ${Simbody_LIBRARIES}_${Simbody_VERSION})
    SET(TEST_STATIC_TARGET ${STATIC_TARGET_VN} ${Simbody_STATIC_LIBRARIES}_${Simbody_VERSION})
ENDIF(BUILD_UNVERSIONED_LIBRARIES)

# These are all the places to search for header files which are
# to be part of the API.
SET(API_INCLUDE_DIRS) # start empty
FOREACH(subdir ${MOLMODEL_SOURCE_SUBDIRS})
    # append
    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}
                         ${subdir}/include
                         ${subdir}/include/molmodel
                         ${subdir}/include/molmodel/internal)
ENDFOREACH(subdir)

# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,
# and absolute pathnames.
SET(API_REL_INCLUDE_FILES)   # start these out empty
SET(API_ABS_INCLUDE_FILES)

FOREACH(dir ${API_INCLUDE_DIRS})
    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames
    SET(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})

    FOREACH(pathname ${fullpaths})
        GET_FILENAME_COMPONENT(filename ${pathname} NAME)
        SET(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})
    ENDFOREACH(pathname)
ENDFOREACH(dir)

# collect up source files
SET(SOURCE_FILES) # empty
SET(SOURCE_INCLUDE_FILES)

FOREACH(subdir ${MOLMODEL_SOURCE_SUBDIRS})
    FILE(GLOB src_files ${subdir}/src/*.cpp ${subdir}/src/*.c
                        ${subdir}/src/*/*.cpp ${subdir}/src/*/*.c)
    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h)


    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append
    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})

    ## Make sure we find these locally before looking in
    ## SimTK/include if Molmodel was previously installed there.
    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/${subdir}/include)
ENDFOREACH(subdir)

INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/src)
INCLUDE_DIRECTORIES(BEFORE ${GEMMI_PATH})

#
# Allow automated build and dashboard.
#

INCLUDE (Dart)
IF (BUILD_TESTING)
    #IF (UNIX AND NOT CYGWIN AND NOT APPLE)
    #  IF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
    #    ADD_DEFINITIONS(-fprofile-arcs -ftest-coverage)
    #    LINK_LIBRARIES(gcov)
    #  ENDIF (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE MATCHES Debug)
    #ENDIF (UNIX AND NOT CYGWIN AND NOT APPLE)
    #
    # Testing
    #
    ENABLE_TESTING()
ENDIF (BUILD_TESTING)

INCLUDE(ApiDoxygen.cmake)
# libraries are installed from their subdirectories; headers here
# install headers
FILE(GLOB CORE_HEADERS     include/*.h                  */include/*.h)
FILE(GLOB TOP_HEADERS      include/molmodel/*.h          */include/molmodel/*.h)
FILE(GLOB INTERNAL_HEADERS include/molmodel/internal/*.h */include/molmodel/internal/*.h)
INSTALL_FILES(/include/                  FILES ${CORE_HEADERS})
INSTALL_FILES(/include/molmodel/         FILES ${TOP_HEADERS})
INSTALL_FILES(/include/molmodel/internal FILES ${INTERNAL_HEADERS})

# Install documents.
FILE(GLOB TOPLEVEL_DOCS doc/*.pdf doc/*.txt)
INSTALL(FILES ${TOPLEVEL_DOCS} DESTINATION doc)

# all the various variables have been SET above.

IF(BUILD_STATIC_LIBRARIES)
    ADD_SUBDIRECTORY(staticTarget)
ENDIF()

IF(BUILD_SHARED_LIBRARIES)
    ADD_SUBDIRECTORY(sharedTarget)
ENDIF(BUILD_SHARED_LIBRARIES)

IF(MOLMODEL_USE_OpenMM)
    ADD_SUBDIRECTORY(OpenMMPlugin)
ENDIF (MOLMODEL_USE_OpenMM)

IF(BUILD_EXAMPLES)
    ADD_SUBDIRECTORY(examples)
ENDIF()

IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(tests)
ENDIF()
