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
#   [ns_]SimTKmolmodel[_d].dll
#   [ns_]SimTKmolmodel[_d].lib
#   [ns_]SimTKmolmodel_static[_d].lib
# Unix:
#   lib[ns_]SimTKmolmodel[_d].so
#   lib[ns_]SimTKmolmodel_static[_d].a
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

CMAKE_MINIMUM_REQUIRED(VERSION 3.2)

PROJECT(Molmodel)

SET(MOLMODEL_MAJOR_VERSION 3)
SET(MOLMODEL_MINOR_VERSION 0)
SET(MOLMODEL_PATCH_VERSION 0)

SET(MOLMODEL_COPYRIGHT_YEARS "2006-12")

# underbar separated list of dotted authors, no spaces or commas
set(MOLMODEL_AUTHORS "Christopher.Bruns_Michael.Sherman")
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

option(ASSUME_VERSIONED_SIMBODY "Assume that Simbody libraries were built with the version tag appended to the libraries name" OFF)

# Report the version number to the CMake UI. Don't include the
# build version if it is zero.
SET(PATCH_VERSION_STRING)
IF(MOLMODEL_PATCH_VERSION)
    SET(PATCH_VERSION_STRING ".${MOLMODEL_PATCH_VERSION}")
ENDIF()

set(
    MOLMODEL_VERSION
    "${MOLMODEL_MAJOR_VERSION}.${MOLMODEL_MINOR_VERSION}${MOLMODEL_PATCH_VERSION}"
)

set(
    MOLMODEL_SONAME_VERSION
    "${MOLMODEL_MAJOR_VERSION}.${MOLMODEL_MINOR_VERSION}"
)

# Permit use of custom FindOpenMM and FindSimbody modules
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")

# Find Simbody using the local FindSimbody.cmake if present.
find_package(Simbody REQUIRED)
find_package(OpenMM REQUIRED)
find_package(ZLIB REQUIRED)

include_directories(SYSTEM ${Simbody_INCLUDE_DIR})
include_directories(SYSTEM ${OpenMM_INCLUDE_DIR})
include_directories(SYSTEM ${ZLIB_INCLUDE_DIRS})

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

SET(NS)
IF (BUILD_USING_NAMESPACE)
    SET(NS "${BUILD_USING_NAMESPACE}_")
ENDIF ()

SET(MOLMODEL_LIBRARY_NAME ${NS}SimTKmolmodel CACHE STRING
    "Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."
    FORCE)

if (ASSUME_VERSIONED_SIMBODY)
    set(SIMBODY_VERTAG _${Simbody_VERSION})
endif ()

SET(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Control volume of build output" )

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
    ADD_DEFINITIONS(-DNOMINMAX -DWIN32_LEAN_AND_MEAN -D_USE_MATH_DEFINES)

	IF(MSVC)
	    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /MP")
	ENDIF()
ENDIF()

SET(BUILD_PLATFORM "${CMAKE_HOST_SYSTEM_NAME}:${PLATFORM_ABI}" CACHE STRING
    "This is the platform and ABI we're building for. Not changeable here; use a different CMake generator instead."
    FORCE)

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

ADD_DEFINITIONS(-DSimTK_MOLMODEL_LIBRARY_NAME=${MOLMODEL_LIBRARY_NAME}
                -DSimTK_MOLMODEL_MAJOR_VERSION=${MOLMODEL_MAJOR_VERSION}
                -DSimTK_MOLMODEL_MINOR_VERSION=${MOLMODEL_MINOR_VERSION}
                -DSimTK_MOLMODEL_PATCH_VERSION=${MOLMODEL_PATCH_VERSION})

ADD_DEFINITIONS(-DSimTK_MOLMODEL_SVN_REVISION="${MOLMODEL_SVN_REVISION}"
                -DSimTK_MOLMODEL_COPYRIGHT_YEARS="${MOLMODEL_COPYRIGHT_YEARS}"
                -DSimTK_MOLMODEL_AUTHORS="${MOLMODEL_AUTHORS}")

# -DSimTK_MOLMODEL_TYPE has to be defined in the target subdirectories.
# -Dmolmodel_EXPORTS defined automatically when Windows DLL build is being done.

SET(SHARED_TARGET ${MOLMODEL_LIBRARY_NAME})
SET(STATIC_TARGET ${MOLMODEL_LIBRARY_NAME}_static)

set(CMAKE_DEBUG_POSTFIX "_d")

SET(TEST_SHARED_TARGET ${SHARED_TARGET} ${Simbody_LIBRARIES})
SET(TEST_STATIC_TARGET ${STATIC_TARGET} ${Simbody_STATIC_LIBRARIES})

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

if (NOT ${GEMMI_DIR} STREQUAL "")
    include_directories("${GEMMI_DIR}/include")
endif ()

#
# Allow automated build and dashboard.
#

INCLUDE (Dart)
if (BUILD_TESTING)
    ENABLE_TESTING()
endif ()

# all the various variables have been SET above.

# Fixup Simbody library names in case we are building against versioned Simbody libraries
if (ASSUME_VERSIONED_SIMBODY)
    foreach(SL ${Simbody_LIBRARIES})
        list(APPEND SIMBODY_FIXED_LIBS ${SL}${SIMBODY_VERTAG})
        list(APPEND SIMBODY_FIXED_STATIC_LIBS ${SL}${SIMBODY_VERTAG}_static)
    endforeach()
else ()
    set(SIMBODY_FIXED_LIBS ${Simbody_LIBRARIES})
    set(SIMBODY_FIXED_STATIC_LIBS ${Simbody_STATIC_LIBRARIES})
endif ()

# Set installation directories
if (UNIX)
    include(GNUInstallDirs)
elseif (WIN32)
    set(CMAKE_INSTALL_LIBDIR "lib")
    set(CMAKE_INSTALL_INCLUDEDIR "include")
    set(CMAKE_INSTALL_BINDIR "bin")
endif ()

if (BUILD_SHARED_LIBRARIES)
    add_library(
        ${SHARED_TARGET} SHARED
        ${SOURCE_FILES}
        ${SOURCE_INCLUDE_FILES}
        ${API_ABS_INCLUDE_FILES}
    )
    set_target_properties(
        ${SHARED_TARGET}
        PROPERTIES
            PROJECT_LABEL "Library - ${SHARED_TARGET}"
    )
    target_compile_definitions(${SHARED_TARGET} PRIVATE -DSimTK_MOLMODEL_BUILDING_SHARED_LIBRARY)

    if (WIN32)
        target_link_libraries(
            ${SHARED_TARGET}
            PRIVATE ${SIMBODY_FIXED_LIBS}
            PRIVATE ${ZLIB_LIBRARIES}
            PRIVATE ws2_32
        )
    else ()
        target_link_libraries(
            ${SHARED_TARGET}
            PRIVATE ${SIMBODY_FIXED_LIBS}
            PRIVATE ${ZLIB_LIBRARIES}
        )
    endif ()

    install(
        TARGETS ${SHARED_TARGET} EXPORT MolmodelTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
endif ()

if (BUILD_STATIC_LIBRARIES)
    add_library(
        ${STATIC_TARGET} STATIC
        ${SOURCE_FILES}
        ${SOURCE_INCLUDE_FILES}
        ${API_ABS_INCLUDE_FILES}
    )
    set_target_properties(
        ${STATIC_TARGET}
        PROPERTIES
            PROJECT_LABEL "Library - ${STATIC_TARGET}"
    )
    target_compile_definitions(${STATIC_TARGET} PRIVATE -DSimTK_MOLMODEL_BUILDING_STATIC_LIBRARY -DSimTK_USE_STATIC_LIBRARIES)

    if (WIN32)
        target_link_libraries(
            ${STATIC_TARGET}
            PRIVATE ${SIMBODY_FIXED_STATIC_LIBS}
            PRIVATE ${ZLIB_LIBRARIES}
            PRIVATE ws2_32
        )
    else ()
        target_link_libraries(
            ${STATIC_TARGET}
            PRIVATE ${SIMBODY_FIXED_STATIC_LIBS}
            PRIVATE ${ZLIB_LIBRARIES}
        )
    endif ()

    install(
        TARGETS ${STATIC_TARGET} EXPORT MolmodelTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    )
endif ()

IF(BUILD_EXAMPLES)
    ADD_SUBDIRECTORY(examples)
ENDIF()

IF(BUILD_TESTING)
    ADD_SUBDIRECTORY(tests)
ENDIF()

if (WIN32)
    # Install include files into base include folder since it's a sandbox
    set(MOLMODEL_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR})
else ()
    # Install include files in simbody subfolder to avoid polluting the
    # global build folder
    set(MOLMODEL_INCLUDE_INSTALL_DIR ${CMAKE_INSTALL_INCLUDEDIR}/molmodel)
endif ()

INCLUDE(ApiDoxygen.cmake)
# libraries are installed from their subdirectories; headers here
# install headers
FILE(GLOB CORE_HEADERS     include/*.h                   */include/*.h)
FILE(GLOB TOP_HEADERS      include/molmodel/*.h          */include/molmodel/*.h)
FILE(GLOB INTERNAL_HEADERS include/molmodel/internal/*.h */include/molmodel/internal/*.h)
install(FILES ${CORE_HEADERS} DESTINATION ${MOLMODEL_INCLUDE_INSTALL_DIR})
install(FILES ${TOP_HEADERS} DESTINATION ${MOLMODEL_INCLUDE_INSTALL_DIR}/molmodel)
install(FILES ${INTERNAL_HEADERS} DESTINATION ${MOLMODEL_INCLUDE_INSTALL_DIR}/molmodel/internal)

# Install documents.
FILE(GLOB TOPLEVEL_DOCS doc/*.pdf doc/*.txt)
INSTALL(FILES ${TOPLEVEL_DOCS} DESTINATION doc)

set(PKG_NAME ${PROJECT_NAME})
set(PKG_LIBRARIES ${MOLMODEL_LIBRARY_NAME})

if (WIN32)
    set(MOLMODEL_CMAKE_DIR cmake)
else ()
    set(MOLMODEL_CMAKE_DIR ${CMAKE_INSTALL_LIBDIR}/cmake/molmodel/)
endif ()

include(CMakePackageConfigHelpers)
configure_package_config_file(
    ${CMAKE_SOURCE_DIR}/cmake/MolmodelConfig.cmake.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/MolmodelConfigForInstall.cmake
    INSTALL_DESTINATION "${MOLMODEL_CMAKE_DIR}"
    PATH_VARS
        CMAKE_INSTALL_PREFIX
        MOLMODEL_INCLUDE_INSTALL_DIR
        CMAKE_INSTALL_LIBDIR
        CMAKE_INSTALL_BINDIR
)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/MolmodelConfigForInstall.cmake
    DESTINATION ${MOLMODEL_CMAKE_DIR}
    RENAME MolmodelConfig.cmake
)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/MolmodelConfigVersion.cmake
    VERSION "${MOLMODEL_VERSION}"
    COMPATIBILITY SameMinorVersion
)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/MolmodelConfigVersion.cmake
    DESTINATION ${MOLMODEL_CMAKE_DIR}
)

install(EXPORT MolmodelTargets DESTINATION "${MOLMODEL_CMAKE_DIR}")

if (WIN32)
    # Nothing specific is required on Win32
elseif (APPLE)
    set(PKGCONFIG_PLATFORM_LIBS "-ldl")
else ()
    set(PKGCONFIG_PLATFORM_LIBS "-lpthread -lrt -ldl -lm")
endif ()

configure_file(
    ${CMAKE_SOURCE_DIR}/cmake/pkgconfig/molmodel.pc.in
    ${CMAKE_CURRENT_BINARY_DIR}/cmake/pkgconfig/molmodel.pc @ONLY
)
install(
    FILES ${CMAKE_CURRENT_BINARY_DIR}/cmake/pkgconfig/molmodel.pc
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/pkgconfig/
)

# Build the OpenMM plugin
add_subdirectory(OpenMMPlugin)
