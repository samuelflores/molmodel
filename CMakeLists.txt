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





cmake_minimum_required(VERSION 2.8)





project(Molmodel)





set(MOLMODEL_MAJOR_VERSION 3)


set(MOLMODEL_MINOR_VERSION 0)


set(MOLMODEL_PATCH_VERSION 0)





set(MOLMODEL_COPYRIGHT_YEARS "2006-12")





# underbar separated list of dotted authors, no spaces or commas


set(MOLMODEL_AUTHORS "Christopher.Bruns_Michael.Sherman")





# Report the version number to the CMake UI. Don't include the 


# build version if it is zero.


set(PATCH_VERSION_STRING)


IF(MOLMODEL_PATCH_VERSION)


    SET(PATCH_VERSION_STRING ".${MOLMODEL_PATCH_VERSION}")


ENDIF()





set(MOLMODEL_VERSION 


    "${MOLMODEL_MAJOR_VERSION}.${MOLMODEL_MINOR_VERSION}${PATCH_VERSION_STRING}"


    CACHE STRING 


    "This is the version that will be built (can't be changed here)." 


    FORCE)


#set(MOLMODEL_SONAME_VERSION

#    "${MOLMODEL_MAJOR_VERSION}.${MOLMODEL_MINOR_VERSION}"

#    CACHE STRING

#    "Soname version; appended to names of shared libs

#    (can't be changed in GUI)."

#    FORCE)



#set(SONAME

#    "${MOLMODEL_MAJOR_VERSION}.${MOLMODEL_MINOR_VERSION}"

#    CACHE STRING

#    "Soname version; appended to names of shared libs

#    (can't be changed in GUI)."

#    FORCE)






# This is the suffix if we're building and depending on versioned libraries.


set(VN "_${MOLMODEL_VERSION}")





set(BUILD_BINARY_DIR ${CMAKE_BINARY_DIR}


    CACHE PATH 


    "The Molmodel build (not the install) puts all the libraries and executables together here (with /Release, etc. appended on some platforms).")





# Make everything go in the same binary directory. (These are CMake-defined


# variables.)


set(EXECUTABLE_OUTPUT_PATH ${BUILD_BINARY_DIR})


set(LIBRARY_OUTPUT_PATH ${BUILD_BINARY_DIR})



SET(BUILD_SHARED_LIBRARIES TRUE CACHE BOOL
    "Build dynamically linked (shared) version of the libraries")


# Static libraries, tests, and examples won't be built unless this 


# is set.


SET(BUILD_STATIC_LIBRARIES FALSE CACHE BOOL
    "Build '_static' versions of the libraries libraries?")

IF((NOT ${BUILD_SHARED_LIBRARIES}) AND (NOT ${BUILD_STATIC_LIBRARIES}))
    MESSAGE(FATAL_ERROR "Neither shared nor static build has been enabled. Aborting")
ENDIF()



# Use this to generate a private set of libraries whose names


# won't conflict with installed versions.


set(BUILD_USING_NAMESPACE "" CACHE STRING


	"All library names will be prefixed with 'xxx_' if this is set to xxx.")





set(BUILD_UNVERSIONED_LIBRARIES TRUE CACHE BOOL


 "Build library names, and assume dependency names, with no version numbers?")





set(BUILD_VERSIONED_LIBRARIES FALSE CACHE BOOL


 "Build library names, and assume dependency names, with version numbers?")





set(NS)


if(BUILD_USING_NAMESPACE)


    set(NS "${BUILD_USING_NAMESPACE}_")


endif()





set(MOLMODEL_LIBRARY_NAME ${NS}SimTKmolmodel CACHE STRING


"Base name of the library being built; can't be changed here; see BUILD_USING_NAMESPACE variable."


FORCE)








# Permit use of custom FindOpenMM and FindSimbody modules


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake_modules")





set(CMAKE_VERBOSE_MAKEFILE OFF CACHE BOOL "Control volume of build output" )





set(OpenMM_ON OFF)


find_package(OpenMM)


if(OpenMM_FOUND)


    set(OpenMM_ON ON)


endif()


set(MOLMODEL_USE_OpenMM ${OpenMM_ON} CACHE BOOL


    "Control whether Molmodel builds the OpenMM Plugin (requires that OpenMM is installed on the build machine)." )





# Caution: this variable is automatically created by the CMake


# ENABLE_TESTING() command, but we'll take it over here for


# our own purposes too.


set(BUILD_TESTING ON CACHE BOOL


	"Control building of Molmodel test programs." )





set(BUILD_EXAMPLES ON CACHE BOOL


	"Control building of Molmodel example programs." )





# Turning this off reduces the build time (and space) substantially,


# but you may miss the occasional odd bug. Also currently on Windows it


# is easier to debug the static tests than the DLL-liked ones.


SET( BUILD_TESTING_STATIC ON CACHE BOOL


    "If building static libraries, build static test and example programs too?" )





SET( BUILD_TESTING_SHARED ON CACHE BOOL


    "If building test or example programs, include dynamically-linked ones?" )





#


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


# 32 bit Intel instruction set, or x64 for 64 bit Intel instruction set.





IF(${CMAKE_SIZEOF_VOID_P} EQUAL 8)


    SET(PLATFORM_ABI x64)


ELSE()


    SET(PLATFORM_ABI x86)


ENDIF()





SET(BUILD_PLATFORM "${CMAKE_HOST_SYSTEM_NAME}:${PLATFORM_ABI}" CACHE STRING


    "This is the platform and ABI we're building for. Not changeable here; use a different CMake generator instead."


    FORCE)





# Find Simbody using the local FindSimbody.cmake if present.


find_package(Simbody REQUIRED)


include_directories(${Simbody_INCLUDE_DIR})


link_directories(${Simbody_LIB_DIR})





# If CMAKE_INSTALL_PREFIX is /usr/local, then the LIBDIR should necessarily be


# lib/. Sometimes (on Linux), LIBDIR is something like x86_64-linux-gnu. The


# linker will search /usr/lib/x86_64-linux-gnu (this path is in


# /etc/ld.so.conf.d), but it will NOT search /usr/local/lib/x86-64-linux-gnu.


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





## Choose the maximum level of x86 instruction set that the compiler is 


## allowed to use. 


## Was using sse2 but changed to let the compilers choose. Most will


## probably use sse2 or later by default.


## On 64 bit MSVC 2013, the default is sse2 and the argument


## isn't recognized so don't specify it.


if (CMAKE_CL_64)


    set(default_build_inst_set)


else()


    set(default_build_inst_set)


endif()





## This can be set to a different value by the person running CMake.


SET(BUILD_INST_SET ""


    CACHE STRING 


    "CPU instruction level compiler is permitted to use (default: let compiler decide).")


MARK_AS_ADVANCED( BUILD_INST_SET )





if (BUILD_INST_SET)


    set(inst_set_to_use ${BUILD_INST_SET})


else()


    set(inst_set_to_use ${default_build_inst_set})


endif()





## When building in any of the Release modes, tell gcc/clang to use 


## not-quite most agressive optimization.  Here we 


## are specifying *all* of the Release flags, overriding CMake's defaults.


## Watch out for optimizer bugs in particular gcc versions!





IF(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU" OR


   ${CMAKE_CXX_COMPILER_ID} MATCHES "Clang")





    # If using either of these compilers, provide the option of 


    # compiling using the c++11 standard.


    OPTION(SIMBODY_STANDARD_11


        "Compile using the C++11 standard, if using GCC or Clang." ON)


    


    IF(${SIMBODY_STANDARD_11})


        # Using C++11 on OSX requires using libc++ instead of libstd++.


        # libc++ is an implementation of the C++ standard library for OSX.


        IF(APPLE)


            IF(XCODE)


                SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LANGUAGE_STANDARD "c++11")


                SET(CMAKE_XCODE_ATTRIBUTE_CLANG_CXX_LIBRARY "libc++")


            ELSE()


                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


                SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")


            ENDIF()


        ELSE() # not APPLE


            SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")


        ENDIF()


    ENDIF()





    if (inst_set_to_use)


        string(TOLOWER ${inst_set_to_use} GCC_INST_SET)


        set(GCC_INST_SET "-m${GCC_INST_SET}")


    else()


        set(GCC_INST_SET)


    endif()





    # Get the gcc or clang version number in major.minor.build format


    execute_process(COMMAND ${CMAKE_C_COMPILER} -dumpversion


                    OUTPUT_VARIABLE GCC_VERSION)





    # Unrolling fixed-count loops was a useful optimization for Simmatrix


    # in earlier gcc versions.


    # Doesn't have a big effect for current compiler crop and may be 


    # pushing our luck with optimizer bugs. So let the compilers decide


    # how to handle loops instead.


    ##SET(GCC_OPT_ENABLE "-funroll-loops")





    # If you know of optimization bugs that affect Simbody in particular


    # gcc versions, this is the place to turn off those optimizations.


    SET(GCC_OPT_DISABLE)


    # We know Gcc 4.4.3 on Ubuntu 10 is buggy.


    IF(${CMAKE_CXX_COMPILER_ID} MATCHES "GNU")


      if (GCC_VERSION VERSION_EQUAL 4.4)


        SET(GCC_OPT_DISABLE 


        "-fno-strict-aliasing -fno-tree-vrp -fno-guess-branch-probability")


      endif()


    ENDIF()





    # C++


    SET(BUILD_CXX_FLAGS_DEBUG          "-g ${GCC_INST_SET}")


    SET(BUILD_CXX_FLAGS_RELEASE        


      "-DNDEBUG -O2 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")


    SET(BUILD_CXX_FLAGS_RELWITHDEBINFO 


    "-DNDEBUG -O2 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")


    SET(BUILD_CXX_FLAGS_MINSIZEREL     "-DNDEBUG -Os ${GCC_INST_SET}")





    # C


    SET(BUILD_C_FLAGS_DEBUG            "-g ${GCC_INST_SET}")


    SET(BUILD_C_FLAGS_RELEASE          


      "-DNDEBUG -O2 ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")


    SET(BUILD_C_FLAGS_RELWITHDEBINFO   


    "-DNDEBUG -O2 -g ${GCC_OPT_ENABLE} ${GCC_OPT_DISABLE} ${GCC_INST_SET}")


    SET(BUILD_C_FLAGS_MINSIZEREL       "-DNDEBUG -Os ${GCC_INST_SET}")





    # C++


    SET(CMAKE_CXX_FLAGS_DEBUG ${BUILD_CXX_FLAGS_DEBUG}


        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)


    SET(CMAKE_CXX_FLAGS_RELEASE ${BUILD_CXX_FLAGS_RELEASE}


        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)


    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO ${BUILD_CXX_FLAGS_RELWITHDEBINFO}


        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)


    SET(CMAKE_CXX_FLAGS_MINSIZEREL ${BUILD_CXX_FLAGS_MINSIZEREL}


        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)





    # C


    SET(CMAKE_C_FLAGS_DEBUG ${BUILD_C_FLAGS_DEBUG}


        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)


    SET(CMAKE_C_FLAGS_RELEASE ${BUILD_C_FLAGS_RELEASE}         


        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)


    SET(CMAKE_C_FLAGS_RELWITHDEBINFO ${BUILD_C_FLAGS_RELWITHDEBINFO}


        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)


    SET(CMAKE_C_FLAGS_MINSIZEREL ${BUILD_C_FLAGS_MINSIZEREL}


        CACHE STRING "Can't change here -- see BUILD_C..." FORCE)





ENDIF()





## When building in any of the Release modes, tell VC++ cl compiler to use 


## intrinsics (i.e. sqrt instruction rather than sqrt subroutine) with 


## flag /Oi.


## Caution: can't use CMAKE_CXX_COMPILER_ID MATCHES MSVC here because


## "MSVC" is a predefined CMAKE variable and will get expanded to 1 or 0


IF(MSVC)


    if (inst_set_to_use)


        string(TOUPPER ${inst_set_to_use} CL_INST_SET)


        set(CL_INST_SET "/arch:${CL_INST_SET}")


    else()


        set(CL_INST_SET)


    endif()





    set(BUILD_LIMIT_PARALLEL_COMPILES "" CACHE STRING


        "Set a maximum number of simultaneous compilations.")


    mark_as_advanced(BUILD_LIMIT_PARALLEL_COMPILES)


    set(mxcpu ${BUILD_LIMIT_PARALLEL_COMPILES}) # abbreviation





    ## C++


    SET(BUILD_CXX_FLAGS_DEBUG        


    "/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- ${CL_INST_SET}")


    SET(BUILD_CXX_FLAGS_RELEASE        


    "/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- ${CL_INST_SET}")


    SET(BUILD_CXX_FLAGS_RELWITHDEBINFO 


    "/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- ${CL_INST_SET}")


    SET(BUILD_CXX_FLAGS_MINSIZEREL 


    "/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- ${CL_INST_SET}")





    ## C


    SET(BUILD_C_FLAGS_DEBUG        


    "/D _DEBUG /MDd /Od /Ob0 /RTC1 /Zi /GS- ${CL_INST_SET}")


    SET(BUILD_C_FLAGS_RELEASE        


    "/D NDEBUG /MD  /O2 /Ob2 /Oi /GS- ${CL_INST_SET}")


    SET(BUILD_C_FLAGS_RELWITHDEBINFO 


    "/D NDEBUG /MD  /O2 /Ob2 /Oi /Zi /GS- ${CL_INST_SET}")


    SET(BUILD_C_FLAGS_MINSIZEREL 


    "/D NDEBUG /MD  /O1 /Ob1 /Oi /GS- ${CL_INST_SET}")





    ## C++


    SET(CMAKE_CXX_FLAGS_DEBUG "/MP${mxcpu} ${BUILD_CXX_FLAGS_DEBUG}"


        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)


    SET(CMAKE_CXX_FLAGS_RELEASE "/MP${mxcpu} ${BUILD_CXX_FLAGS_RELEASE}"


        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)


    SET(CMAKE_CXX_FLAGS_RELWITHDEBINFO "/MP${mxcpu} ${BUILD_CXX_FLAGS_RELWITHDEBINFO}"


        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)


    SET(CMAKE_CXX_FLAGS_MINSIZEREL "/MP${mxcpu} ${BUILD_CXX_FLAGS_MINSIZEREL}"


        CACHE STRING "Can't change here -- see BUILD_CXX..." FORCE)





    ## C


    SET(CMAKE_C_FLAGS_DEBUG "/MP${mxcpu} ${BUILD_C_FLAGS_DEBUG}"


        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)


    SET(CMAKE_C_FLAGS_RELEASE "/MP${mxcpu} ${BUILD_C_FLAGS_RELEASE}"


        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)


    SET(CMAKE_C_FLAGS_RELWITHDEBINFO "/MP${mxcpu} ${BUILD_C_FLAGS_RELWITHDEBINFO}"


        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)


    SET(CMAKE_C_FLAGS_MINSIZEREL "/MP${mxcpu} ${BUILD_C_FLAGS_MINSIZEREL}"


        CACHE STRING "Can't change here -- see BUILD_C_..." FORCE)





ENDIF()





# The source is organized into subdirectories, but we handle them all from


# this CMakeLists file rather than letting CMake visit them as SUBDIRS.


set(MOLMODEL_SOURCE_SUBDIRS .)








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


find_program(SVNVERSION_EXECUTABLE svnversion PATHS "C:/cygwin/bin")


if(SVNVERSION_EXECUTABLE)


    exec_program(${SVNVERSION_EXECUTABLE}


              # Works better on Windows to set working directory rather than


              # passing argument to svnversion


              ${CMAKE_SOURCE_DIR} # cwd for run


              OUTPUT_VARIABLE OUT)


    set(MOLMODEL_SVN_REVISION "${OUT}" CACHE STRING "Molmodel svn revision number" FORCE)


else(SVNVERSION_EXECUTABLE)


    message(STATUS 


         "Could not find 'svnversion' executable; 'about' will be wrong. (Cygwin provides one on Windows.)"


        )


    set(MOLMODEL_SVN_REVISION "unknown" CACHE STRING "Molmodel svn revision number")


endif(SVNVERSION_EXECUTABLE)


mark_as_advanced(MOLMODEL_SVN_REVISION)





# Remove colon from build version, for easier placement in directory names


STRING(REPLACE ":" "_" MOLMODEL_SVN_REVISION ${MOLMODEL_SVN_REVISION})





ADD_DEFINITIONS(-DSimTK_MOLMODEL_LIBRARY_NAME=${MOLMODEL_LIBRARY_NAME}


                -DSimTK_MOLMODEL_MAJOR_VERSION=${MOLMODEL_MAJOR_VERSION}


                -DSimTK_MOLMODEL_MINOR_VERSION=${MOLMODEL_MINOR_VERSION}


		-DSimTK_MOLMODEL_PATCH_VERSION=${MOLMODEL_PATCH_VERSION})





# CMake quotes automatically when building Visual Studio projects but we need


# to add them ourselves for Linux or Cygwin. Two cases to avoid duplicate quotes


# in Visual Studio which end up in the binary.





IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")


   SET(NEED_QUOTES FALSE)


ELSE (${CMAKE_GENERATOR} MATCHES "Visual Studio")


   SET(NEED_QUOTES TRUE)


ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")





##TODO: doesn't work without quotes in nightly build


set(NEED_QUOTES TRUE)





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





set(SHARED_TARGET ${MOLMODEL_LIBRARY_NAME})


set(STATIC_TARGET ${MOLMODEL_LIBRARY_NAME}_static)


set(SHARED_TARGET_VN ${MOLMODEL_LIBRARY_NAME}${VN})


set(STATIC_TARGET_VN ${MOLMODEL_LIBRARY_NAME}${VN}_static)





# Ensure that debug libraries have "_d" appended to their names.


# CMake gets this right on Windows automatically with this definition.


IF (${CMAKE_GENERATOR} MATCHES "Visual Studio")


    SET(CMAKE_DEBUG_POSTFIX "_d" CACHE INTERNAL "" FORCE)


ENDIF (${CMAKE_GENERATOR} MATCHES "Visual Studio")





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


	SET(TEST_SHARED_TARGET ${SHARED_TARGET})


	SET(TEST_STATIC_TARGET ${STATIC_TARGET})


ELSE(BUILD_UNVERSIONED_LIBRARIES)


	SET(TEST_SHARED_TARGET ${SHARED_TARGET_VN})


	SET(TEST_STATIC_TARGET ${STATIC_TARGET_VN})


ENDIF(BUILD_UNVERSIONED_LIBRARIES)





IF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")


    SET(ADDITIONAL_LINK_LIBRARIES)


ELSE(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")


    ## Assume Microsoft Visual Studio


	## ws2_32 is for socket calls in VMD connection


    SET(ADDITIONAL_LINK_LIBRARIES ws2_32)


ENDIF(${CMAKE_C_COMPILER} MATCHES "cc" OR ${CMAKE_C_COMPILER} MATCHES "clang")





SET ( USE_GEMMI  TRUE  CACHE BOOL "Should the Gemmi library be linked to the build?" )
SET ( GEMMI_PATH ""    CACHE PATH "Where is the Gemmi library installed?" )

IF ( USE_GEMMI )


	ADD_DEFINITIONS     ( -DGEMMI_USAGE )
	INCLUDE_DIRECTORIES ( ${GEMMI_PATH} )
	SET ( ADDITIONAL_LINK_LIBRARIES ${ADDITIONAL_LINK_LIBRARIES} "-lz" )


ENDIF ( USE_GEMMI )





# These are all the places to search for header files which are


# to be part of the API.


set(API_INCLUDE_DIRS) # start empty


FOREACH(subdir ${MOLMODEL_SOURCE_SUBDIRS})


    # append


    SET(API_INCLUDE_DIRS ${API_INCLUDE_DIRS}


                         ${subdir}/include 


                         ${subdir}/include/molmodel 


                         ${subdir}/include/molmodel/internal)


ENDFOREACH(subdir)





# We'll need both *relative* path names, starting with their API_INCLUDE_DIRS,


# and absolute pathnames.


set(API_REL_INCLUDE_FILES)   # start these out empty


set(API_ABS_INCLUDE_FILES)





FOREACH(dir ${API_INCLUDE_DIRS})


    FILE(GLOB fullpaths ${dir}/*.h)	# returns full pathnames


    SET(API_ABS_INCLUDE_FILES ${API_ABS_INCLUDE_FILES} ${fullpaths})





    FOREACH(pathname ${fullpaths})


        GET_FILENAME_COMPONENT(filename ${pathname} NAME)


        SET(API_REL_INCLUDE_FILES ${API_REL_INCLUDE_FILES} ${dir}/${filename})


    ENDFOREACH(pathname)


ENDFOREACH(dir)





# collect up source files


set(SOURCE_FILES) # empty


set(SOURCE_INCLUDE_FILES)





FOREACH(subdir ${MOLMODEL_SOURCE_SUBDIRS})


    FILE(GLOB src_files  ${subdir}/src/*.cpp ${subdir}/src/*.c


                         ${subdir}/src/*/*.cpp ${subdir}/src/*/*.c)


    FILE(GLOB incl_files ${subdir}/src/*.h ${subdir}/src/*/*.h )


    SET(SOURCE_FILES         ${SOURCE_FILES}         ${src_files})   #append


    SET(SOURCE_INCLUDE_FILES ${SOURCE_INCLUDE_FILES} ${incl_files})





    ## Make sure we find these locally before looking in


    ## SimTK/include if Molmodel was previously installed there.


    INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/${subdir}/include)


ENDFOREACH(subdir)





INCLUDE_DIRECTORIES(BEFORE ${PROJECT_SOURCE_DIR}/src)





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


INSTALL_FILES(/include/                 FILES ${CORE_HEADERS})


INSTALL_FILES(/include/molmodel/         FILES ${TOP_HEADERS})


INSTALL_FILES(/include/molmodel/internal FILES ${INTERNAL_HEADERS})





# Install documents.


FILE(GLOB TOPLEVEL_DOCS doc/*.pdf doc/*.txt)


INSTALL(FILES ${TOPLEVEL_DOCS} DESTINATION doc)





# These are at the end because we want them processed after


# all the various variables have been set above.





IF(BUILD_STATIC_LIBRARIES)


        ADD_SUBDIRECTORY( staticTarget )


ENDIF()

IF(BUILD_SHARED_LIBRARIES)
    ADD_SUBDIRECTORY( sharedTarget )
ENDIF(BUILD_SHARED_LIBRARIES)





IF (MOLMODEL_USE_OpenMM)


    ADD_SUBDIRECTORY( OpenMMPlugin )


ENDIF (MOLMODEL_USE_OpenMM)





IF( BUILD_EXAMPLES )


  ADD_SUBDIRECTORY( examples )


ENDIF()





IF( BUILD_TESTING )


  ADD_SUBDIRECTORY( tests )


ENDIF()








