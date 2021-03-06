# Build the Molmodel OpenMM Plugin

find_package(OpenMM REQUIRED)
INCLUDE_DIRECTORIES(${OpenMM_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

# The library is dependent on its local .cpp file and some of
# the Molmodel classes. We'll list everything to be safe.

IF (OpenMM_STATIC_LIBRARY)
    ADD_LIBRARY(OpenMMPlugin STATIC
        OpenMMPlugin.cpp
        ${SOURCE_INCLUDE_FILES}
        ${API_ABS_INCLUDE_FILES})
ELSE ()
    ADD_LIBRARY(OpenMMPlugin SHARED
        OpenMMPlugin.cpp
        ${SOURCE_INCLUDE_FILES}
        ${API_ABS_INCLUDE_FILES})
ENDIF ()

IF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)
    set_target_properties(OpenMMPlugin PROPERTIES DEBUG_POSTFIX "_d")
ENDIF (UNIX AND CMAKE_BUILD_TYPE MATCHES Debug)

# The plugin depends on the OpenMM.dll that is part of
# the OpenMM distribution.

IF (BUILD_UNVERSIONED_LIBRARIES)
    TARGET_LINK_LIBRARIES(OpenMMPlugin
        ${SHARED_TARGET}			   # SimTKMolmodel.dll
		SimTKcommon
        ${OpenMM_LIBRARIES})
ELSEIF(BUILD_VERSIONED_LIBRARIES)
    TARGET_LINK_LIBRARIES(OpenMMPlugin
        ${SHARED_TARGET_VN}			   # SimTKMolmodel.dll
		SimTKcommon_${Simbody_VERSION}
        ${OpenMM_LIBRARIES})
ELSE()
    MESSAGE(FATAL_ERROR "SimTKmolmodel was not configured to build")
ENDIF()
    
SET_TARGET_PROPERTIES(OpenMMPlugin PROPERTIES
	PROJECT_LABEL "Plugin - OpenMMPlugin")

# Install the plugin with the Molmodel distribution; 
# on Windows both .lib and .dll go in the lib/plugins
# directory although the ".lib" should not be needed for
# a runtime-linked DLL.
INSTALL_TARGETS(/lib${LIB64}/plugins 
                RUNTIME_DIRECTORY /lib${LIB64}/plugins
                OpenMMPlugin)
