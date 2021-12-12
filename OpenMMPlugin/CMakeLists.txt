# Build the Molmodel OpenMM Plugin

find_package(OpenMM REQUIRED)
INCLUDE_DIRECTORIES(${OpenMM_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

set(CMAKE_DEBUG_POSTFIX "_d")

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

SET_TARGET_PROPERTIES(OpenMMPlugin PROPERTIES
	PROJECT_LABEL "Plugin - OpenMMPlugin")

# The plugin depends on the OpenMM.dll that is part of
# the OpenMM distribution.

TARGET_LINK_LIBRARIES(
    OpenMMPlugin
    PRIVATE ${SHARED_TARGET}        # SimTKMolmodel.dll
    PRIVATE ${SIMBODY_FIXED_LIBS}   # We actually just need SimTKcommon but there is no easy way how to extract it from the list of all libs
    PRIVATE ${OpenMM_LIBRARIES}
)


install(
    TARGETS OpenMMPlugin EXPORT MolmodelTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/plugins
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/plugins
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/plugins
)
