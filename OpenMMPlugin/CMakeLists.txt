# Build the Molmodel OpenMM Plugin

find_package(OpenMM REQUIRED)
INCLUDE_DIRECTORIES(${OpenMM_INCLUDE_DIR})
INCLUDE_DIRECTORIES(${CMAKE_SOURCE_DIR}/src)

set(CMAKE_DEBUG_POSTFIX "_d")

# The library is dependent on its local .cpp file and some of
# the Molmodel classes. We'll list everything to be safe.

if (BUILD_STATIC_LIBRARIES)
    ADD_LIBRARY(OpenMMPlugin_static STATIC
        OpenMMPlugin.cpp
        ${SOURCE_INCLUDE_FILES}
        ${API_ABS_INCLUDE_FILES})

    SET_TARGET_PROPERTIES(OpenMMPlugin_static
                          PROPERTIES
                          OUTPUT_NAME OpenMMPlugin
                              PROJECT_LABEL "Plugin - OpenMMPlugin (static)")

    TARGET_LINK_LIBRARIES(
        OpenMMPlugin_static
        PRIVATE ${STATIC_TARGET}                # SimTKMolmodel.dll
        PRIVATE ${SIMBODY_FIXED_STATIC_LIBS}    # We actually just need SimTKcommon but there is no easy way how to extract it from the list of all libs
        PRIVATE ${OpenMM_LIBRARIES}
    )

    install(
        TARGETS OpenMMPlugin_static EXPORT MolmodelTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/plugins
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/plugins
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/plugins
    )
endif ()
if (BUILD_SHARED_LIBRARIES)
    ADD_LIBRARY(OpenMMPlugin_shared SHARED
        OpenMMPlugin.cpp
        ${SOURCE_INCLUDE_FILES}
        ${API_ABS_INCLUDE_FILES})

    SET_TARGET_PROPERTIES(OpenMMPlugin_shared
                          PROPERTIES
                              OUTPUT_NAME OpenMMPlugin
                              PROJECT_LABEL "Plugin - OpenMMPlugin")

    TARGET_LINK_LIBRARIES(
        OpenMMPlugin_shared
        PRIVATE ${SHARED_TARGET}         # SimTKMolmodel.dll
        PRIVATE ${SIMBODY_FIXED_LIBS}    # We actually just need SimTKcommon but there is no easy way how to extract it from the list of all libs
        PRIVATE ${OpenMM_LIBRARIES}
    )

    install(
        TARGETS OpenMMPlugin_shared EXPORT MolmodelTargets
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}/plugins
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}/plugins
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}/plugins
    )
endif ()
