## This whole directory exists just so I could define this extra 
## preprocessor value.
ADD_DEFINITIONS(-DSimTK_MOLMODEL_BUILDING_SHARED_LIBRARY)

IF(BUILD_UNVERSIONED_LIBRARIES)

    # These are dependencies so we know when we have to rebuild.
    ADD_LIBRARY(${SHARED_TARGET} SHARED 
                ${SOURCE_FILES} 
                ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
	IF (WIN32)
        TARGET_LINK_LIBRARIES(${SHARED_TARGET} 
            PRIVATE ${Simbody_LIBRARIES}
            PRIVATE zdll
			PRIVATE ws2_32)
	ELSE()
        TARGET_LINK_LIBRARIES(${SHARED_TARGET} 
            PRIVATE ${Simbody_LIBRARIES}
            z)
	ENDIF()
			
    SET_TARGET_PROPERTIES(${SHARED_TARGET} PROPERTIES
    	PROJECT_LABEL "Library - ${SHARED_TARGET}"
        SOVERSION "${MOLMODEL_MAJOR_VERSION}.${MOLMODEL_MINOR_VERSION}")
    
    # install library; on Windows .dll goes in the bin directory.
    INSTALL_TARGETS(/lib${LIB64} 
    		    RUNTIME_DIRECTORY /bin
                    ${SHARED_TARGET})

ENDIF(BUILD_UNVERSIONED_LIBRARIES)



IF(BUILD_VERSIONED_LIBRARIES)

    # These are dependencies so we know when we have to rebuild.
    ADD_LIBRARY(${SHARED_TARGET_VN} SHARED 
                ${SOURCE_FILES} 
                ${SOURCE_INCLUDE_FILES} 
                ${API_ABS_INCLUDE_FILES})
    
    IF (WIN32)
        TARGET_LINK_LIBRARIES(${SHARED_TARGET_VN} 
            PRIVATE ${Simbody_LIBRARIES}_${Simbody_VERSION}
            PRIVATE zdll
			PRIVATE ws2_32)
    ELSE()
            TARGET_LINK_LIBRARIES(${SHARED_TARGET_VN} 
            PRIVATE ${Simbody_LIBRARIES}_${Simbody_VERSION}
            PRIVATE z)
    ENDIF()

    SET_TARGET_PROPERTIES(${SHARED_TARGET_VN} PROPERTIES
    	PROJECT_LABEL "Library - ${SHARED_TARGET_VN}")
    
    # install library; on Windows .dll goes in the bin directory.
    INSTALL_TARGETS(/lib${LIB64} 
    		    RUNTIME_DIRECTORY /bin
                    ${SHARED_TARGET_VN})

ENDIF(BUILD_VERSIONED_LIBRARIES)
