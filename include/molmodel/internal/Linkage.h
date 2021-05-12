#ifndef _MOLMODEL_LINKAGE_H
#define _MOLMODEL_LINKAGE_H

#ifdef _WIN32
    #if defined(SimTK_MOLMODEL_BUILDING_SHARED_LIBRARY)
        #define SimTK_MOLMODEL_EXPORT __declspec(dllexport)
    #elif defined(SimTK_MOLMODEL_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SimTK_MOLMODEL_EXPORT
    #else
        #define SimTK_MOLMODEL_EXPORT __declspec(dllimport)   // i.e., a client of a shared library
    #endif
#else
    #if defined(SimTK_MOLMODEL_BUILDING_SHARED_LIBRARY)
        #define SimTK_MOLMODEL_EXPORT __attribute__ ((visibility ("default")))
    #else
        #define SimTK_MOLMODEL_EXPORT
    #endif
#endif

#endif // _MOLMODEL_LINKAGE_H
