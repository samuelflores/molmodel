#ifndef SimTK_SIMBODY_VERSION_CHECK_H_
#define SimTK_SIMBODY_VERSION_CHECK_H_

#define SIMBODY_VERSION_CHECK(major, minor, patch) \
    ((major * 100000) + (minor * 1000) + patch)

#define SIMBODY_CURRENT_VERSION \
    ((SimTK_SimTKCOMMON_MAJOR_VERSION * 100000) + (SimTK_SimTKCOMMON_MINOR_VERSION * 1000) + SimTK_SimTKCOMMON_PATCH_VERSION)

#define SIMBODY_VERSION_STRING(major, minor, patch) \
    #major "." #minor "." #patch

#define SIMBODY_CURRENT_VERSION_STRING \
    SimTK_SimTKCOMMON_MAJOR_VERSION "." SimTK_SimTKCOMMON_MINOR_VERSION "." SimTK_SimTKCOMMON_PATCH_VERSION

#endif // SimTK_SIMBODY_VERSION_CHECK_H_
