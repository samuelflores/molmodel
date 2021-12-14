# Molmodel
Molmodel is the C++ molecular-modeling front-end for Simbody.

## Build instructions
Molmodel uses [CMake](https://cmake.org/) build system to configure and build itself. Molmodel provides various build-time configuration options. Depending on how your particular system is set up, some of these options may have to be adjusted in order to build Molmodel.

Molmodel depends on a series of libraries that must be installed in order to build and use Molmodel. Current list of Molmodel dependencies is as follows:

- [zlib](https://zlib.net/)
- [Simbody](https://simtk.org/projects/simbody)
- [OpenMM](https://openmm.org/)
- [Gemmi](https://gemmi.readthedocs.io/en/latest/)

### List of configuration options
- `BUILD_SHARED_LIBRARIES` [`TRUE`|`FALSE`] - Whether to build Molmodel as a dynamically linked (shared) library. Defaults to `TRUE`
- `BUILD_STATIC_LIBRARIES` [`TRUE`|`FALSE`] - Whether to build Molmodel as a statically linked library. Defaults to `FALSE`

Note that at least one of the options above shall be enabled

- `BUILD_TESTING` [`ON`|`OFF`] - Whether to build unit tests. Defaults to `ON`
- `BUILD_EXAMPLES` [`ON`|`OFF`] - Whether to build example code. Defaults to `ON`
- `BUILD_TESTING_STATIC` [`ON`|`OFF`] - Build tests that link against shared Molmodel library. This requires `BUILD_SHARED_LIBRARIES` to be set to `TRUE`
- `BUILD_TESTING_STATIC` [`ON`|`OFF`] - Build tests that link against static Molmodel library. This requires `BUILD_STATIC_LIBRARIES` to be set to `TRUE`

- `ASSUME_VERSIONED_SIMBODY` [`ON|OFF`] - Simbody libraries may be built with the Simbody version tag appended to the file names of the libraries. If your Simbody installation is built like this, enable this option to allow MMB to link properly against Simbody. Defaults to `OFF`

- `GEMMI_DIR` - Path to [Gemmi](https://gemmi.readthedocs.io/en/latest/) library installation. This shall be specified only if Gemmi is not installed as a system-wide library.

### Build and installation

Under most circumstances, Molmodel should be able to set itself up for build automatically. It may be necessary to pass some additional configuration parameters to CMake if some of Molmodel dependencies are not installed in standard system paths.

#### UNIX systems

To build Simbody on a UNIX system, `cd` into the directory with the Molmodel source code and issue the following commands

    mkdir build
    cd build
    cmake .. -DCMAKE_BUILD_TYPE=Release
    make
    make install

note that the `make install` command may require root privileges.

If either Simbody or OpenMM are not installed in standard system paths, it is necessary to tell CMake where to look for them. This can be done with the `CMAKE_PREFIX_PATH` variable. For example, if your Simbody is installed in `/opt/simbody`, use the following command to point CMake to the correct directory

    cmake .. -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH=/opt/simbody

If you need to specify more than one path, use the `;` (semicolon) to input multiple paths, i.e.

    -DCMAKE_PREFIX_PATH=/path/to/simbody;/path/to/openmm

#### Windows systems

Configuration and build process on Windows is similar to that on UNIX systems. Since there are no standard system paths on Windows, it is necessary to specify paths to all Molmodel dependencies. Note that CMake for Windows comes with a graphical user interface to input the configuration options.
