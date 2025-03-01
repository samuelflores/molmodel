

#############
# New Dockerfile for MolModel        
# 
#############

#Download base image 
FROM ubuntu:22.04 as OSSetup

#############
#Get some libraries and executables
#############
# cmake requires build-essential libssl-dev
RUN apt update
#RUN apt upgrade
# This repo has python:
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt install -y git wget swig doxygen libblas-dev liblapack-dev  cmake-curses-gui zlib1g zlib1g-dev apt-utils snapd build-essential libssl-dev software-properties-common lsb-release vim python3.9
RUN mkdir /github
#############

#############
#for cmake:
#############
WORKDIR /github
RUN wget https://github.com/Kitware/CMake/archive/refs/tags/v3.25.0-rc3.tar.gz
RUN tar -zxvf v3.25.0-rc3.tar.gz
WORKDIR /github/CMake-3.25.0-rc3
# Build cmake
RUN ./bootstrap
RUN gmake -j8
RUN gmake install
#############

#############
# clone git repos
#############
WORKDIR /github
RUN git clone https://github.com/project-gemmi/gemmi.git /github/gemmi
#get a specific release of openmm:
RUN git clone https://github.com/pandegroup/openmm.git   /github/openmm --branch 7.7.0
RUN git clone  -b simbody-3.7 --single-branch https://github.com/simbody/simbody.git /github/simbody
RUN git clone https://github.com/seqan/seqan.git /github/seqan
# Clone specifically molmodel 3.1.0:
RUN git clone -b v3.1.0  https://github.com/samuelflores/molmodel.git /github/molmodel
# make build directories
RUN mkdir /github/gemmi/build
RUN mkdir /github/simbody/build
RUN mkdir /github/molmodel/build
RUN mkdir /github/openmm/build
#############

#############
# build gemmi    
#############
WORKDIR /github/gemmi/build
#latest build is not compatible with MMB. checkout version 5.6, 21 sept 2022:
RUN git checkout cca2ac864f47c9690aec73c7702bb5247d665f5a
WORKDIR /github/gemmi/build
RUN cmake ..
RUN make -j8
RUN make install
#############

#############
# Build simbody 
#############
WORKDIR /github/simbody/build
RUN cmake -DUSE_GEMMI=TRUE -DGEMMI_PATH=/github/gemmi/include -DSimbody_DIR=/usr/local/lib/cmake/simbody/ -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DBUILD_TESTING_SHARED=OFF -DBUILD_TESTING_STATIC=OFF ..
RUN make -j8
RUN make install
#############

#############
# build openmm   
#############
WORKDIR /github/openmm/build
RUN cmake ..
RUN make -j8
RUN make install
#############

#############
# Build molmodel
#############
WORKDIR  /github/molmodel/build
RUN cmake -DUSE_GEMMI=TRUE -DGEMMI_PATH=/github/gemmi/include -DSimbody_DIR=/usr/local/lib/cmake/simbody/ -DCMAKE_INSTALL_PREFIX=/usr/local -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=ON  -DBUILD_TESTING_SHARED=ON  -DBUILD_TESTING_STATIC=ON.  -DCMAKE_BUILD_TYPE=Release ..
RUN make -j8
RUN make install
#############


#############
# Set up environment
#############
ENV LD_LIBRARY_PATH=/usr/local/lib:/usr/local/openmm/lib
ENV PATH "/usr/local/bin:$PATH"
WORKDIR /work                  
#############
# Done!
#############
