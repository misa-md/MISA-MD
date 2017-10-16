cmake_minimum_required(VERSION 3.6)
project(CrystalMDproject)

set(CMAKE_CXX_STANDARD 11)

set(HEADER_FILES
        atom.h
        createatom.h
        domain.h
        domaindecomposition.h
        eam.h
        input.h
        integrator.h
        InterpolationObject.h
        latparticledata.h
        particledata.h
        simulation.h
        )

set(SOURCE_FILES
        atom.cpp
        createatom.cpp
        domain.cpp
        domaindecomposition.cpp
        eam.cpp
        input.cpp
        integrator.cpp
        InterpolationObject.cpp
        latparticledata.cpp
        particledata.cpp
        simulation.cpp
        )

set(MAIN_FILE
        main.cpp
        )

add_executable(CrystalMD ${MAIN_FILE} ${HEADER_FILES} ${SOURCE_FILES})