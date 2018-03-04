
set(INCLUDE_FILES
        include/args.hpp
        include/toml.hpp)


set(HEADER_FILES
        atom.h
        createatom.h
        crystal_md.h
        config.h
        domain.h
        domaindecomposition.h
        eam.h
        input.h
        integrator.h
        InterpolationObject.h
        latparticledata.h
        mpi_utils.h
        particledata.h
        pre_config.h
        simulation.h
        )

set(SOURCE_FILES
        atom.cpp
        createatom.cpp
        crystal_md.cpp
        config.cpp
        domain.cpp
        domaindecomposition.cpp
        eam.cpp
        input.cpp
        integrator.cpp
        InterpolationObject.cpp
        latparticledata.cpp
        mpi_utils.cpp
        particledata.cpp
        simulation.cpp
        )

set(MAIN_FILE main.cpp)

add_executable(CrystalMD ${MAIN_FILE} ${INCLUDE_FILES} ${HEADER_FILES} ${SOURCE_FILES})

target_link_libraries(CrystalMD ${MPI_LIBRARIES}) # mpi lib
