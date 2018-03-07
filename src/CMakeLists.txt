
set(INCLUDE_FILES
        include/args.hpp
        include/toml.hpp)


set(HEADER_FILES
        arch_env.hpp
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
        hardware_accelerate.hpp)

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

# codes runs on sunway.
set(ARCH_SUNWAY_FILES
        arch/sunway/athread_accelerate.h
        arch/sunway/athread_accelerate.cpp
        arch/sunway/sunway_env.h
        arch/sunway/sunway_env.cpp
        arch/sunway/sunway_slave.c
        arch/sunway/sunway_slave2.c
        )

set(ARCH_FILES ${ARCH_SUNWAY_FILES})

set(MAIN_FILE main.cpp)

add_executable(CrystalMD ${MAIN_FILE} ${INCLUDE_FILES} ${ARCH_FILES} ${HEADER_FILES} ${SOURCE_FILES})

target_link_libraries(CrystalMD ${MPI_LIBRARIES}) # mpi lib
