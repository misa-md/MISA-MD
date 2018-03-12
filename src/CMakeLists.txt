# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${PROJECT_SOURCE_DIR}/src/pre_config.h.in"
        "${PROJECT_SOURCE_DIR}/src/pre_config.h"
)

include_directories(${PROJECT_SOURCE_DIR}/vendor/include)

#set(INCLUDE_FILES
#        args.hpp
#        toml.hpp)

set(HEADER_FILES
        arch_env.hpp
        atom.h
        createatom.h
        crystal_md.h
        config.h
        config_values.h
        domain.h
        domaindecomposition.h
        eam.h
        hardware_accelerate.hpp
        input.h
        integrator.h
        InterpolationObject.h
        latparticledata.h
        particledata.h
        pre_config.h
        simulation.h
        utils/data_pack.h
        utils/mpi_utils.h)

set(SOURCE_FILES
        atom.cpp
        createatom.cpp
        crystal_md.cpp
        config.cpp
        config_values.cpp
        domain.cpp
        domaindecomposition.cpp
        eam.cpp
        input.cpp
        integrator.cpp
        InterpolationObject.cpp
        latparticledata.cpp
        particledata.cpp
        simulation.cpp
        utils/data_pack.cpp
        utils/mpi_utils.cpp)

# codes runs on sunway.
set(ARCH_SUNWAY_FILES
        arch/sunway/athread_accelerate.h
        arch/sunway/athread_accelerate.cpp
        arch/sunway/sunway_env.h
        arch/sunway/sunway_env.cpp
        arch/sunway/sunway_slave.c
        arch/sunway/sunway_slave2.c
        )

if (ARCH_SW)
    set(ARCH_FILES ${ARCH_FILES} ${ARCH_SUNWAY_FILES}) # add sunway architecture files to cmake.
endif ()

set(MAIN_FILE main.cpp)

# set binary path and lib storage path wile compiling process.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# make crystal-md lib
set(MD_LIB_NAME crystalmd)
add_library(${MD_LIB_NAME} STATIC ${ARCH_FILES} ${HEADER_FILES} ${SOURCE_FILES})
set(EXTRA_LIBS ${EXTRA_LIBS} ${MD_LIB_NAME})

set(PROGRAM_NAME ${PROJECT_NAME})
add_executable(${PROGRAM_NAME} ${MAIN_FILE})

target_link_libraries(${PROGRAM_NAME} ${EXTRA_LIBS})


##########################
##### install files
##########################
install(TARGETS ${PROGRAM_NAME} ${MD_LIB_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

install(FILES ${HEADER_FILES} DESTINATION include) #todo head file
