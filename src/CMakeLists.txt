# include directory
include_directories(${PROJECT_SOURCE_DIR}/vendor/pkg/include)

# include kiwi framework.
add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/src/kiwi ${PROJECT_BINARY_DIR}/vendor/kiwi)
include_directories(${KIWI_EXPORT_INCLUDE_PATH})
set(EXTRA_LIBS fmt kiwi ${EXTRA_LIBS})

# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${PROJECT_SOURCE_DIR}/src/pre_config.h.in"
        "${PROJECT_SOURCE_DIR}/src/pre_config.h"
)

#set(INCLUDE_FILES
#        args.hpp
#        toml.hpp)

set(HEADER_FILES
        arch_env.hpp
        atom.h
        create_atom.h
        crystal_md.h
        toml_config.h
        config_values.h
        domain.h
        eam.h
        hardware_accelerate.hpp
        input.h
        integrator.h
        interpolation_object.h
        lat_particle_data.h
        particledata.h
        pre_config.h
        simulation.h
        )

set(SOURCE_FILES
        atom.cpp
        create_atom.cpp
        crystal_md.cpp
        toml_config.cpp
        config_values.cpp
        domain.cpp
        eam.cpp
        input.cpp
        integrator.cpp
        interpolation_object.cpp
        lat_particle_data.cpp
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

if (ARCH_SW)
    set(ARCH_FILES ${ARCH_FILES} ${ARCH_SUNWAY_FILES}) # add sunway architecture files to cmake.
endif ()

set(MAIN_FILE main.cpp)

# set binary path and lib storage path for compiling process.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# make crystal-md lib
add_library(${MD_LIB_NAME} STATIC ${ARCH_FILES} ${HEADER_FILES} ${SOURCE_FILES})
set(EXTRA_LIBS ${MD_LIB_NAME} ${EXTRA_LIBS})

add_executable(${EXECUTE_BIN_NAME} ${MAIN_FILE})

target_link_libraries(${EXECUTE_BIN_NAME} ${EXTRA_LIBS})


##########################
##### install files
##########################
install(TARGETS ${EXECUTE_BIN_NAME} ${MD_LIB_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

install(FILES ${HEADER_FILES} DESTINATION include) # todo head file
