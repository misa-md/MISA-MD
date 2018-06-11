# include directory
include_directories(${PROJECT_SOURCE_DIR}/vendor/pkg/include)

# include kiwi framework.
include_directories(${KIWI_EXPORT_INCLUDE_PATH})

# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${PROJECT_SOURCE_DIR}/src/types/pre_config.h.in"
        "${PROJECT_SOURCE_DIR}/src/types/pre_config.h"
)

set(HEADER_FILES
        arch_env.hpp
        atom.h
        atom_dump.h
        crystal_md.h
        toml_config.h
        config_values.h
        def_config_values.h
        domain.h
        hardware_accelerate.hpp
        input.h
        newton_motion.h
        simulation.h
        world_builder.h
        atom/atom_element.h
        atom/atom_list.h
        atom/inter_atom_list.h
        types/atom_types.h
        types/pre_config.h
        types/pre_define.h
        types/atom_info_dump.h
        pack/lat_particle_data.h
        pack/pack.h
        pack/particledata.h
        potential/eam.h
        potential/eam_one_way.h
        potential/eam_phi.h
        potential/eam_parser.h
        potential/interpolation_object.h
        utils/rpcc.hpp
        )

set(SOURCE_FILES
        atom.cpp
        atom_dump.cpp
        crystal_md.cpp
        toml_config.cpp
        config_values.cpp
        domain.cpp
        input.cpp
        newton_motion.cpp
        simulation.cpp
        world_builder.cpp
        atom/atom_element.cpp
        atom/atom_list.cpp
        atom/inter_atom_list.cpp
        types/atom_info_dump.cpp
        pack/lat_particle_data.cpp
        pack/pack.cpp
        pack/particledata.cpp
        potential/eam.cpp
        potential/eam_one_way.cpp
        potential/eam_phi.cpp
        potential/eam_parser.cpp
        potential/interpolation_object.cpp
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
