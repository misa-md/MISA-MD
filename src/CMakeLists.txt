# vendor include directory
include_directories(${VENDOR_PATH}/include)

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
        hardware_accelerate.hpp
        input.h
        newton_motion.h
        simulation.h
        world_builder.h
        atom/atom_element.h
        atom/atom_list.h
        atom/inter_atom_list.h
        domain/domain.h
        domain/region.hpp
        types/atom_types.h
        types/pre_config.h
        types/pre_define.h
        types/atom_info_dump.h
        pack/lat_particle_data.h
        pack/pack.h
        pack/particledata.h
        utils/mpi_data_types.h
        utils/mpi_domain.h
        utils/rpcc.hpp
        )

set(SOURCE_FILES
        atom.cpp
        atom_dump.cpp
        crystal_md.cpp
        toml_config.cpp
        config_values.cpp
        input.cpp
        newton_motion.cpp
        simulation.cpp
        world_builder.cpp
        atom/atom_element.cpp
        atom/atom_list.cpp
        atom/inter_atom_list.cpp
        atom/inter_atom_list_minor.cpp
        domain/domain.cpp
        types/atom_info_dump.cpp
        pack/lat_particle_data.cpp
        pack/pack.cpp
        pack/particledata.cpp
        utils/mpi_domain.cpp
        utils/mpi_data_types.cpp
        )


# set binary path and lib storage path for compiling process.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# make crystal-md lib
add_library(${MD_LIB_NAME} STATIC ${HEADER_FILES} ${SOURCE_FILES})
target_include_directories(
        ${MD_LIB_NAME}
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(
        ${MD_LIB_NAME}
        PUBLIC pot::pot
        # todo kiwi, fmt.
)

add_executable(${EXECUTE_BIN_NAME} main.cpp)

target_link_libraries(${EXECUTE_BIN_NAME}
        PRIVATE ${MD_LIB_NAME} ${EXTRA_LIBS}
        )


##########################
##### install files
##########################
# install bin and libs
install(TARGETS ${EXECUTE_BIN_NAME} ${MD_LIB_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )

# install header.
install(DIRECTORY "./"
        DESTINATION "include"
        FILES_MATCHING
        PATTERN "*.h"
        PATTERN "*.hpp"
        )