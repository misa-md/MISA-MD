# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${PROJECT_SOURCE_DIR}/src/types/pre_config.h.in"
        "${PROJECT_SOURCE_DIR}/src/types/pre_config.h"
)

set(HEADER_FILES
        arch_env.hpp
        atom.h
        atom_dump.h
        config_values.h
        def_config_values.h
        hardware_accelerate.hpp
        input.h
        newton_motion.h
        simulation.h
        system_configuration.h
        world_builder.h
        atom/atom_element.h
        atom/atom_list.h
        atom/atom_set.h
        atom/inter_atom_list.h
        atom/neighbour_index.h
        atom/neighbour_iterator.h
        lattice/lattice.h
        lattice/ws_utils.h
        types/atom_types.h
        types/pre_config.h
        types/pre_define.h
        types/atom_info_dump.h
        pack/df_embed_packer.h
        pack/force_packer.h
        pack/inter_border_packer.h
        pack/inter_particle_packer.h
        pack/lat_particle_data.h
        pack/lat_particle_packer.h
        pack/particledata.h
        pack/rho_packer.h
        utils/mpi_data_types.h
        utils/mpi_domain.h
        utils/rpcc.hpp
        )

set(SOURCE_FILES
        atom.cpp
        atom_dump.cpp
        config_values.cpp
        input.cpp
        newton_motion.cpp
        simulation.cpp
        system_configuration.cpp
        world_builder.cpp
        atom/atom_element.cpp
        atom/atom_list.cpp
        atom/atom_set.cpp
        lattice/box.h
        atom/inter_atom_list.cpp
        atom/neighbour_index.inl
        atom/neighbour_iterator.inl
        lattice/lattice.cpp
        lattice/ws_utils.cpp
        types/atom_info_dump.cpp
        pack/df_embed_packer.cpp
        pack/force_packer.cpp
        pack/inter_border_packer.cpp
        pack/inter_particle_packer.cpp
        pack/lat_particle_data.cpp
        pack/lat_particle_packer.cpp
        pack/particledata.cpp
        pack/rho_packer.cpp
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
        PUBLIC
        pot::pot
        comm::comm
        kiwi::kiwi
        # fmt is included in kiwi.
)

##########################
##### install files
##########################
# install bin and libs
install(TARGETS ${MD_LIB_NAME}
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