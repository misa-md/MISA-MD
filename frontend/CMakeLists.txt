# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/frontend_config.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/frontend_config.h"
)

# set binary path and lib storage path for compiled files.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_executable(
        ${EXECUTE_BIN_NAME}
        atom_dump.h
        atom_dump.cpp
        atom_info_dump.h
        atom_info_dump.cpp
        buffered_io.h
        buffered_io.cpp
        crystal_md.h
        crystal_md.cpp
        toml_config.h
        toml_config.cpp
        main.cpp
        device.h
        md_simulation.cpp
        md_simulation.h
        output_interface.h
        output_interface.cpp
)

target_link_libraries(${EXECUTE_BIN_NAME}
        PRIVATE ${MD_LIB_NAME} ${EXTRA_LIBS}
        # todo dependency kiwi
        )

# install bin and libs
install(TARGETS ${EXECUTE_BIN_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
