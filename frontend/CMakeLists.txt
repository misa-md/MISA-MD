# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/frontend_config.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/frontend_config.h"
)

# set binary path and lib storage path for compiled files.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_executable(${EXECUTE_BIN_NAME}
        crystal_md.h
        crystal_md.cpp
        toml_config.h
        toml_config.cpp
        main.cpp
        device.h
        md_simulation.cpp
        md_simulation.h
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
