set(SUNWAY_ARCH_FILES
        athread_accelerate.h
        athread_accelerate.cpp
        sunway_env.h
        sunway_env.cpp
        sunway_slave.c
        # sunway_slave2.c
        )

set_source_files_properties(sunway_slave.c
        PROPERTIES COMPILE_FLAGS
        "-mslave"
        )

set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

# make md-sunway lib
set(SUNWAY_LIB_NAME md_sunway)
add_library(${SUNWAY_LIB_NAME} STATIC ${SUNWAY_ARCH_FILES})
target_include_directories(
        ${SUNWAY_LIB_NAME}
        PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}
)
target_link_libraries(
        ${SUNWAY_LIB_NAME}
        PUBLIC
        m_slave # math lib for sunway slave cores.
)
