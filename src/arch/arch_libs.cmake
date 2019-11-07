# this file in included in src directory.

set(CURRENT_ARCH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/arch)

include(${CURRENT_ARCH_SOURCE_DIR}/arch_configure.cmake)

# ARCH FILES
set(ARCH_FILES
        ${CURRENT_ARCH_SOURCE_DIR}/arch_building_config.h
        ${CURRENT_ARCH_SOURCE_DIR}/arch_env.hpp
        ${CURRENT_ARCH_SOURCE_DIR}/hardware_accelerate.hpp
        )

# ARCH_LIBS
if (SUNWAY_ARCH_ENABLE_FLAG)
    # if it is sunway architecture.
    add_subdirectory(${CURRENT_ARCH_SOURCE_DIR}/sunway)
    # may link sunway lib to libmd
    set(ARCH_LIBS md_sunway)
endif ()
