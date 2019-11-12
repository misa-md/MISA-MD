# this file in included in src directory.

set(CURRENT_ARCH_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/arch)

include(${CURRENT_ARCH_SOURCE_DIR}/arch_configure.cmake)

# set ARCH_SRC_PATH and ARCH_LIBS
if (SUNWAY_ARCH_ENABLE_FLAG) #  Sunway
    set(ARCH_LIBS md_sunway)  # may link sunway lib to libmd
    if (SUNWAY_ARCH_SRC_PATH STREQUAL "")
        set(ARCH_SRC_PATH "${CURRENT_ARCH_SOURCE_DIR}/sunway")
    else ()
        set(ARCH_SRC_PATH ${SUNWAY_ARCH_SRC_PATH})
    endif ()
elseif (CUDA_ARCH_ENABLE_FLAG) # CUDA/HIP
#    set(ARCH_LIBS md_cuda)
    if (CUDA_ARCH_SRC_PATH STREQUAL "")
        set(ARCH_SRC_PATH "${CURRENT_ARCH_SOURCE_DIR}/cuda")
    else ()
        set(ARCH_SRC_PATH ${CUDA_ARCH_SRC_PATH})
    endif ()
endif ()

# check ARCH_SRC_PATH dir
if (SUNWAY_ARCH_ENABLE_FLAG OR CUDA_ARCH_ENABLE_FLAG)
    if (NOT EXISTS "${ARCH_SRC_PATH}" OR NOT IS_DIRECTORY "${ARCH_SRC_PATH}")
        message(FATAL_ERROR "Architecture source files directory not found: ${ARCH_SRC_PATH}")
    else ()
        MESSAGE(STATUS "Arch source files is ${ARCH_SRC_PATH}")
        add_subdirectory(${ARCH_SRC_PATH})
    endif ()
endif ()

# ARCH FILES
set(ARCH_FILES
        ${CURRENT_ARCH_SOURCE_DIR}/arch_building_config.h
        ${CURRENT_ARCH_SOURCE_DIR}/arch_env.hpp
        ${CURRENT_ARCH_SOURCE_DIR}/hardware_accelerate.hpp
        )
