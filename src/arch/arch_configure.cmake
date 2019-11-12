# set arch relative macor and write to a header file 'arch_building_config.h'.

# set arch flags, and source files dir.
set(ACCELERATE_ENABLE OFF)
if (SUNWAY_ARCH_ENABLE_FLAG)
    set(ACCELERATE_ENABLE ON)
    set(ARCH_SUNWAY ON)
elseif (CUDA_ARCH_ENABLE_FLAG)
    set(ACCELERATE_ENABLE ON)
    set(ARCH_CUDA ON)
    set(ARCH_HIP ON)
endif ()

configure_file(
        "${CURRENT_ARCH_SOURCE_DIR}/arch_building_config.h.in"
        "${CURRENT_ARCH_SOURCE_DIR}/arch_building_config.h"
)
