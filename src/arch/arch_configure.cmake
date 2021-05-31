# set arch relative macor and write to a header file 'arch_building_config.h'.

# set arch flags, and source files dir.
set(ACCELERATE_ENABLED OFF)
if (MD_SUNWAY_ARCH_ENABLE_FLAG)
    set(ARCH_NAME sunway)
    set(ACCELERATE_ENABLED ON)
    set(ARCH_SUNWAY ON)
elseif (MD_CUDA_ARCH_ENABLE_FLAG)
    set(ARCH_NAME cuda)
    set(ACCELERATE_ENABLED ON)
    set(ARCH_CUDA ON)
elseif (MD_HIP_ARCH_ENABLE_FLAG)
    set(ARCH_NAME hip)
    set(ACCELERATE_ENABLED ON)
    set(ARCH_HIP ON)
endif ()

configure_file(
        "${CURRENT_ARCH_SOURCE_DIR}/arch_building_config.h.in"
        "${CONFIGURE_GENERATED_PATH}/arch_building_config.h"
)

install(FILES "${CONFIGURE_GENERATED_PATH}/arch_building_config.h"
        DESTINATION "include"
        )
