# set arch relative macor and write to a header file 'arch_building_config.h'.

if (SUNWAY_ARCH_ENABLE_FLAG)
    set(ARCH_SUNWAY ON)
endif ()

configure_file(
        "${CURRENT_ARCH_SOURCE_DIR}/arch_building_config.h.in"
        "${CURRENT_ARCH_SOURCE_DIR}/arch_building_config.h"
)
