# configure a header file to pass some of the CMake settings to the source code

if (MD_RAND MATCHES "LCG")
    set(RAND_LCG ON)
elseif (MD_RAND MATCHES "MT")
    set(RAND_MT ON)
elseif (MD_RAND MATCHES "STC")
    set(RAND_STC ON)
elseif (MD_RAND MATCHES "xoshiro")
    set(RAND_XOSHIRO ON)
elseif (MD_RAND MATCHES "LEGACY")
    set(RAND_LEGACY ON)
elseif (MD_RAND MATCHES "REAL")
    #    set(RAND_LINUX_REAL TRUE)
    MESSAGE(SEND_ERROR "real rand number is not currently supported")
else ()
    MESSAGE(SEND_ERROR "unsupported random number generation method ${KMC_RAND}")
endif ()


configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/md_building_config.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/md_building_config.h"
)
