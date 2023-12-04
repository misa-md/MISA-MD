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
    MESSAGE(SEND_ERROR "unsupported random number generation method ${MD_RAND}")
endif ()

if (MD_RUNTIME_CHECKING_FLAG)
    set(MD_RUNTIME_CHECKING ON)
endif ()

# set kernel strategy.
string(TOLOWER ${MD_ATOMS_MEMORY_LAYOUT} MD_ATOMS_MEMORY_LAYOUT_LOWER)
if (MD_ATOMS_MEMORY_LAYOUT_LOWER MATCHES "default")
    set(MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS ON) # default memory layout is AoS
elseif (MD_ATOMS_MEMORY_LAYOUT_LOWER MATCHES "soa")
    set(MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_SOA ON)
elseif (MD_ATOMS_MEMORY_LAYOUT_LOWER MATCHES "aos")
    set(MD_ATOM_HASH_ARRAY_MEMORY_LAYOUT_AOS ON)
else ()
    MESSAGE(FATAL_ERROR "unsupported kernel strategy ${MD_ATOMS_MEMORY_LAYOUT}")
endif ()
MESSAGE(STATUS "current kernel strategy is: ${MD_ATOMS_MEMORY_LAYOUT}")

configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/md_building_config.h.in"
        "${CONFIGURE_GENERATED_PATH}/md_building_config.h"
)

install(FILES "${CONFIGURE_GENERATED_PATH}/md_building_config.h"
        DESTINATION "include"
        )
