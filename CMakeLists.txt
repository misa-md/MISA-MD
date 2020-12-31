cmake_minimum_required(VERSION 3.6)
project(MISA-MD)

set(CMAKE_CXX_STANDARD 11)

include(version.cmake)
include(config.cmake)
include(dependency.cmake)

add_subdirectory(src)
add_subdirectory(frontend)

##############################
## test configure
##############################
if (MD_TEST_ENABLE_FLAG)
    enable_testing()
    MESSAGE(STATUS "TEST is enabled")
    add_subdirectory(tests)
endif ()

##############################
## tools
##############################
if (MD_TOOLS_BUILD_ENABLE_FLAG)
    MESSAGE(STATUS "Tools Building is enabled")
    add_subdirectory(tools)
endif ()


#############################
## install files
#############################
install(FILES README.md DESTINATION share/docs/)
install(DIRECTORY docs/ DESTINATION share/docs/d)
