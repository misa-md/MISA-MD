# mpi compiler toolchain for sunway architecture of taihulight supercomputer.

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR alpha)

set(TOOLCHAIN_PATH  "/usr/sw-mpp/mpi2/mpiswgcc/bin")
#set(CMAKE_C_COMPILER "${TOOLCHAIN_PATH}/mpiswgcc")
set(CMAKE_C_COMPILER "${CMAKE_CURRENT_LIST_DIR}/mpisw5gcc_wrapper.sh") # mpiswgcc wrapper: remove -mhost when compiling slave code.
set(CMAKE_CXX_COMPILER "${TOOLCHAIN_PATH}/mpiswg++")

# spesific cc and cxx link flag
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -mhybrid" CACHE STRING "" FORCE)
set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -mhybrid" CACHE STRING "" FORCE)
