# mpi compiler toolchain for sunway architecture of taihulight supercomputer.

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR alpha)

set(TOOLCHAIN_PATH  "/usr/sw-mpp/mpi2/mpiswgcc/bin")
set(CMAKE_C_COMPILER "${TOOLCHAIN_PATH}/mpiswgcc")
set(CMAKE_CXX_COMPILER "${TOOLCHAIN_PATH}/mpiswg++")
