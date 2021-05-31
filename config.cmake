#############
## options ##
#############
# change to mpicc and mpicxx
#set(CMAKE_C_COMPILER mpicc -cc=gcc -cxx=g++)
#set(CMAKE_CXX_COMPILER mpicxx -cc=gcc -cxx=g++)

option(MD_OpenMP_ENABLE_FLAG "Use OpenMP" OFF) #change this flag to OFF to disable OpenMP
option(MD_MPI_ENABLE_FLAG "Use MPI library" ON) #change this flag to false to disable mpi
option(MD_TEST_ENABLE_FLAG "Enable test" ON) # enable test
option(MD_TEST_MPI_ENABLE_FLAG "Enable MPI in test" ON) # enable mpi in test, its value depends on option MPI_ENABLE_FLAG.
option(MD_TOOLS_BUILD_ENABLE_FLAG "Enable tools building" ON) # enable tools building (in tools directory) binary.(tools example: convert simulation result binary file to text file)

## architecture related values (only used in src/arch dir).
option(MD_SUNWAY_ARCH_ENABLE_FLAG "Enable sunway athread" OFF) # enable sunway athread if its running on sunway system.
option(MD_CUDA_ARCH_ENABLE_FLAG "Enable GPU using CUDA" OFF) # enable GPU if its running on nvidia GPU.
option(MD_HIP_ARCH_ENABLE_FLAG "Enable DCU using HIP" OFF) # enable DCU if its running on sugon DCU.
set(MD_SUNWAY_ARCH_SRC_PATH "" CACHE PATH "Source files directory of sunway architecture") # source file directory of sunway arch code.
set(MD_CUDA_ARCH_SRC_PATH "" CACHE PATH "Source directory of CUDA architecture code") # source file directory of cuda arch code.
set(MD_HIP_ARCH_SRC_PATH "" CACHE PATH "Source directory of HIP architecture code") # source file directory of hip arch code.

set(MD_RAND "MT" CACHE STRING "random number generating algorithm") # random number generating
# options are:
#   LCG: linear congruential
#   MT: mersenne twister
#   STC: subtract with carry
#   xoshiro: http://xoshiro.di.unimi.it
#   LEGACY: lammps legacy
#   REAL: real random number privided by linux OS.

if (CMAKE_BUILD_TYPE MATCHES "^(Debug|DEBUG|debug)$")
    set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
    set(CMAKE_C_FLAGS_DEBUG "${CMAKE_C_FLAGS_DEBUG} -Wall")

    set(MD_DEV_MODE ON)
endif ()


#############
## const ##
#############
# all variables here start with "MD_"
set(EXECUTE_BIN_NAME misamd)
set(MD_LIB_NAME md) # use PARENT_SCOPE to modify globle variable.
set(MD_FRONTEND_LIB_NAME frontend)

set(CONFIGURE_GENERATED_PATH ${CMAKE_BINARY_DIR}/generated)
# test
set(MD_TEST_NAME "md-test")
