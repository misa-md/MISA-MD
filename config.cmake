# all variables here start with "MD_"
set(MD_VERSION "0.2.0")

#############
## options ##
#############
# change to mpicc and mpicxx
#set(CMAKE_C_COMPILER mpicc -cc=gcc -cxx=g++)
#set(CMAKE_CXX_COMPILER mpicxx -cc=gcc -cxx=g++)

option(OpenMP_ENABLE_FLAG "Use OpenMP" OFF) #change this flag to OFF to disable OpenMP
option(MPI_ENABLE_FLAG "Use MPI library" ON) #change this flag to false to disable mpi
option(TEST_ENABLE_FLAG "Enable test" ON) # enable test
option(TEST_MPI_ENABLE_FLAG "Enable MPI in test" ON) # enable mpi in test, its value depends on option MPI_ENABLE_FLAG.
option(TOOLS_BUILD_ENABLE_FLAG "Enable tools building" ON) # enable tools building (in tools directory) binary.(tools example: convert simulation result binary file to text file)

## architecture ralated values.
option(ARCH_SW "Enable sunway athread" OFF) # enable sunway athread if its running on sunway system.

## debug
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")

#############
## const ##
#############
set(EXECUTE_BIN_NAME CrystalMD)
set(MD_LIB_NAME md) # use PARENT_SCOPE to modify globle variable.

# test
set(MD_TEST_NAME "md-test")
