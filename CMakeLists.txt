cmake_minimum_required(VERSION 3.6)
project(CrystalMD)

# change to mpicc and mpicxx
#set(CMAKE_C_COMPILER mpicc -cc=gcc -cxx=g++)
#set(CMAKE_CXX_COMPILER mpicxx -cc=gcc -cxx=g++)

set(CMAKE_CXX_STANDARD 11)

################################
# MPI and OpenMP
################################

#find_package(OpenMP REQUIRED)
#
#if(OPENMP_FOUND)
#    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
#    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
#    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
#endif()

find_package(MPI REQUIRED)
#MESSAGE(STATUS "This is MPI_INCLUDE dir" ${MPI_INCLUDE_PATH})
#MESSAGE(STATUS "This is MPI_LIBRARIES dir" ${MPI_LIBRARIES})

if (MPI_COMPILE_FLAGS)
    set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}")
endif ()

if (MPI_LINK_FLAGS)
    set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
endif ()

include_directories(${MPI_INCLUDE_PATH})

add_subdirectory(src)

##### mpi and openmp end