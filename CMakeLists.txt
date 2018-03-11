cmake_minimum_required(VERSION 3.6)
project(CrystalMD)

# change to mpicc and mpicxx
#set(CMAKE_C_COMPILER mpicc -cc=gcc -cxx=g++)
#set(CMAKE_CXX_COMPILER mpicxx -cc=gcc -cxx=g++)

set(CMAKE_CXX_STANDARD 11)

option(OpenMP_ENABLE_FLAG "Use OpenMP" OFF) #change this flag to OFF to disable OpenMP
option(MPI_ENABLE_FLAG "Use MPI library" ON) #change this flag to false to disable mpi
option(TEST_ENABLE_FLAG "Enable test" ON) # enable test
option(TEST_MPI_ENABLE_FLAG "Enable MPI in test" ON) # enable mpi in test, its value depends on option MPI_ENABLE_FLAG.

## architecture ralated values.
option(ARCH_SW "Enable sunway athread" OFF) # enable sunway athread if its running on sunway system.

################################
# MPI and OpenMP
################################
if (OpenMP_ENABLE_FLAG)
    find_package(OpenMP REQUIRED)

    if (OPENMP_FOUND)
        set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
        set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    endif ()
endif ()

if (MPI_ENABLE_FLAG)
    find_package(MPI REQUIRED)
    MESSAGE(STATUS "MPI_INCLUDE dir:" ${MPI_INCLUDE_PATH})
    MESSAGE(STATUS "MPI_LIBRARIES dir:" ${MPI_LIBRARIES})

    if (MPI_COMPILE_FLAGS)
        set(COMPILE_FLAGS "${COMPILE_FLAGS} ${MPI_COMPILE_FLAGS}") # todo
    endif ()

    if (MPI_LINK_FLAGS)
        set(LINK_FLAGS "${LINK_FLAGS} ${MPI_LINK_FLAGS}")
    endif ()

    include_directories(${MPI_INCLUDE_PATH})

    set(EXTRA_LIBS ${EXTRA_LIBS} ${MPI_LIBRARIES}) #add mpi lib
endif ()
##### mpi and openmp end

add_subdirectory(src)

##############################
## test configure
##############################
if (TEST_ENABLE_FLAG)
    MESSAGE(STATUS "TEST is enabled")
    add_subdirectory(test)
endif ()


#############################
## install files
#############################
install(FILES README.md DESTINATION share/docs/)
install(DIRECTORY docs/ DESTINATION share/docs/d)
