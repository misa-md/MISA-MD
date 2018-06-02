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

##############################
# add kiwi framework globally
##############################
add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/src/kiwi ${PROJECT_BINARY_DIR}/vendor/kiwi)
set(EXTRA_LIBS ${KIWI_EXPORT_LINK_LIBS} ${EXTRA_LIBS})

##########################################
# add googletest framework for test only
##########################################
if (TEST_ENABLE_FLAG)
    # or using ExternalProject_Add, see:https://github.com/kaizouman/gtest-cmake-example

    # download googletest from https://github.com/google/googletest/,and copy it to "external" directory.
    # Include the gtest library. gtest_SOURCE_DIR is available due to 'project(gtest)' above.
    # an issue for mingw on winodws, see: https://github.com/google/googletest/issues/606#issuecomment-234733757
    if (NOT TARGET gtest)
        add_subdirectory(${PROJECT_SOURCE_DIR}/vendor/src/googletest/googletest ${PROJECT_BINARY_DIR}/vendor/googletest)
    endif ()

    #googlemock
    #add_subdirectory("${source_dir}/googlemock")
    #include_directories(${gmock_SOURCE_DIR}/include ${gmock_SOURCE_DIR})

    # Standard linking to gtest stuff.
    set(EXTRA_LIBS ${EXTRA_LIBS} gtest gtest_main)
    MESSAGE(STATUS "Set googletest as test library.")
endif ()