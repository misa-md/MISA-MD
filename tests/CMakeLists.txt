if (NOT (MD_MPI_ENABLE_FLAG))
    set(MD_TEST_MPI_ENABLE_FLAG OFF)  # if we disabled MPI_ENABLE_FLAG,we cann't also use MD_TEST_MPI_ENABLE_FLAG.
endif ()
MESSAGE(STATUS "MPI in test: ${MD_TEST_MPI_ENABLE_FLAG}")

add_subdirectory(unit)
