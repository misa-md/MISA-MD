enable_testing()

# Extra linking for the project.
set(EXTRA_LIBS gtest gtest_main ${MD_LIB_NAME} ${EXTRA_LIBS})  # add crystal-md lib and gtest lib
MESSAGE(STATUS "Set googletest as test library.")

if (NOT (MPI_ENABLE_FLAG))
    set(TEST_MPI_ENABLE_FLAG OFF)  # if we disabled MPI_ENABLE_FLAG,we cann't also use TEST_MPI_ENABLE_FLAG.
endif ()
MESSAGE(STATUS "MPI in test: ${TEST_MPI_ENABLE_FLAG}")

add_subdirectory(unit)
