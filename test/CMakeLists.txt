enable_testing()

# Extra linking for the project.
set(EXTRA_LIBS gtest gtest_main ${MD_LIB_NAME} ${EXTRA_LIBS})  # add crystal-md lib and gtest lib
MESSAGE(STATUS "Set googletest as test library.")

if (NOT (MPI_ENABLE_FLAG))
    set(TEST_MPI_ENABLE_FLAG OFF)  # if we disabled MPI_ENABLE_FLAG,we cann't also use TEST_MPI_ENABLE_FLAG.
endif ()
MESSAGE(STATUS "MPI in test: ${TEST_MPI_ENABLE_FLAG}")

# configure a header file to pass some of the CMake settings to the source code
configure_file(
        "${PROJECT_SOURCE_DIR}/test/test_config.h.in"
        "${PROJECT_SOURCE_DIR}/test/test_config.h"
)

#########################
# Unit Tests Settings
#########################
set(TEST_SOURCE_FILES
        test_config.h
        domain_test_utils.h
        main.cpp
        f_test.cpp
        domain_test.cpp
        atom_iterator.cpp
        world_builder_test.cpp
        domain_test_utils.cpp
        eam_test.cpp
        )

# set binary path and lib storage path wile compiling process.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_executable(${MD_TEST_NAME} ${TEST_SOURCE_FILES})

# add crystldmd lib and gtest, mpi lib.
target_link_libraries(${MD_TEST_NAME}
        PRIVATE ${EXTRA_LIBS}
        )

# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable runUnitTests to see those specific tests.
add_test(NAME ${MD_TEST_NAME} COMMAND ${MD_TEST_NAME})
#if (MPI_ENABLE_FLAG)
#    set(MPI_PROC_NO 4)
#    add_test(NAME crystalMD_Tests COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} ${MPI_PROC_NO} ${MD_TEST_NAME}) # run with 4 processors.
#else ()
#    add_test(NAME crystalMD_Tests COMMAND ${MD_TEST_NAME})
#endif ()


##########################
##### install files
##########################
install(TARGETS ${MD_TEST_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
