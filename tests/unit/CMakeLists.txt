
# configure a header file to pass some of the CMake settings to the source code
set(TEST_TEMP_FILES_DIR ${CMAKE_CURRENT_BINARY_DIR})
configure_file(
        "${CMAKE_CURRENT_SOURCE_DIR}/test_config.h.in"
        "${CMAKE_CURRENT_SOURCE_DIR}/test_config.h"
)

#########################
# Unit Tests Settings
#########################
set(TEST_SOURCE_FILES
        test_config.h
        domain_test_utils.h
        main.cpp
        f_test.cpp
        atom_iterator.cpp
        inter_atom_test.cpp
        lattice/ws_utils_test.cpp
        world_builder_test.cpp
        domain_test_utils.cpp
        neighbour_index_test.cpp
        neighbour_index_iterator_test.cpp
        particledata_test.cpp
        )

# set binary path and lib storage path wile compiling process.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

add_executable(${MD_TEST_NAME} ${TEST_SOURCE_FILES})

# add crystldmd lib and gtest, mpi lib.
target_link_libraries(${MD_TEST_NAME}
        PRIVATE GTest::gtest GTest::gtest_main ${MD_LIB_NAME} ${MD_FRONTEND_LIB_NAME} ${EXTRA_LIBS}
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