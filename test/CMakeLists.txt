enable_testing()

# we use catch2 as test lib.
include_directories(${PROJECT_SOURCE_DIR}/vendor/include)


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
        catch2_test_env.hpp
        main.cpp
#        app_test.cpp
        f_test.cpp)

include_directories(${PROJECT_SOURCE_DIR}/src)

# set binary path and lib storage path wile compiling process.
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)

set(MD_TEST_NAME "md-test")
add_executable(${MD_TEST_NAME} ${TEST_SOURCE_FILES})

# Extra linking for the project.
target_link_libraries(${MD_TEST_NAME} ${EXTRA_LIBS}) # including crystldmd lib and mpi lib.

# This is so you can do 'make test' to see all your tests run, instead of
# manually running the executable runUnitTests to see those specific tests.
add_test(NAME crystalMD_Tests COMMAND ${MD_TEST_NAME})


##########################
##### install files
##########################
install(TARGETS ${MD_LIB_NAME}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib
        )
