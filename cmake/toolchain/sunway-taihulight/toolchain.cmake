# cc/cxx toolchain for sunway architecture of taihulight supercomputer.

set(CMAKE_SYSTEM_NAME Linux)
set(CMAKE_SYSTEM_PROCESSOR alpha)


set(CMAKE_C_COMPILER "sw5gcc")
set(CMAKE_CXX_COMPILER "sw5g++")

# see also https://stackoverflow.com/questions/11423313/cmake-cross-compiling-c-flags-from-toolchain-file-ignored
#set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -mhost -std=c++11" CACHE STRING "" FORCE)
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -mhost -std=c++11" CACHE STRING "" FORCE)
add_compile_options(-mhost -std=c++11)

# spesific cc and cxx link flag
set(CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -mhybrid -std=c++11" CACHE STRING "" FORCE)
set(CMAKE_C_LINK_FLAGS "${CMAKE_C_LINK_FLAGS} -mhybrid -std=c++11" CACHE STRING "" FORCE)
