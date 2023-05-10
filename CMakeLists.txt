cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(gmsh LANGUAGES CXX)
# set(CMAKE_CXX_STANDARD 11)

set(gmsh_src src/hangingMeshv2.cpp)

set(gmsh_inc /home/gaurav/gmsh-4.11.1-Linux64-sdk/include/gmsh.h )

add_executable(hangingMeshv2 ${gmsh_inc} ${gmsh_src})
target_include_directories(hangingMeshv2 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/include)
target_link_libraries(hangingMeshv2 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)