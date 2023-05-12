cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(gmsh LANGUAGES CXX)
# set(CMAKE_CXX_STANDARD 11)

set(gmsh_src_v1 src/hangingMeshv1.cpp src/gmshUtils.cpp )
set(gmsh_src_v2 src/hangingMeshv2.cpp src/gmshUtils.cpp )
set(gmsh_src_v3 src/hangingMeshv3.cpp src/gmshUtils.cpp )
set(gmsh_src_v4 src/hangingMeshv4.cpp src/gmshUtils.cpp )

set(gmsh_inc /home/gaurav/gmsh-4.11.1-Linux64-sdk/include/gmsh.h include/gmshUtils.hpp )

include_directories( include )

add_executable(hangingMeshv4 ${gmsh_inc} ${gmsh_src_v4})
add_executable(hangingMeshv3 ${gmsh_inc} ${gmsh_src_v3})
add_executable(hangingMeshv2 ${gmsh_inc} ${gmsh_src_v2})
add_executable(hangingMeshv1 ${gmsh_inc} ${gmsh_src_v1})

target_include_directories(hangingMeshv4 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv4 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv3 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv3 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv2 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv2 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv1 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv1 PUBLIC /home/gaurav/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)