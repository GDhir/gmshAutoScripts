cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(gmsh LANGUAGES CXX)
# set(CMAKE_CXX_STANDARD 11)

set(gmsh_src_v1 src/hangingMeshv1.cpp src/gmshUtils.cpp )
set(gmsh_src_v2 src/hangingMeshv2.cpp src/gmshUtils.cpp )
set(gmsh_src_v3 src/hangingMeshv3.cpp src/gmshUtils.cpp )
set(gmsh_src_v4 src/hangingMeshv4.cpp src/gmshUtils.cpp )
set(gmsh_src_v5 src/hangingMeshv5.cpp src/gmshUtils.cpp )
set(gmsh_src_v6 src/hangingMeshv6.cpp src/gmshUtils.cpp )
set(gmsh_src_v7 src/hangingMeshv7.cpp src/gmshUtils.cpp )
set(gmsh_src_v8 src/hangingMeshv8.cpp src/gmshUtils.cpp )
set(gmsh_src_regularv1 src/regularMeshv1.cpp src/gmshUtils.cpp )
set(gmsh_src_regularv2 src/regularMeshv2.cpp src/gmshUtils.cpp )
set(gmsh_src_regularv3 src/regularMeshv3.cpp src/gmshUtils.cpp )
set(gmsh_src_import src/hangingMeshImport.cpp src/gmshUtils.cpp )
set(gmsh_src_order src/hangingMeshSetOrder.cpp src/gmshUtils.cpp )

set(gmsh_inc $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include/gmsh.h include/gmshUtils.hpp )

include_directories( include )

add_executable(hangingMeshImport ${gmsh_inc} ${gmsh_src_import})
add_executable(hangingMeshSetOrder ${gmsh_inc} ${gmsh_src_order})
add_executable(regularMeshv1 ${gmsh_inc} ${gmsh_src_regularv1})
add_executable(regularMeshv2 ${gmsh_inc} ${gmsh_src_regularv2})
add_executable(regularMeshv3 ${gmsh_inc} ${gmsh_src_regularv3})
add_executable(hangingMeshv8 ${gmsh_inc} ${gmsh_src_v8})
add_executable(hangingMeshv7 ${gmsh_inc} ${gmsh_src_v7})
add_executable(hangingMeshv6 ${gmsh_inc} ${gmsh_src_v6})
add_executable(hangingMeshv5 ${gmsh_inc} ${gmsh_src_v5})
add_executable(hangingMeshv4 ${gmsh_inc} ${gmsh_src_v4})
add_executable(hangingMeshv3 ${gmsh_inc} ${gmsh_src_v3})
add_executable(hangingMeshv2 ${gmsh_inc} ${gmsh_src_v2})
add_executable(hangingMeshv1 ${gmsh_inc} ${gmsh_src_v1})

target_include_directories(regularMeshv3 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(regularMeshv3 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(regularMeshv2 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(regularMeshv2 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(regularMeshv1 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(regularMeshv1 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshImport PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshImport PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshSetOrder PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshSetOrder PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv8 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv8 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv7 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv7 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv6 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv6 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv5 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv5 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv4 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv4 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv3 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv3 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv2 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv2 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)

target_include_directories(hangingMeshv1 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/include include/)
target_link_libraries(hangingMeshv1 PUBLIC $ENV{HOME}/gmsh-4.11.1-Linux64-sdk/lib/libgmsh.so)