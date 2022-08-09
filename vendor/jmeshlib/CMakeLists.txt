cmake_minimum_required(VERSION 3.1)
project(jmeshlib)

set (CMAKE_CXX_STANDARD 11)

include_directories(./include)

if (NOT MSVC)
add_compile_options(-fpermissive)
endif()

add_definitions(-DIS64BITPLATFORM)

set(HEADERS
    ./include/binTree.h
    ./include/clusterGraph.h
    ./include/dijkstraGraph.h
    ./include/edge.h
    ./include/graph.h
    ./include/heap.h
    ./include/j_mesh.h
    ./include/jmesh.h
    ./include/jqsort.h
    ./include/list.h
    ./include/matrix.h
    ./include/point.h
    ./include/tin.h
    ./include/triangle.h
    ./include/vertex.h
)

set(SOURCES
    ./src/MESH_STRUCTURE/io.cpp
    ./src/MESH_STRUCTURE/point.cpp
    ./src/MESH_STRUCTURE/tin.cpp
    ./src/MESH_STRUCTURE/triangle.cpp
    ./src/MESH_STRUCTURE/vertex.cpp
    ./src/MESH_STRUCTURE/edge.cpp
    ./src/MESH_STRUCTURE/checkAndRepair.cpp
    ./src/JMESH/jmesh.cpp
    ./src/PRIMITIVES/dijkstraGraph.cpp
    ./src/PRIMITIVES/graph.cpp
    ./src/PRIMITIVES/list.cpp
    ./src/PRIMITIVES/matrix.cpp
    ./src/PRIMITIVES/binTree.cpp
    ./src/PRIMITIVES/clusterGraph.cpp
    ./src/PRIMITIVES/jqsort.cpp
    ./src/PRIMITIVES/heap.cpp
)

add_library(jmeshlib STATIC ${SOURCES} ${HEADERS})