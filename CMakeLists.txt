project(dem2mesh)
cmake_minimum_required(VERSION 2.8)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif()

set (CMAKE_CXX_STANDARD 11)
find_package(GDAL REQUIRED)
find_package(OpenMP REQUIRED)
if(OPENMP_FOUND)
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()
include_directories(${GDAL_INCLUDE_DIR})

# Add compiler options.
if (NOT WIN32)
add_definitions(-Wall -Wextra)
endif()

# Add source directory
aux_source_directory("./src" SRC_LIST)

add_subdirectory("vendor/jmeshlib")
include_directories("vendor/jmeshlib/include")

# Add exectuteable
add_executable(${PROJECT_NAME} ${SRC_LIST})

target_link_libraries(${PROJECT_NAME} ${GDAL_LIBRARY} jmeshlib)

install(TARGETS ${PROJECT_NAME} RUNTIME DESTINATION bin)
