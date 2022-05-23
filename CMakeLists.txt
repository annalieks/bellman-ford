cmake_minimum_required(VERSION 3.22)
project(bellman_ford)

set(CMAKE_CXX_STANDARD 14)

find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set (CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

include_directories(headers)

add_executable(open-mp open-mp.cpp)
add_executable(mpi mpi.cpp)
add_executable(sequential sequential.cpp)

find_package(MPI)
if (MPI_FOUND)
    include_directories(SYSTEM ${MPI_INCLUDE_PATH})
    target_link_libraries(mpi ${MPI_LIBRARIES})
endif()
