cmake_minimum_required(VERSION 3.12)

project(YourProjectName)

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Access the value of the environment variable
file(TO_CMAKE_PATH "$ENV{OPPIC_PATH}" OPPIC_PATH)
if(NOT OPPIC_PATH STREQUAL "")
    message("OPPIC_PATH: ${OPPIC_PATH}")
else()
    message("Environment variable OPPIC_PATH is not set or empty.")
endif()

file(TO_CMAKE_PATH "$ENV{PETSC_INSTALL_PATH}" PETSC_INSTALL_PATH)
if(NOT PETSC_INSTALL_PATH STREQUAL "")
    message("PETSC_INSTALL_PATH: ${PETSC_INSTALL_PATH}")
else()
    message("Environment variable PETSC_INSTALL_PATH is not set or empty.")
endif()

# Set compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -fPIC")

# Add optional debug flags
option(DEBUG "Enable debugging" OFF)
if(DEBUG)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O0 -g")
    add_definitions(-DDEBUG_LOG)
else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

# Add optional IEEE flags
option(IEEE "Enable IEEE flags" OFF)
if(IEEE)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fp-model strict -fp-model source -prec-div -prec-sqrt")
endif()

# Add optional PETSC flags
option(PETSC "Use PETSC" ON)
if(PETSC)
    add_definitions(-DUSE_PETSC)

#    find_package(PkgConfig REQUIRED)
#    pkg_check_modules(PETSc REQUIRED PETSc)
#
#    include_directories(${PETSc_INCLUDE_DIRS})
#    link_directories(${PETSc_LIBRARY_DIRS})
#    add_definitions(${PETSc_CFLAGS_OTHER})
endif()

# Add optional PARMETIS flags
option(PARMETIS "Use PARMETIS" OFF)
if(PARMETIS)
    add_definitions(-DHAVE_PARMETIS)
    find_package(Parmetis REQUIRED)
    include_directories(${PARMETIS_INCLUDE_DIRS})
    link_libraries(${PARMETIS_LIBRARIES})
endif()

include_directories(
    ${OPPIC_PATH}/include
    ${PETSC_INSTALL_PATH}/include
)

link_directories(
    ${OPPIC_PATH}/lib
    ${PETSC_INSTALL_PATH}/lib
)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR})
set(CPP_SOURCES
    ${SRC_DIR}/fempic_op.cpp
    ${SRC_DIR}/fempic_ori/meshes.cpp
    ${SRC_DIR}/fempic_ori/particles.cpp
    ${SRC_DIR}/fempic_ori/maths.cpp
    ${SRC_DIR}/FESolver.cpp
)

# Create the common library
add_library(common OBJECT ${CPP_SOURCES})

# Create the seq target
add_executable(seq ${SRC_DIR}/seq/seqkernels.cpp $<TARGET_OBJECTS:common>)
target_link_libraries(seq opp_seq petsc)

set(CMAKE_CXX_COMPILER "mpicxx")
add_definitions(-DUSE_MPI)
# Create the mpi_common library
add_library(mpi_common OBJECT ${SRC_DIR}/mpi/mpikernels.cpp)
target_compile_definitions(mpi_common PRIVATE -DUSE_MPI)

# Create the mpi target
add_executable(mpi ${SRC_DIR}/mpi/mpikernels.cpp $<TARGET_OBJECTS:common> $<TARGET_OBJECTS:mpi_common>)
target_link_libraries(mpi opp_mpi)
