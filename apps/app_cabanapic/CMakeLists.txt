cmake_minimum_required(VERSION 3.10)

project(Cabana_OPP_Project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG_LOG")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

option(DEBUG_LOG "Enable Debug Logging" OFF)
if (DEBUG_LOG)
    set(CPPFLAGS "${CPPFLAGS} -DDEBUG_LOG")
endif()

include_directories($ENV{OPP_PATH}/include)
link_directories($ENV{OPP_PATH}/lib)

find_package(MPI)
if (MPI_FOUND)
    include_directories(${MPI_INCLUDE_PATH})
endif()

find_package(CUDA)
if(CUDA_FOUND)
    enable_language(CUDA)
    set(CMAKE_CUDA_ARCHITECTURES $ENV{CUDA_ARCH})
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda -DUSE_THRUST -m64 -Xptxas -dlcm=ca")
endif()

# HIP support
option(HIP "Enable HIP support" OFF)
if (HIP)
    set(HIPCCFLAGS "-std=c++17 -Wall -D__HIP_PLATFORM_AMD__ -march=native")
endif()

# Parmetis support
option(PARMETIS "Enable Parmetis support" OFF)
if (PARMETIS)
    set(ALL_LIBS "${ALL_LIBS} -lparmetis -lmetis")
    set(ALL_INC "${ALL_INC} -DHAVE_PARMETIS")
endif()



# ----------------------------------------------------------------------------------------------------------------------
add_executable(seq
    cabana_opp.cpp
    seq/opp_kernels.cpp
)
target_compile_options(seq PRIVATE ${CPPFLAGS})
target_link_libraries(seq PRIVATE opp_seq)
message(STATUS "Configured CabanPIC seq backend")

# ----------------------------------------------------------------------------------------------------------------------
find_package(OpenMP)
if (OpenMP_CXX_FOUND)
    add_executable(omp
        cabana_opp.cpp
        omp/opp_kernels.cpp
    )
    target_compile_options(omp PRIVATE ${CPPFLAGS} ${OpenMP_CXX_FLAGS} -DUSE_OMP)
    target_link_libraries(omp PRIVATE opp_omp ${OpenMP_CXX_LIBRARIES})
    message(STATUS "Configured CabanPIC omp backend")
endif()

# ----------------------------------------------------------------------------------------------------------------------
if (MPI_FOUND)
    add_executable(mpi
        cabana_opp.cpp
        mpi/opp_kernels.cpp
    )
    target_compile_options(mpi PRIVATE ${CPPFLAGS} ${MPI_CXX_COMPILE_FLAGS} -DUSE_MPI)
    target_link_libraries(mpi PRIVATE opp_mpi ${MPI_CXX_LIBRARIES})
    message(STATUS "Configured CabanPIC mpi backend")
endif()

# ----------------------------------------------------------------------------------------------------------------------
if (CUDA_FOUND)
    add_executable(cuda
        cabana_opp.cpp
        cuda/opp_kernels.cu
    )
    set_target_properties(cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(cuda PRIVATE ${CPPFLAGS} ${CUDA_NVCC_FLAGS})
    target_link_libraries(cuda PRIVATE opp_cuda ${CUDA_LIB})
    message(STATUS "Configured CabanPIC cuda backend")
endif()

# ----------------------------------------------------------------------------------------------------------------------
if (MPI_FOUND AND CUDA_FOUND)
    add_executable(cuda_mpi
        cabana_opp.cpp
        cuda/opp_kernels.cu
    )
    set_target_properties(cuda_mpi PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(cuda_mpi PRIVATE ${CPPFLAGS} ${CUDA_NVCC_FLAGS} -DUSE_MPI)
    target_link_libraries(cuda_mpi PRIVATE opp_cuda_mpi ${CUDA_LIB})
    message(STATUS "Configured CabanPIC cuda_mpi backend")
endif()

# # ----------------------------------------------------------------------------------------------------------------------
# # HIP build
# if (HIP)
#     add_executable(hip
#         src/cabana_opp.cpp
#         src/hip/opp_kernels.cpp
#     )
#     target_compile_options(hip PRIVATE ${HIPCCFLAGS})
#     target_link_libraries(hip PRIVATE opp_hip)
# endif()

# # ----------------------------------------------------------------------------------------------------------------------
# # Build options
# add_custom_target(build_all
#     DEPENDS seq omp mpi cuda hip
# )

