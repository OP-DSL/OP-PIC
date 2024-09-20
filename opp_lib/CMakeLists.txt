cmake_minimum_required(VERSION 3.10)

project(OPP_Project LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG_LOG")

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE "Release")
endif()

find_package(OpenMP)
find_package(CUDA)
find_package(MPI)

if(CUDA_FOUND)
    option(CUDA "Enable CUDA support" ON)
    set(CMAKE_CUDA_ARCHITECTURES $ENV{CUDA_ARCH})
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --extended-lambda -DUSE_THRUST -Xptxas -dlcm=ca")
endif()

include_directories($ENV{OPP_PATH}/include)
link_directories($ENV{OPP_PATH}/lib)

# ------------------------------------------------------------------------------------------
set(SRC_SEQ_CPP
    "src/core/opp_lib_core.cpp"
    "src/core/opp_util.cpp"
    "src/core/opp_params.cpp"
    "src/seq/opp_seq.cpp"
)

add_library(opp_seq STATIC ${SRC_SEQ_CPP})
# target_compile_options(opp_seq PRIVATE ${CPPFLAGS})
# target_include_directories(opp_seq PRIVATE "$ENV{OPP_PATH}/include")

add_custom_target(seq ALL DEPENDS opp_seq)
message(STATUS "Configured OP-PIC seq backend")

# ------------------------------------------------------------------------------------------
if (OpenMP_CXX_FOUND)
    set(SRC_OMP_CPP 
        "src/core/opp_lib_core.cpp"
        "src/core/opp_util.cpp"
        "src/core/opp_params.cpp"
        "src/omp/opp_omp.cpp"
    )
    add_library(opp_omp STATIC ${SRC_OMP_CPP})
        # target_compile_options(opp_omp PRIVATE ${CPPFLAGS} ${OpenMP_CXX_FLAGS})
    target_compile_options(opp_omp PRIVATE ${OpenMP_CXX_FLAGS})
    # target_link_libraries(opp_omp PRIVATE OpenMP::OpenMP_CXX)
        # target_include_directories(opp_omp PRIVATE "$ENV{OPP_PATH}/include")

    add_custom_target(omp ALL DEPENDS opp_omp)
    message(STATUS "Configured OP-PIC omp backend")
endif()

# ------------------------------------------------------------------------------------------
if(MPI_FOUND)
    set(SRC_MPI_CPP 
        "src/core/opp_lib_core.cpp"
		"src/core/opp_util.cpp"
		"src/core/opp_params.cpp"
		"src/mpi/opp_mpi_halo.cpp"
		"src/mpi/opp_mpi_halo_core.cpp"
		"src/mpi/opp_mpi_partition.cpp"
		"src/mpi/opp_mpi.cpp"
		"src/mpi/opp_mpi_core.cpp"
		"src/mpi/opp_mpi_utils.cpp"
		"src/mpi/opp_mpi_particle_comm.cpp"
		"src/mpi/opp_mpi_double_ind_reducs.cpp"
    )
    add_library(opp_mpi STATIC ${SRC_MPI_CPP})
    target_compile_options(opp_mpi PRIVATE ${CPPFLAGS} -DUSE_MPI ${MPI_CXX_COMPILE_FLAGS})
    # set_source_files_properties(${SRC_MPI_CPP} PROPERTIES COMPILE_FLAGS "${MPI_CXX_COMPILE_FLAGS}")
    # target_link_libraries(opp_mpi PRIVATE MPI::MPI_CXX)
    target_include_directories(opp_mpi PRIVATE ${MPI_INCLUDE_PATH})

    add_custom_target(mpi ALL DEPENDS opp_mpi)
    message(STATUS "Configured OP-PIC mpi backend")
endif()

# ------------------------------------------------------------------------------------------
if (CUDA_FOUND)
    enable_language(CUDA)

    set(SRC_CUDA_CPP
        "src/core/opp_lib_core.cpp"
        "src/core/opp_util.cpp"
        "src/core/opp_params.cpp"
        "src/cuda/opp_cuda.cu"
        "src/cuda/opp_particle_sorter.cu"
        "src/cuda/opp_cuda_halo.cu"
    )

    add_library(opp_cuda STATIC ${SRC_CUDA_CPP})
    target_compile_options(opp_cuda PRIVATE ${CPPFLAGS} -DUSE_THRUST)
    set_target_properties(opp_cuda PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON
    )
    target_include_directories(opp_cuda PRIVATE 
        $ENV{OPP_PATH}/include 
        $ENV{CUDA_INSTALL_PATH}/include
    )

    add_custom_target(cuda ALL DEPENDS opp_cuda)
    message(STATUS "Configured OP-PIC cuda backend")
endif()

if (MPI_FOUND AND CUDA_FOUND)
    set(SRC_CUDA_MPI_CPP
        "src/core/opp_lib_core.cpp"
        "src/core/opp_util.cpp"
        "src/core/opp_params.cpp"
        "src/mpi/opp_mpi_halo_core.cpp"
        "src/mpi/opp_mpi_partition.cpp"
        "src/mpi/opp_mpi_core.cpp"
        "src/mpi/opp_mpi_particle_comm.cpp"
        "src/mpi/opp_mpi_double_ind_reducs.cpp"
        "src/mpi/opp_mpi_utils.cpp"
        "src/cuda/opp_cuda.cu"
        "src/cuda/opp_particle_sorter.cu"
        "src/cuda/opp_cuda_halo.cu"
    )

    add_library(opp_cuda_mpi STATIC ${SRC_CUDA_MPI_CPP})
    target_compile_options(opp_cuda_mpi PRIVATE ${CPPFLAGS} -DUSE_THRUST -DUSE_MPI)
    set_target_properties(opp_cuda_mpi PROPERTIES
        CUDA_SEPARABLE_COMPILATION ON
        POSITION_INDEPENDENT_CODE ON
    )
    target_include_directories(opp_cuda_mpi PRIVATE 
        $ENV{OPP_PATH}/include 
        $ENV{CUDA_INSTALL_PATH}/include 
        ${MPI_INCLUDE_PATH}
    )

    add_custom_target(cuda_mpi ALL DEPENDS opp_cuda_mpi)
    message(STATUS "Configured OP-PIC cuda_mpi backend")
endif()

# ------------------------------------------------------------------------------------------
