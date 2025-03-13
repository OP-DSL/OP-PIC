
SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = $(CC_COMPILER)
MPICPP	   = $(MPI_COMPILER)
NVCC       = nvcc
HIPCC      = $(HIP_COMPILER)
SYCLCC     = icpx

CPPFLAGS   = -std=c++17 -Wall -fPIC -march=native $(CPPFLAGS_ADD)
NVCCFLAGS  = -std=c++17 -m64 -Xptxas -dlcm=ca $(NVCCFLAGS_ADD) -DUSE_THRUST 
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -Wall -isystem -D__HIP_PLATFORM_AMD__ -munsafe-fp-atomics -DUSE_THRUST -march=native
SYCLFLAGS  = -w -fsycl -std=c++17 $(SYCLFLAGS_ADD) -DUSE_SYCL 

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g -lineinfo
	SYCLFLAGS  += -O0 -g
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g -lineinfo
	HIPCCFLAGS += -O3 -g
	SYCLFLAGS  += -O3 -g
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
	HIPCCFLAGS += -DDEBUG_LOG
	SYCLFLAGS  += -DDEBUG_LOG
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
	HIPCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
# for demos with gcc
# CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include -I$(MPI_INSTALL_PATH)/include
# CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda -I$(MPI_INSTALL_PATH)/lib -lmpi

HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include -D__HIP_PLATFORM_AMD__
HIP_LIB   = -L$(ROCM_INSTALL_DIR)/lib -lamdhip64
# # for archer2 gpus
# HIP_INC   = -I${CRAY_MPICH_DIR}/include
# HIP_LIB   = -L$(CRAY_MPICH_DIR)/lib -lmpi ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a} -lamdhip64

PARMETIS ?= 0
ifeq ($(PARMETIS), 1)
	ALL_LIBS += -lparmetis -lmetis
	ALL_INC  += -DHAVE_PARMETIS
endif

HDF5 ?= 0
ifeq ($(HDF5), 1)
	ALL_INC  += -I$(HDF5_INSTALL_PATH)/include 
	ALL_LIBS += -L$(HDF5_INSTALL_PATH)/lib
endif

ALL_INC  += -I$(OPP_PATH)/include
ALL_LIBS += -L$(OPP_PATH)/lib

.PHONY: clean mklib

mklib:
	@mkdir -p $(OBJ) $(BIN)

# ----------------------------------------------------------------------------------------------------------------------
cabana_opp_seq: 
		 $(CPP) $(CPPFLAGS) -c cabana_opp.cpp -o $(OBJ)/cabana_opp_seq.o $(ALL_INC)
opp_kernels_seq: 
	     $(CPP) $(CPPFLAGS) -c seq/opp_kernels.cpp -o $(OBJ)/opp_kernels_seq.o $(ALL_INC)
seq: mklib opp_kernels_seq cabana_opp_seq
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq \
			$(OBJ)/cabana_opp_seq.o \
			$(OBJ)/opp_kernels_seq.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_seq -lz

# ----------------------------------------------------------------------------------------------------------------------
cabana_opp_omp:
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c cabana_opp.cpp -o $(OBJ)/cabana_opp_omp.o $(ALL_INC)
opp_kernels_omp:
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp.o $(ALL_INC)
omp: mklib opp_kernels_omp cabana_opp_omp
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp \
			$(OBJ)/cabana_opp_omp.o \
			$(OBJ)/opp_kernels_omp.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_omp -lz

# ----------------------------------------------------------------------------------------------------------------------
cabana_opp_mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c cabana_opp.cpp -o $(OBJ)/cabana_opp_mpi.o $(ALL_INC)
opp_kernels_mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/opp_kernels.cpp -o $(OBJ)/opp_kernels_mpi.o $(ALL_INC)
mpi: mklib cabana_opp_mpi opp_kernels_mpi
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi \
		 	$(OBJ)/cabana_opp_mpi.o \
			$(OBJ)/opp_kernels_mpi.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_mpi -lz

# ----------------------------------------------------------------------------------------------------------------------
cabana_opp_omp_mpi:
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -DUSE_OMP -c cabana_opp.cpp -o $(OBJ)/cabana_opp_omp_mpi.o $(ALL_INC)
opp_kernels_omp_mpi:
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp_mpi.o $(ALL_INC)
omp_mpi: mklib cabana_opp_omp_mpi opp_kernels_omp_mpi
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -o $(BIN)/omp_mpi \
			$(OBJ)/cabana_opp_omp_mpi.o \
			$(OBJ)/opp_kernels_omp_mpi.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_omp_mpi -lz

# ----------------------------------------------------------------------------------------------------------------------
opp_kernels_cuda:
		 $(NVCC) $(NVCCFLAGS) -DUSE_CUDA -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda.o $(ALL_INC) $(CUDA_INC)
cabana_opp_cuda:
		 $(CPP) $(CPPFLAGS) -DUSE_CUDA -c cabana_opp.cpp -o $(OBJ)/cabana_opp_cuda.o $(ALL_INC)
cuda: mklib cabana_opp_cuda opp_kernels_cuda
	     $(CPP) $(CPPFLAGS) -DUSE_CUDA -o $(BIN)/cuda \
		 	$(OBJ)/cabana_opp_cuda.o \
			$(OBJ)/opp_kernels_cuda.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda -lz

# ----------------------------------------------------------------------------------------------------------------------
opp_kernels_cuda_mpi:
		 $(NVCC) $(NVCCFLAGS) -DUSE_CUDA -DUSE_MPI -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
cabana_opp_cuda_mpi:
	     $(MPICPP) $(CPPFLAGS) -DUSE_CUDA -DUSE_MPI -c cabana_opp.cpp -o $(OBJ)/cabana_opp_cuda_mpi.o $(ALL_INC)
cuda_mpi: mklib cabana_opp_cuda_mpi opp_kernels_cuda_mpi
	     $(MPICPP) $(CPPFLAGS) -DUSE_CUDA -DUSE_MPI -o $(BIN)/cuda_mpi \
			$(OBJ)/cabana_opp_cuda_mpi.o \
			$(OBJ)/opp_kernels_cuda_mpi.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi -lz

# ----------------------------------------------------------------------------------------------------------------------
opp_kernels_sycl:
		 $(SYCLCC) $(SYCLFLAGS) -c sycl/opp_kernels.cpp -o $(OBJ)/opp_kernels_sycl.o $(ALL_INC)
cabana_opp_sycl:
		 $(SYCLCC) $(SYCLFLAGS) -c cabana_opp.cpp -o $(OBJ)/cabana_opp_sycl.o $(ALL_INC)
sycl: mklib cabana_opp_sycl opp_kernels_sycl
	     $(SYCLCC) $(SYCLFLAGS) -o $(BIN)/sycl \
		 	$(OBJ)/cabana_opp_sycl.o \
			$(OBJ)/opp_kernels_sycl.o \
			$(ALL_LIBS) $(CUDA_LIB) -lopp_sycl -lz

# ----------------------------------------------------------------------------------------------------------------------
cabana_opp_sycl_mpi:
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c cabana_opp.cpp -o $(OBJ)/cabana_opp_sycl_mpi.o $(ALL_INC)
opp_kernels_sycl_mpi:
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c sycl/opp_kernels.cpp -o $(OBJ)/opp_kernels_sycl_mpi.o $(ALL_INC)
sycl_mpi: mklib opp_kernels_sycl_mpi cabana_opp_sycl_mpi
	     $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -o $(BIN)/sycl_mpi \
		 	$(OBJ)/cabana_opp_sycl_mpi.o \
			$(OBJ)/opp_kernels_sycl_mpi.o \
			$(ALL_LIBS) $(CUDA_LIB) -lopp_sycl_mpi -lz

# ----------------------------------------------------------------------------------------------------------------------
opp_kernels_hip:
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip.o $(ALL_INC) $(HIP_INC)
cabana_opp_hip:
	     $(MPICPP) $(CPPFLAGS) -DUSE_HIP -c cabana_opp.cpp -o $(OBJ)/cabana_opp_hip.o $(ALL_INC) $(HIP_INC)
hip: mklib cabana_opp_hip opp_kernels_hip
	     $(MPICPP) $(CPPFLAGS) -DUSE_HIP -o $(BIN)/hip \
			$(OBJ)/cabana_opp_hip.o \
			$(OBJ)/opp_kernels_hip.o \
			$(ALL_LIBS) $(HIP_LIB) -lopp_hip -lz

# ----------------------------------------------------------------------------------------------------------------------
opp_kernels_hip_mpi:
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -DUSE_MPI -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip_mpi.o $(ALL_INC) $(HIP_INC)
cabana_opp_hip_mpi:
	     $(MPICPP) $(CPPFLAGS) -DUSE_HIP -DUSE_MPI -c cabana_opp.cpp -o $(OBJ)/cabana_opp_hip_mpi.o $(ALL_INC) $(HIP_INC)
hip_mpi: mklib opp_kernels_hip_mpi cabana_opp_hip_mpi
	     $(MPICPP) $(CPPFLAGS) -DUSE_HIP -o $(BIN)/hip_mpi \
			$(OBJ)/cabana_opp_hip_mpi.o \
			$(OBJ)/opp_kernels_hip_mpi.o \
			$(ALL_LIBS) $(HIP_LIB) -lopp_hip_mpi -lz

# ----------------------------------------------------------------------------------------------------------------------
mpi_hdf5: mklib
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c cabana_op_hdf5.cpp -o $(OBJ)/cabana_op_hdf5.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/mpikernels.cpp -o $(OBJ)/mpikernels.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c $(OPP_PATH)/src/hdf5/opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi_hdf5 \
		 	$(OBJ)/cabana_op_hdf5.o \
			$(OBJ)/opp_mpi_hdf5.o \
			$(OBJ)/mpikernels.o \
			-lopp_mpi $(ALL_INC) $(ALL_LIBS) -lz -lhdf5

# ----------------------------------------------------------------------------------------------------------------------
mpi_hdf5_conv: mklib
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c cabana_convert_hdf5.cpp -o $(OBJ)/cabana_convert_hdf5.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi_hdf5_conv \
		 	$(OBJ)/cabana_convert_hdf5.o \
			$(OBJ)/opp_mpi_hdf5.o \
			-lopp_mpi $(ALL_INC) $(ALL_LIBS) -lz -lhdf5

# ----------------------------------------------------------------------------------------------------------------------
clean: 
	     rm -rf *.o *.d *.a
	     rm -rf $(OBJ)/*
	     rm -rf $(BIN)/*
	     rm -rf $(LIB)/*
		 
