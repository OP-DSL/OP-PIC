
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

CPPFLAGS   = -std=c++17 -Wall -fPIC $(CPPFLAGS_ADD)
NVCCFLAGS  = -std=c++17 -m64 -Xptxas -dlcm=ca $(NVCCFLAGS_ADD)
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -Wall -isystem -D__HIP_PLATFORM_AMD__
SYCLFLAGS  = -w -fsycl -std=c++17 $(SYCLFLAGS_ADD)

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -G -g -lineinfo
	HIPCCFLAGS += -O0 -g
	SYCLFLAGS  += -O0 -g
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g -DLOG_HOPS
	NVCCFLAGS  += -O3 -g -DLOG_HOPS
	HIPCCFLAGS += -O3 -g -DLOG_HOPS
	SYCLFLAGS  += -O3 -g -DLOG_HOPS
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
	HIPCCFLAGS += -DDEBUG_LOG
	SYCLFLAGS += -DDEBUG_LOG
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
	HIPCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include -I$(MPI_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda -I$(MPI_INSTALL_PATH)/lib -lmpi

HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include
HIP_LIB   = -L$(ROCM_INSTALL_DIR)/lib -lamdhip64

ALL_INC    += -I$(OPP_PATH)/include
ALL_LIBS   += -L$(OPP_PATH)/lib

.PHONY: clean mklib

mklib:
	@mkdir -p $(OBJ) $(BIN)

# ------------------------------------------------------------------------------------------
advec_opp+seq:
	 $(CPP) $(CPPFLAGS) -c advec_opp.cpp -o $(OBJ)/advec_opp_seq.o $(ALL_INC)
opp_kernels+seq:
	 $(CPP) $(CPPFLAGS) -c seq/opp_kernels.cpp -o $(OBJ)/opp_kernels_seq.o $(ALL_INC)
seq: mklib advec_opp+seq opp_kernels+seq
	 $(CPP) $(CPPFLAGS) -o $(BIN)/seq \
		$(OBJ)/advec_opp_seq.o \
		$(OBJ)/opp_kernels_seq.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_seq -lz

# ------------------------------------------------------------------------------------------
advec_opp+mpi:
	 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c advec_opp.cpp -o $(OBJ)/advec_opp_mpi.o $(ALL_INC)
opp_kernels+mpi:
	 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/opp_kernels.cpp -o $(OBJ)/opp_kernels_mpi.o $(ALL_INC)
mpi: mklib advec_opp+mpi opp_kernels+mpi
	 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi \
		$(OBJ)/advec_opp_mpi.o \
		$(OBJ)/opp_kernels_mpi.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_mpi -lz

# ------------------------------------------------------------------------------------------
advec_opp+omp:
	 $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c advec_opp.cpp -o $(OBJ)/advec_opp_omp.o $(ALL_INC)
opp_kernels+omp:
	 $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp.o $(ALL_INC)
omp: mklib advec_opp+omp opp_kernels+omp
	 $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp \
		$(OBJ)/advec_opp_omp.o \
		$(OBJ)/opp_kernels_omp.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_omp -lz

# ------------------------------------------------------------------------------------------
advec_opp+omp_mpi:
	 $(MPICPP) $(CPPFLAGS) -fopenmp -DUSE_MPI -DUSE_OMP -c advec_opp.cpp -o $(OBJ)/advec_opp_omp_mpi.o $(ALL_INC)
opp_kernels+omp_mpi:
	 $(MPICPP) $(CPPFLAGS) -fopenmp -DUSE_MPI -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp_mpi.o $(ALL_INC)
omp_mpi: mklib advec_opp+omp_mpi opp_kernels+omp_mpi
	 $(MPICPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp_mpi \
		$(OBJ)/advec_opp_omp_mpi.o \
		$(OBJ)/opp_kernels_omp_mpi.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_omp_mpi -lz

# ------------------------------------------------------------------------------------------
advec_opp+cuda:
	 $(CPP) $(CPPFLAGS) -DUSE_CUDA -c advec_opp.cpp -o $(OBJ)/advec_opp_cuda.o $(ALL_INC)
opp_kernels+cuda:
	 $(NVCC) $(NVCCFLAGS) -DUSE_CUDA -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda.o $(ALL_INC) $(CUDA_INC)
cuda: mklib advec_opp+cuda opp_kernels+cuda
	 $(CPP) $(CPPFLAGS) -o $(BIN)/cuda \
		$(OBJ)/advec_opp_cuda.o \
		$(OBJ)/opp_kernels_cuda.o \
		$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda -lz

# ------------------------------------------------------------------------------------------
advec_opp+cuda_mpi:
	 $(MPICPP) $(CPPFLAGS) -DUSE_CUDA -DUSE_MPI -c advec_opp.cpp -o $(OBJ)/advec_opp_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
opp_kernels+cuda_mpi:
	 $(NVCC) $(NVCCFLAGS) -DUSE_CUDA -DUSE_MPI -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
cuda_mpi: mklib advec_opp+cuda_mpi opp_kernels+cuda_mpi
	 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi \
		$(OBJ)/advec_opp_cuda_mpi.o \
		$(OBJ)/opp_kernels_cuda_mpi.o \
		$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi -lz

# ------------------------------------------------------------------------------------------
advec_opp+sycl:
	 $(SYCLCC) $(SYCLFLAGS) -DUSE_SYCL -c advec_opp.cpp -o $(OBJ)/advec_opp_sycl.o $(ALL_INC)
opp_kernels+sycl:
	 $(SYCLCC) $(SYCLFLAGS) -DUSE_SYCL -c sycl/opp_kernels.cpp -o $(OBJ)/opp_kernels_sycl.o $(ALL_INC)
sycl: mklib advec_opp+sycl opp_kernels+sycl
	 $(SYCLCC) $(SYCLFLAGS) -o $(BIN)/sycl \
		$(OBJ)/advec_opp_sycl.o \
		$(OBJ)/opp_kernels_sycl.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_sycl -lz

# ------------------------------------------------------------------------------------------
advec_opp+sycl_mpi:
	 $(MPICPP) $(SYCLFLAGS) -DUSE_SYCL -DUSE_MPI -c advec_opp.cpp -o $(OBJ)/advec_opp_sycl_mpi.o $(ALL_INC)
opp_kernels+sycl_mpi:
	 $(MPICPP) $(SYCLFLAGS) -DUSE_SYCL -DUSE_MPI -c sycl/opp_kernels.cpp -o $(OBJ)/opp_kernels_sycl_mpi.o $(ALL_INC)
sycl_mpi: mklib advec_opp+sycl_mpi opp_kernels+sycl_mpi
	 $(MPICPP) $(SYCLFLAGS) -o $(BIN)/sycl_mpi \
		$(OBJ)/advec_opp_sycl_mpi.o \
		$(OBJ)/opp_kernels_sycl_mpi.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_sycl_mpi -lz

# ------------------------------------------------------------------------------------------
advec_opp+hip:
	 $(MPICPP) $(CPPFLAGS) -DUSE_HIP -c advec_opp.cpp -o $(OBJ)/advec_opp_hip.o $(ALL_INC) $(HIP_INC)
opp_kernels+hip:
	 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip.o $(ALL_INC) $(HIP_INC)
hip: mklib advec_opp+hip opp_kernels+hip
	 $(MPICPP) $(CPPFLAGS) -o $(BIN)/hip \
		$(OBJ)/advec_opp_hip.o \
		$(OBJ)/opp_kernels_hip.o \
		$(ALL_INC) $(HIP_INC) $(ALL_LIBS) -lopp_hip -lz $(HIP_LIB)

# ------------------------------------------------------------------------------------------
advec_opp+hip_mpi:
	 $(MPICPP) $(CPPFLAGS) -DUSE_HIP -DUSE_MPI -c advec_opp.cpp -o $(OBJ)/advec_opp_hip_mpi.o $(ALL_INC) $(HIP_INC)
opp_kernels+hip_mpi:
	 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -DUSE_MPI -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip_mpi.o $(ALL_INC) $(HIP_INC)
hip_mpi: mklib advec_opp+hip_mpi opp_kernels+hip_mpi
	 $(MPICPP) $(CPPFLAGS) -o $(BIN)/hip_mpi \
		$(OBJ)/advec_opp_hip_mpi.o \
		$(OBJ)/opp_kernels_hip_mpi.o \
		$(ALL_INC) $(HIP_INC) $(ALL_LIBS) -lopp_hip_mpi -lmpi -lz $(HIP_LIB)

# ------------------------------------------------------------------------------------------
clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
