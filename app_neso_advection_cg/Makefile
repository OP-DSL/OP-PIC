
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
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -Wall -isystem -D__HIP_PLATFORM_AMD__
SYCLFLAGS  = -w -fsycl -std=c++17 $(SYCLFLAGS_ADD)

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g
	SYCLFLAGS  += -O0 -g
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g
	HIPCCFLAGS += -O3 -g
	SYCLFLAGS  += -O3 -g
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

CUDA ?= 0
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
	CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
endif

HIP ?= 0
ifeq ($(HIP), 1)
	HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include
	HIP_LIB   = -L$(ROCM_INSTALL_DIR)/lib -lamdhip64
endif

ALL_INC    += -I$(OPP_PATH)/include
ALL_LIBS   += -L$(OPP_PATH)/lib

.PHONY: clean mklib

mklib:
	@mkdir -p $(OBJ) $(BIN)

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_seq.o: advec_opp.cpp
	 $(CPP) $(CPPFLAGS) -c advec_opp.cpp -o $@ $(ALL_INC)
$(OBJ)/opp_kernels_seq.o: seq/opp_kernels.cpp
	 $(CPP) $(CPPFLAGS) -c $< -o $@ $(ALL_INC)
seq: mklib $(OBJ)/advec_opp_seq.o $(OBJ)/opp_kernels_seq.o
	 $(CPP) $(CPPFLAGS) -o $(BIN)/seq \
		$(OBJ)/advec_opp_seq.o \
		$(OBJ)/opp_kernels_seq.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_seq

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_mpi.o: advec_opp.cpp
	 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c advec_opp.cpp -o $@ $(ALL_INC)
$(OBJ)/opp_kernels_mpi.o: mpi/opp_kernels.cpp
	 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c $< -o $@ $(ALL_INC)
mpi: mklib $(OBJ)/advec_opp_mpi.o $(OBJ)/opp_kernels_mpi.o
	 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi \
		$(OBJ)/advec_opp_mpi.o \
		$(OBJ)/opp_kernels_mpi.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_mpi

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_omp.o: advec_opp.cpp
	 $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c advec_opp.cpp -o $@ $(ALL_INC)
$(OBJ)/opp_kernels_omp.o: omp/opp_kernels.cpp
	 $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c $< -o $@ $(ALL_INC)
omp: mklib $(OBJ)/advec_opp_omp.o $(OBJ)/opp_kernels_omp.o
	 $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp \
		$(OBJ)/advec_opp_omp.o \
		$(OBJ)/opp_kernels_omp.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_omp

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_omp_mpi.o: advec_opp.cpp
	 $(MPICPP) $(CPPFLAGS) -fopenmp -DUSE_MPI -DUSE_OMP -c advec_opp.cpp -o $@ $(ALL_INC)
$(OBJ)/opp_kernels_omp_mpi.o: omp/opp_kernels.cpp
	 $(MPICPP) $(CPPFLAGS) -fopenmp -DUSE_MPI -DUSE_OMP -c $< -o $@ $(ALL_INC)
omp_mpi: mklib $(OBJ)/advec_opp_omp_mpi.o $(OBJ)/opp_kernels_omp_mpi.o
	 $(MPICPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp_mpi \
		$(OBJ)/advec_opp_omp_mpi.o \
		$(OBJ)/opp_kernels_omp_mpi.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_omp_mpi

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_cuda.o: advec_opp.cpp
	 $(CPP) $(CPPFLAGS) -DUSE_CUDA -c advec_opp.cpp -o $@ $(ALL_INC)
$(OBJ)/opp_kernels_cuda.o: cuda/opp_kernels.cpp
	 $(NVCC) $(NVCCFLAGS) -DUSE_CUDA -c $< -o $@ $(ALL_INC) $(CUDA_INC)
cuda: mklib $(OBJ)/advec_opp_cuda.o $(OBJ)/opp_kernels_cuda.o
	 $(CPP) $(CPPFLAGS) -o $(BIN)/cuda \
		$(OBJ)/advec_opp_cuda.o \
		$(OBJ)/opp_kernels_cuda.o \
		$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_cuda_mpi.o: advec_opp.cpp
	 $(MPICPP) $(CPPFLAGS) -DUSE_CUDA -DUSE_MPI -c advec_opp.cpp -o $@ $(ALL_INC) $(CUDA_INC)
$(OBJ)/opp_kernels_cuda_mpi.o: cuda/opp_kernels.cpp
	 $(NVCC) $(NVCCFLAGS) -DUSE_CUDA -DUSE_MPI -c $< -o $@ $(ALL_INC) $(CUDA_INC)
cuda_mpi: mklib $(OBJ)/advec_opp_cuda_mpi.o $(OBJ)/opp_kernels_cuda_mpi.o
	 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi \
		$(OBJ)/advec_opp_cuda_mpi.o \
		$(OBJ)/opp_kernels_cuda_mpi.o \
		$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_sycl.o: advec_opp.cpp
	 $(SYCLCC) $(SYCLFLAGS) -DUSE_SYCL -c advec_opp.cpp -o $@ $(ALL_INC)
$(OBJ)/opp_kernels_sycl.o: sycl/opp_kernels.cpp
	 $(SYCLCC) $(SYCLFLAGS) -DUSE_SYCL -c $< -o $@ $(ALL_INC)
sycl: mklib $(OBJ)/advec_opp_sycl.o $(OBJ)/opp_kernels_sycl.o
	 $(SYCLCC) $(SYCLFLAGS) -o $(BIN)/sycl \
		$(OBJ)/advec_opp_sycl.o \
		$(OBJ)/opp_kernels_sycl.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_sycl

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_sycl_mpi.o: advec_opp.cpp
	 $(MPICPP) $(SYCLFLAGS) -DUSE_SYCL -DUSE_MPI -c advec_opp.cpp -o $@ $(ALL_INC)
$(OBJ)/opp_kernels_sycl_mpi.o: sycl/opp_kernels.cpp
	 $(MPICPP) $(SYCLFLAGS) -DUSE_SYCL -DUSE_MPI -c $< -o $@ $(ALL_INC)
sycl_mpi: mklib $(OBJ)/advec_opp_sycl_mpi.o $(OBJ)/opp_kernels_sycl_mpi.o
	 $(MPICPP) $(SYCLFLAGS) -o $(BIN)/sycl_mpi \
		$(OBJ)/advec_opp_sycl_mpi.o \
		$(OBJ)/opp_kernels_sycl_mpi.o \
		$(ALL_INC) $(ALL_LIBS) -lopp_sycl_mpi

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_hip.o: advec_opp.cpp
	 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -c advec_opp.cpp -o $@ $(ALL_INC) $(HIP_INC)
$(OBJ)/opp_kernels_hip.o: hip/opp_kernels.cpp
	 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -c $< -o $@ $(ALL_INC) $(HIP_INC)
hip: mklib $(OBJ)/advec_opp_hip.o $(OBJ)/opp_kernels_hip.o
	 $(HIPCC) $(HIPCCFLAGS) -o $(BIN)/hip \
		$(OBJ)/advec_opp_hip.o \
		$(OBJ)/opp_kernels_hip.o \
		$(ALL_INC) $(HIP_INC) $(ALL_LIBS) -lopp_hip

# ------------------------------------------------------------------------------------------
$(OBJ)/advec_opp_hip_mpi.o: advec_opp.cpp
	 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -DUSE_MPI -c advec_opp.cpp -o $@ $(ALL_INC) $(HIP_INC)
$(OBJ)/opp_kernels_hip_mpi.o: hip/opp_kernels.cpp
	 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -DUSE_MPI -c $< -o $@ $(ALL_INC) $(HIP_INC)
hip_mpi: mklib $(OBJ)/advec_opp_hip_mpi.o $(OBJ)/opp_kernels_hip_mpi.o
	 $(HIPCC) $(HIPCCFLAGS) -o $(BIN)/hip_mpi \
		$(OBJ)/advec_opp_hip_mpi.o \
		$(OBJ)/opp_kernels_hip_mpi.o \
		$(ALL_INC) $(HIP_INC) $(ALL_LIBS) -lopp_hip_mpi -lmpi

# ------------------------------------------------------------------------------------------
clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
