# export OPPIC_PATH=/ext-home/zl/phd/neso_test/OP-PIC/lib_oppic

SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = $(CC_COMPILER)
MPICPP	   = $(MPI_COMPILER)
NVCC       = nvcc
HIPCC      = $(HIP_COMPILER)

CPPFLAGS   = -std=c++14 -Wall -fPIC 
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -Wall -isystem -D__HIP_PLATFORM_AMD__

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g
	HIPCCFLAGS += -O3 -g
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
	HIPCCFLAGS += -DDEBUG_LOG
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
	HIPCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

CUDA ?= 1
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
	CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
endif

HIP ?= 1
ifeq ($(HIP), 1)
	HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include
endif

ALL_INC    += -I$(OPP_PATH)/include
ALL_LIBS   += -L$(OPP_PATH)/lib

.PHONY: clean mklib

ALL_TARGETS = clean seq omp cuda mpi

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

common: mklib 

seq: common
		 $(CPP) $(CPPFLAGS) -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c seq/seqkernels.cpp -o $(OBJ)/seqkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq $(OBJ)/advec_op.o $(OBJ)/seqkernels.o $(ALL_INC) $(ALL_LIBS) -lopp_seq

mpi: common
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/mpikernels.cpp -o $(OBJ)/mpikernels.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi $(OBJ)/advec_op.o $(OBJ)/mpikernels.o $(ALL_INC) $(ALL_LIBS) -lopp_mpi

omp: common
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp $(OBJ)/advec_op.o $(OBJ)/ompkernels.o $(ALL_INC) $(ALL_LIBS) -lopp_omp

omp_mpi: common
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_MPI -DUSE_OMP -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_MPI -DUSE_OMP -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp_mpi $(OBJ)/advec_op.o $(OBJ)/ompkernels.o $(ALL_INC) $(ALL_LIBS) -lopp_omp_mpi

cuda: common
		 $(CPP) $(CPPFLAGS) -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC)
		 $(NVCC) $(NVCCFLAGS) -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/cuda $(OBJ)/advec_op.o $(OBJ)/cudakernels.o $(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda

cuda_mpi: common
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC) $(CUDA_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -c cuda/cudakernels.cu -o $(OBJ)/cudakernels.o $(ALL_INC) $(CUDA_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi $(OBJ)/advec_op.o $(OBJ)/cudakernels.o $(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi

hip: mklib
		 $(HIPCC) $(HIPCCFLAGS) -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -c hip/hipkernels.cpp -o $(OBJ)/hipkernels.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -o $(BIN)/hip $(OBJ)/advec_op.o $(OBJ)/hipkernels.o $(ALL_INC) $(HIP_INC) $(ALL_LIBS) -lopp_hip

hip_mpi: common
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c advec_op.cpp -o $(OBJ)/advec_op.o $(ALL_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c hip/hipkernels.cpp -o $(OBJ)/hipkernels.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -o $(BIN)/hip_mpi $(OBJ)/advec_op.o $(OBJ)/hipkernels.o $(ALL_INC) $(HIP_INC) $(ALL_LIBS) -lopp_hip_mpi -lmpi

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
