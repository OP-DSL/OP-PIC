# source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi

SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = $(CC_COMPILER)
MPICPP	   = $(MPI_COMPILER)
NVCC       = nvcc
HIPCC      = $(HIP_COMPILER)
# -qopt-zmm-usage=high

CPPFLAGS   = -std=c++14 -Wall -fPIC -march=native 
NVCCFLAGS  = --extended-lambda -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD) -DUSE_THRUST 
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -isystem -D__HIP_PLATFORM_AMD__ -munsafe-fp-atomics -DUSE_THRUST -march=native

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g -lineinfo
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g
	HIPCCFLAGS += -O3 -g
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
endif

CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda

HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include

ALL_INC    = -I$(OPP_PATH)/include $(PETSC_INC) $(CUDA_INC)
ALL_LIBS   = -L$(OPP_PATH)/lib $(PETSC_LIB) $(CUDA_LIB)

.PHONY: clean mklib

ALL_TARGETS = clean omp seq cuda

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

cg: mklib 
	     $(CPP) $(CPPFLAGS) -c simpic_opp.cpp -o $(OBJ)/simpic_opp_cg.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic_cg.o $(ALL_INC)

seq: mklib 
	     $(CPP) $(CPPFLAGS) -c simpic_opp.cpp -o $(OBJ)/simpic_opp_seq.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic_seq.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c seq/opp_kernels.cpp -o $(OBJ)/opp_kernels_seq.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq \
			$(OBJ)/simpic_opp_seq.o \
			$(OBJ)/funcs_from_simpic_seq.o \
			$(OBJ)/opp_kernels_seq.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_seq

omp: mklib
	     $(CPP) $(CPPFLAGS) -fopenmp -c simpic_opp.cpp -o $(OBJ)/simpic_opp_omp.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic_omp.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -fopenmp -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp \
		 	$(OBJ)/simpic_opp_omp.o \
			$(OBJ)/funcs_from_simpic_omp.o \
			$(OBJ)/opp_kernels_omp.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_omp

mpi: mklib
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c simpic_opp.cpp -o $(OBJ)/simpic_opp_mpi.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/opp_kernels.cpp -o $(OBJ)/opp_kernels_mpi.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic_mpi.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi \
			$(OBJ)/simpic_opp_mpi.o \
			$(OBJ)/opp_kernels_mpi.o \
			$(OBJ)/field_solver_cpu_mpi.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_mpi 

cuda: mklib
		 $(CPP) $(CPPFLAGS) -c simpic_opp.cpp -o $(OBJ)/simpic_opp_cuda.o $(ALL_INC)
		 $(NVCC) $(NVCCFLAGS) -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda.o $(ALL_INC) $(CUDA_INC)	     
		 $(NVCC) $(NVCCFLAGS) -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic_cuda.o $(ALL_INC) $(CUDA_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/cuda \
		 	$(OBJ)/simpic_opp_cuda.o \
			$(OBJ)/opp_kernels_cuda.o \
			$(OBJ)/funcs_from_simpic_cuda.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda

cuda_mpi: mklib	
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c simpic_opp.cpp -o $(OBJ)/simpic_opp_cuda_mpi.o $(ALL_INC)
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda_mpi.o $(ALL_INC) $(CUDA_INC) 
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi \
		 	$(OBJ)/simpic_opp_cuda_mpi.o \
			$(OBJ)/opp_kernels_cuda_mpi.o \
			$(OBJ)/funcs_from_simpic_cuda_mpi.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi

hip: mklib		 
	     $(HIPCC) $(HIPCCFLAGS) -c simpic_opp.cpp -o $(OBJ)/simpic_opp_hip.o $(ALL_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip.o $(ALL_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic_hip.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) -o $(BIN)/hip \
		 	$(OBJ)/simpic_opp_hip.o \
			$(OBJ)/opp_kernels_hip.o \
			$(OBJ)/funcs_from_simpic_hip.o \
			$(ALL_LIBS) -lopp_hip

hip_mpi: mklib
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c simpic_opp.cpp -o $(OBJ)/simpic_opp_hip_mpi.o $(ALL_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip_mpi.o $(ALL_INC) $(HIP_INC)
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic_hip_mpi.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) -DUSE_MPI -o $(BIN)/hip_mpi \
		 	$(OBJ)/simpic_opp_hip_mpi.o \
			$(OBJ)/opp_kernels_hip_mpi.o \
			$(OBJ)/funcs_from_simpic_hip_mpi.o \
			$(ALL_LIBS) -lopp_hip_mpi -lmpi

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
