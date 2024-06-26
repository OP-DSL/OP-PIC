
SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = $(CC_COMPILER)
MPICPP	   = $(MPI_COMPILER)
NVCC       = nvcc
HIPCC      = $(HIP_COMPILER)

CPPFLAGS   = -std=c++17 -Wall -fPIC -march=native
NVCCFLAGS  = -std=c++14 -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -Wall -isystem -D__HIP_PLATFORM_AMD__ -munsafe-fp-atomics -march=native

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g -lineinfo
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g -lineinfo
	HIPCCFLAGS += -O3 -g
endif

IEEE ?= 0
ifeq ($(IEEE), 1)
	CPPFLAGS   += -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
	HIPCCFLAGS  += -fmad=false -ftz=false -prec-div=true -prec-sqrt=true
endif

DEBUG_LOG ?= 0
ifeq ($(DEBUG_LOG), 1)
	CPPFLAGS   += -DDEBUG_LOG
	NVCCFLAGS  += -DDEBUG_LOG
	HIPCCFLAGS += -DDEBUG_LOG
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

PARMETIS ?= 0
ifeq ($(PARMETIS), 1)
	ALL_LIBS = -lparmetis -lmetis
	ALL_INC  += -DHAVE_PARMETIS
endif

ALL_INC    += -I$(OPP_PATH)/include -I$(HDF5_INSTALL_PATH)/include 
ALL_LIBS   += -L$(OPP_PATH)/lib -L$(HDF5_INSTALL_PATH)/lib

.PHONY: clean mklib

ALL_TARGETS = clean seq omp mpi

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

seq: mklib
		 $(CPP) $(CPPFLAGS) -c cabana_opp.cpp -o $(OBJ)/cabana_opp_seq.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c seq/opp_kernels.cpp -o $(OBJ)/opp_kernels_seq.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq \
			$(OBJ)/cabana_opp_seq.o \
			$(OBJ)/opp_kernels_seq.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_seq
omp: mklib
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c cabana_opp.cpp -o $(OBJ)/cabana_opp_omp.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp \
			$(OBJ)/cabana_opp_omp.o \
			$(OBJ)/opp_kernels_omp.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_omp

mpi: mklib
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c cabana_opp.cpp -o $(OBJ)/cabana_opp_mpi.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/opp_kernels.cpp -o $(OBJ)/opp_kernels_mpi.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi \
		 	$(OBJ)/cabana_opp_mpi.o \
			$(OBJ)/opp_kernels_mpi.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_mpi

omp_mpi: mklib
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -DUSE_OMP -c cabana_opp.cpp -o $(OBJ)/cabana_opp_omp_mpi.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp_mpi.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -o $(BIN)/omp_mpi \
			$(OBJ)/cabana_opp_omp_mpi.o \
			$(OBJ)/opp_kernels_omp_mpi.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_omp_mpi

mpi_hdf5: mklib
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c cabana_op_hdf5.cpp -o $(OBJ)/cabana_op_hdf5.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/mpikernels.cpp -o $(OBJ)/mpikernels.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c $(OPP_PATH)/src/hdf5/opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi_hdf5 \
		 	$(OBJ)/cabana_op_hdf5.o \
			$(OBJ)/opp_mpi_hdf5.o \
			$(OBJ)/mpikernels.o \
			-lopp_mpi $(ALL_INC) $(ALL_LIBS) -lhdf5

mpi_hdf5_conv: mklib
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c cabana_convert_hdf5.cpp -o $(OBJ)/cabana_convert_hdf5.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi_hdf5_conv \
		 	$(OBJ)/cabana_convert_hdf5.o \
			$(OBJ)/opp_mpi_hdf5.o \
			-lopp_mpi $(ALL_INC) $(ALL_LIBS) -lhdf5

cuda: mklib
		 $(NVCC) $(NVCCFLAGS) -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda.o $(ALL_INC) $(CUDA_INC)
		 $(CPP) $(CPPFLAGS) -c cabana_opp.cpp -o $(OBJ)/cabana_opp_cuda.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/cuda \
		 	$(OBJ)/cabana_opp_cuda.o \
			$(OBJ)/opp_kernels_cuda.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda

cuda_mpi: mklib
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c cabana_opp.cpp -o $(OBJ)/cabana_opp_cuda_mpi.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi \
			$(OBJ)/cabana_opp_cuda_mpi.o \
			$(OBJ)/opp_kernels_cuda_mpi.o \
			$(ALL_INC) $(CUDA_INC) $(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi

hip: mklib
		 $(HIPCC) $(HIPCCFLAGS) -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -c cabana_opp.cpp -o $(OBJ)/cabana_opp_hip.o $(ALL_INC)
	     $(HIPCC) -o $(BIN)/hip \
			$(OBJ)/cabana_opp_hip.o \
			$(OBJ)/opp_kernels_hip.o \
			$(ALL_LIBS) $(HIP_LIB) -lopp_hip

hip_mpi: mklib
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip_mpi.o $(ALL_INC) $(HIP_INC)
	     $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -c cabana_opp.cpp -o $(OBJ)/cabana_opp_hip_mpi.o $(ALL_INC)
	     $(HIPCC) -DUSE_MPI -o $(BIN)/hip_mpi \
			$(OBJ)/cabana_opp_hip_mpi.o \
			$(OBJ)/opp_kernels_hip_mpi.o \
			$(ALL_LIBS) $(HIP_LIB) -lopp_hip_mpi -lmpi

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
