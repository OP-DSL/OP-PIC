
# module load gnu-6.3.0/gnu-6.3.0
# module load gnu-6.3.0/ompi-1.10.6
# module load cuda/toolkit-10.2.89

# module load intel-oneapi-2021.3/intel-classic-2021.3 intel-oneapi-2021.3/impi-classic-2021.3
# module load cuda/toolkit-10.2.89

# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3-release/lib:$LD_LIBRARY_PATH

SRC = src
INC = include
LIB = lib
OBJ = obj

CPP        = mpicxx
MPICPP	   = mpicxx

CPPFLAGS   = -std=c++11 -Wall -O3 -fopenmp -fPIC -fp-model strict -fp-model source -prec-div -prec-sqrt
NVCCFLAGS  = -O3 -m64 -Xptxas -dlcm=ca -Xptxas=-v -arch=compute_70 -code=sm_70 -prec-div=true -prec-sqrt=true -ftz=false

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   = -std=c++11 -Wall -O0 -g -fopenmp -fPIC -debug all -traceback -fp-model strict -fp-model source -prec-div -prec-sqrt
	NVCCFLAGS  = -O0 -g -m64 -Xptxas -dlcm=ca -Xptxas=-v -arch=compute_70 -code=sm_70 -prec-div=true -prec-sqrt=true -ftz=false
endif

CUDA_INSTALL_PATH = /opt/cuda/10.2.89/toolkit
CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart

# PETSC_PATH = /ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3-release
PETSC_PATH = /ext-home/zl/lib_install/petsc-3.18.3-oneapi-2021.3-cuda-rel
PETSC_INC  = -I$(PETSC_PATH)/include
PETSC_LIB  = -L$(PETSC_PATH)/lib -lpetsc

OP2_PATH   = /ext-home/zl/phd/OP2-Common/op2
OP2_INC    = # -I$(OP2_PATH)/include
OP2_LIB    = # -L$(OP2_PATH)/lib

OPPIC_PATH = /ext-home/zl/phd/OP-PIC/lib_oppic
ALL_INC        = -I$(OPPIC_PATH)/$(INC) $(PETSC_INC) $(CUDA_INC)
ALL_LIBS       = $(PETSC_LIB) $(CUDA_LIB)

.PHONY: clean mklib

ALL_TARGETS = clean oppic_seq oppic_omp oppic_cuda

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(LIB) $(OBJ)


oppic_core: mklib
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/core/oppic_lib_core.cpp -o $(OBJ)/oppic_lib_core.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/core/oppic_util.cpp -o $(OBJ)/oppic_util.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/core/trace.cpp -o $(OBJ)/trace.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/core/Params.cpp -o $(OBJ)/Params.o $(ALL_INC)

oppic_seq: mklib oppic_core
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/seq/oppic_seq.cpp -o $(OBJ)/oppic_seq.o $(ALL_INC)
	     ar -r $(LIB)/liboppic_seq.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/Params.o $(OBJ)/oppic_seq.o

oppic_omp: mklib oppic_core
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/omp/oppic_omp.cpp -o $(OBJ)/oppic_omp.o $(ALL_INC)
	     ar -r $(LIB)/liboppic_omp.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/oppic_omp.o

oppic_cuda: mklib oppic_core
	     $(CPP) $(CPPFLAGS) -c $(OPPIC_PATH)/$(SRC)/cuda/oppic_cuda.cpp -o $(OBJ)/oppic_cuda.o $(ALL_INC)
		 nvcc $(NVCCFLAGS) -c $(OPPIC_PATH)/$(SRC)/cuda/particle_sorter.cu -o $(OBJ)/particle_sorter.o $(ALL_INC)
	     ar -r $(LIB)/liboppic_cuda.a $(OBJ)/oppic_lib_core.o $(OBJ)/oppic_util.o $(OBJ)/trace.o $(OBJ)/oppic_cuda.o $(OBJ)/particle_sorter.o

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(LIB)/*
