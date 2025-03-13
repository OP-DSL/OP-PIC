# source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi

SRC = .
INC = .
LIB = lib
OBJ = obj
BIN = bin

CPP        = mpicxx
MPICPP	   = mpicxx

CPPFLAGS   = -std=c++17 -Wall -fPIC
NVCCFLAGS  = -m64 -Xptxas -dlcm=ca -Xptxas=-v $(NVCCFLAGS_ADD)

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g -debug all
	NVCCFLAGS  += -O0 -g -lineinfo
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3
	NVCCFLAGS  += -O3
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

PETSC ?= 0
ifeq ($(PETSC), 1)
	CPPFLAGS   += -DUSE_PETSC
	PETSC_INC  = -I$(PETSC_INSTALL_PATH)/include
	PETSC_LIB  = -L$(PETSC_INSTALL_PATH)/lib -lpetsc
endif

CUDA ?= 0
ifeq ($(CUDA), 1)
	CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
	CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
endif

ALL_INC    = -I$(OPP_PATH)/include $(PETSC_INC) $(CUDA_INC)
ALL_LIBS   = -L$(OPP_PATH)/lib $(PETSC_LIB) $(CUDA_LIB)

.PHONY: clean mklib

ALL_TARGETS = clean omp genseq seq

all: $(ALL_TARGETS)

mklib:
	@mkdir -p $(OBJ) $(BIN)

seq: mklib 
	     $(CPP) $(CPPFLAGS) -c simpic_op.cpp -o $(OBJ)/simpic_op.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -c seq/seqkernels.cpp -o $(OBJ)/seqkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq \
			$(OBJ)/simpic_op.o \
			$(OBJ)/funcs_from_simpic.o \
			$(OBJ)/seqkernels.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_seq

omp: mklib
	     $(CPP) $(CPPFLAGS) -fopenmp -c simpic_op.cpp -o $(OBJ)/simpic_op.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -c funcs_from_simpic.cpp -o $(OBJ)/funcs_from_simpic.o $(ALL_INC)
		 $(CPP) $(CPPFLAGS) -fopenmp -c omp/ompkernels.cpp -o $(OBJ)/ompkernels.o $(ALL_INC)
	     $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp \
		 	$(OBJ)/simpic_op.o \
			$(OBJ)/funcs_from_simpic.o \
			$(OBJ)/ompkernels.o \
			$(ALL_INC) $(ALL_LIBS) -lopp_omp

seq1: mklib
		$(CPP) $(CPPFLAGS) -c simpic.cpp -o $(OBJ)/simpic.o $(ALL_INC)
		$(CPP) $(CPPFLAGS) -o $(BIN)/seq $(OBJ)/simpic.o $(OBJ)/funcs_from_simpic.o $(ALL_INC) $(ALL_LIBS) -lopp_seq

clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
