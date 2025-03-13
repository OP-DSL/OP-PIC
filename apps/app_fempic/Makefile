
# source /ext-home/zl/phd/OP-PIC/scripts/source_oneapi
# -qopt-zmm-usage=high
# scorep-mpiicpc  -qopt-report=2 -qopt-report-phase=vec -xHOST          --generate-line-info  -fsanitize=address -lasan

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

CPPFLAGS   = -std=c++17 -Wall -fPIC -march=native -DUSE_PETSC $(CPPFLAGS_ADD) 
NVCCFLAGS  = -std=c++17 --extended-lambda -m64 -Xptxas -dlcm=ca -DUSE_PETSC $(NVCCFLAGS_ADD) -DUSE_THRUST 
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -isystem -D__HIP_PLATFORM_AMD__ -munsafe-fp-atomics -DUSE_THRUST -DUSE_PETSC -march=native
SYCLFLAGS  = -w -fsycl -std=c++17 $(SYCLFLAGS_ADD) -DUSE_PETSC -DUSE_SYCL 

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g -lineinfo
	SYCLFLAGS  += -O0 -g
	DEBUG_LOG   = 1
else
	CPPFLAGS   += -O3 -g
	NVCCFLAGS  += -O3 -g
	HIPCCFLAGS += -O3 -g
	SYCLFLAGS  += -O3 -g
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
	SYCLFLAGS  += -DDEBUG_LOG
endif

PETSC_INC  = -I$(PETSC_INSTALL_PATH)/include
PETSC_LIB  = -Wl,-rpath,$(PETSC_INSTALL_PATH)/lib -L$(PETSC_INSTALL_PATH)/lib -lpetsc

CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include
CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda
# for demos with gcc
# CUDA_INC   = -I$(CUDA_INSTALL_PATH)/include -I$(MPI_INSTALL_PATH)/include
# CUDA_LIB   = -L$(CUDA_INSTALL_PATH)/lib64 -lcudart -lcuda -I$(MPI_INSTALL_PATH)/lib -lmpi

HIP_INC   = -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include
HIP_LIB   = -L$(ROCM_INSTALL_DIR)/lib -lamdhip64
# for archer2 gpus
# HIP_INC   = -I${CRAY_MPICH_DIR}/include
# HIP_LIB   = -L$(CRAY_MPICH_DIR)/lib -lmpi ${PE_MPICH_GTL_DIR_amd_gfx90a} ${PE_MPICH_GTL_LIBS_amd_gfx90a} -lamdhip64

PARMETIS ?= 0
ifeq ($(PARMETIS), 1)
	PARMETIS_LIB = -L$(PARMETIS_INSTALL_PATH)/lib -lparmetis -lmetis
	ALL_INC  += -I$(PARMETIS_INSTALL_PATH)/include -DHAVE_PARMETIS
endif

HDF5 ?= 0
ifeq ($(HDF5), 1)
	ALL_INC  += -I$(HDF5_INSTALL_PATH)/include 
	ALL_LIBS += -L$(HDF5_INSTALL_PATH)/lib
endif

ALL_INC    += -I$(OPP_PATH)/include $(PETSC_INC)
ALL_LIBS   += -L$(OPP_PATH)/lib $(PETSC_LIB) -lz

.PHONY: clean mklib

mklib:
	@mkdir -p $(OBJ) $(BIN)

# ----------------------------------------------------------------------------------------------------------------------
fempic_opp+seq: 
	     $(CPP) $(CPPFLAGS) -c fempic_opp.cpp -o $(OBJ)/fempic_opp_seq.o $(ALL_INC)
opp_kernels+seq: 
		 $(CPP) $(CPPFLAGS) -c seq/opp_kernels.cpp -o $(OBJ)/opp_kernels_seq.o $(ALL_INC)
field_solver+seq: 
		 $(CPP) $(CPPFLAGS) -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_cpu_seq.o $(ALL_INC)

seq: mklib fempic_opp+seq opp_kernels+seq field_solver+seq
	     $(CPP) $(CPPFLAGS) -o $(BIN)/seq \
		 	$(OBJ)/fempic_opp_seq.o \
			$(OBJ)/opp_kernels_seq.o \
			$(OBJ)/field_solver_cpu_seq.o \
			$(ALL_LIBS) -lopp_seq_pet -lz

# ----------------------------------------------------------------------------------------------------------------------
opp_kernels+omp: 
		 $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp.o $(ALL_INC)
fempic_opp+omp: 
	     $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c fempic_opp.cpp -o $(OBJ)/fempic_opp_omp.o $(ALL_INC)
field_solver+omp: 
		 $(CPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_cpu_omp.o $(ALL_INC)		 

omp: mklib opp_kernels+omp fempic_opp+omp field_solver+omp
		 $(CPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp \
			$(OBJ)/fempic_opp_omp.o \
			$(OBJ)/opp_kernels_omp.o \
			$(OBJ)/field_solver_cpu_omp.o \
			$(ALL_LIBS) -lopp_omp_pet -lz

# ----------------------------------------------------------------------------------------------------------------------
fempic_opp+omp_mpi: 
	     $(MPICPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -DUSE_MPI -c fempic_opp.cpp -o $(OBJ)/fempic_opp_omp_mpi.o $(ALL_INC)
field_solver+omp_mpi: 
		 $(MPICPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -DUSE_MPI -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_cpu_omp_mpi.o $(ALL_INC)
opp_kernels+omp_mpi: 
		 $(MPICPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -DUSE_MPI -c omp/opp_kernels.cpp -o $(OBJ)/opp_kernels_omp_mpi.o $(ALL_INC)
opp_mpi_hdf5+omp_mpi:
		 $(MPICPP) $(CPPFLAGS) -fopenmp -DUSE_OMP -DUSE_MPI -c $(OPP_PATH)/src/hdf5/opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5_omp_mpi.o $(ALL_INC)

omp_mpi: mklib opp_kernels+omp_mpi field_solver+omp_mpi fempic_opp+omp_mpi
		 $(MPICPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp_mpi \
		 	$(OBJ)/fempic_opp_omp_mpi.o \
			$(OBJ)/opp_kernels_omp_mpi.o \
			$(OBJ)/field_solver_cpu_omp_mpi.o \
			$(ALL_LIBS) -lopp_omp_mpi_pet -lz

omp_mpi_hdf5: mklib opp_kernels+omp_mpi field_solver+omp_mpi fempic_opp+omp_mpi opp_mpi_hdf5+omp_mpi
	     $(MPICPP) $(CPPFLAGS) -fopenmp -o $(BIN)/omp_mpi_hdf5 \
		 	$(OBJ)/fempic_opp_omp_mpi.o \
			$(OBJ)/opp_kernels_omp_mpi.o \
			$(OBJ)/field_solver_cpu_omp_mpi.o \
			$(OBJ)/opp_mpi_hdf5_omp_mpi.o \
			$(ALL_LIBS) -lopp_omp_mpi_pet -lz -lhdf5

# ----------------------------------------------------------------------------------------------------------------------
fempic_opp+mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c fempic_opp.cpp -o $(OBJ)/fempic_opp_mpi.o $(ALL_INC)
opp_kernels+mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c mpi/opp_kernels.cpp -o $(OBJ)/opp_kernels_mpi.o $(ALL_INC)
field_solver+mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_cpu_mpi.o $(ALL_INC)
fempic_hdf5_opp+mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c fempic_hdf5_opp.cpp -o $(OBJ)/fempic_hdf5_opp_mpi.o $(ALL_INC)
opp_mpi_hdf5+mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c $(OPP_PATH)/src/hdf5/opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5_mpi_hdf5.o $(ALL_INC)

mpi: mklib fempic_opp+mpi opp_kernels+mpi field_solver+mpi
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi \
			$(OBJ)/fempic_opp_mpi.o \
			$(OBJ)/opp_kernels_mpi.o \
			$(OBJ)/field_solver_cpu_mpi.o \
			$(ALL_LIBS) -lopp_mpi_pet -lz $(PARMETIS_LIB)

mpi_hdf5: mklib fempic_hdf5_opp+mpi opp_kernels+mpi field_solver+mpi opp_mpi_hdf5+mpi
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi_hdf5 \
		 	$(OBJ)/fempic_hdf5_opp_mpi.o  \
		 	$(OBJ)/opp_mpi_hdf5_mpi_hdf5.o  \
		 	$(OBJ)/opp_kernels_mpi.o \
		 	$(OBJ)/field_solver_cpu_mpi.o \
		 	$(ALL_LIBS) -lopp_mpi_pet -lz -lhdf5 $(PARMETIS_LIB)

# ----------------------------------------------------------------------------------------------------------------------
opp_kernels+cuda:
		 $(NVCC) $(NVCCFLAGS) -DUSE_CUDA -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda.o $(ALL_INC) $(CUDA_INC)
fempic_opp+cuda:
	     $(CPP) $(CPPFLAGS) -DUSE_CUDA -c fempic_opp.cpp -o $(OBJ)/fempic_opp_cuda.o $(ALL_INC)
field_solver+cuda:
		 $(NVCC) -x cu $(NVCCFLAGS) -DUSE_CUDA -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_cuda.o $(ALL_INC) $(CUDA_INC)

cuda: mklib opp_kernels+cuda fempic_opp+cuda field_solver+cuda
	     $(CPP) $(CPPFLAGS) -o $(BIN)/cuda \
		 	$(OBJ)/fempic_opp_cuda.o \
			$(OBJ)/opp_kernels_cuda.o \
			$(OBJ)/field_solver_cuda.o \
			$(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_pet -lz

# ----------------------------------------------------------------------------------------------------------------------
opp_kernels+cuda_mpi:	
		 $(NVCC) $(NVCCFLAGS) -DUSE_MPI -DUSE_CUDA -c cuda/opp_kernels.cu -o $(OBJ)/opp_kernels_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
fempic_opp+cuda_mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -DUSE_CUDA -c fempic_opp.cpp -o $(OBJ)/fempic_opp_cuda_mpi.o $(ALL_INC)
field_solver+cuda_mpi:
		 $(NVCC) -x cu $(NVCCFLAGS) -DUSE_MPI -DUSE_CUDA -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_cuda_mpi.o $(ALL_INC) $(CUDA_INC)
opp_mpi_hdf5+cuda_mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -DUSE_CUDA -c $(OPP_PATH)/src/hdf5/opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5_cuda_mpi_hdf5.o $(ALL_INC)

cuda_mpi: mklib opp_kernels+cuda_mpi fempic_opp+cuda_mpi field_solver+cuda_mpi
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi \
		 	$(OBJ)/fempic_opp_cuda_mpi.o \
			$(OBJ)/opp_kernels_cuda_mpi.o \
			$(OBJ)/field_solver_cuda_mpi.o \
			$(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi_pet -lz

cuda_mpi_hdf5: mklib opp_kernels+cuda_mpi fempic_opp+cuda_mpi field_solver+cuda_mpi opp_mpi_hdf5+cuda_mpi
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/cuda_mpi_hdf5 \
		 	$(OBJ)/field_solver_cuda_mpi_hdf5.o \
			$(OBJ)/fempic_hdf5_opp_cuda_mpi_hdf5.o \
			$(OBJ)/opp_mpi_hdf5_cuda_mpi_hdf5.o \
			$(OBJ)/opp_kernels_cuda_mpi_hdf5.o \
			$(ALL_LIBS) $(CUDA_LIB) -lopp_cuda_mpi_pet -lz -lhdf5

# ----------------------------------------------------------------------------------------------------------------------
fempic_opp+hip:
	     $(MPICPP) $(CPPFLAGS) -DUSE_HIP -c fempic_opp.cpp -o $(OBJ)/fempic_opp_hip.o $(ALL_INC) $(HIP_INC)
opp_kernels+hip:
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip.o $(ALL_INC) $(HIP_INC)
field_solver+hip:
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_HIP -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_hip.o $(ALL_INC) $(HIP_INC)

hip: mklib fempic_opp+hip opp_kernels+hip field_solver+hip
		 $(MPICPP) $(CPPFLAGS) -o $(BIN)/hip \
		 	$(OBJ)/fempic_opp_hip.o \
			$(OBJ)/opp_kernels_hip.o \
			$(OBJ)/field_solver_hip.o \
			$(ALL_LIBS) $(HIP_LIB) -lopp_hip_pet -lz

# ----------------------------------------------------------------------------------------------------------------------
fempic_opp+hip_mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -DUSE_HIP -c fempic_opp.cpp -o $(OBJ)/fempic_opp_hip_mpi.o $(ALL_INC) $(HIP_INC)
opp_kernels+hip_mpi:
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -DUSE_HIP -c hip/opp_kernels.cpp -o $(OBJ)/opp_kernels_hip_mpi.o $(ALL_INC) $(HIP_INC)
field_solver+hip_mpi:
		 $(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -DUSE_HIP -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_hip_mpi.o $(ALL_INC) $(HIP_INC)
opp_mpi_hdf5+hip_mpi:
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -DUSE_HIP -c $(OPP_PATH)/src/hdf5/opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5_hip_mpi_hdf5.o $(ALL_INC) $(HIP_INC)

hip_mpi: mklib fempic_opp+hip_mpi opp_kernels+hip_mpi field_solver+hip_mpi
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/hip_mpi \
		 	$(OBJ)/fempic_opp_hip_mpi.o \
			$(OBJ)/opp_kernels_hip_mpi.o \
			$(OBJ)/field_solver_hip_mpi.o \
			$(ALL_LIBS) $(HIP_LIB) -lopp_hip_mpi_pet -lz

hip_mpi_hdf5: mklib fempic_opp+hip_mpi opp_kernels+hip_mpi field_solver+hip_mpi opp_mpi_hdf5+hip_mpi
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/hip_mpi_hdf5 \
		 	$(OBJ)/fempic_hdf5_opp_hip_mpi_hdf5.o \
			$(OBJ)/opp_kernels_hip_mpi_hdf5.o \
			$(OBJ)/opp_mpi_hdf5_hip_mpi_hdf5.o \
			$(OBJ)/field_solver_hip_mpi_hdf5.o \
			$(ALL_LIBS) $(HIP_LIB) -lopp_hip_mpi_pet -lz -lhdf5

# ----------------------------------------------------------------------------------------------------------------------
fempic_opp+sycl:
	     $(SYCLCC) $(SYCLFLAGS) -c fempic_opp.cpp -o $(OBJ)/fempic_opp_sycl.o $(ALL_INC)
opp_kernels+sycl:
		 $(SYCLCC) $(SYCLFLAGS) -c sycl/opp_kernels.cpp -o $(OBJ)/opp_kernels_sycl.o $(ALL_INC)
field_solver+sycl:
		 $(SYCLCC) $(SYCLFLAGS) -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_sycl.o $(ALL_INC)

sycl: mklib fempic_opp+sycl opp_kernels+sycl field_solver+sycl
		 $(SYCLCC) $(SYCLFLAGS) -o $(BIN)/sycl \
		 	$(OBJ)/fempic_opp_sycl.o \
			$(OBJ)/opp_kernels_sycl.o \
			$(OBJ)/field_solver_sycl.o \
			$(ALL_LIBS) $(SYCL_LIB) -lopp_sycl_pet -lz

# ----------------------------------------------------------------------------------------------------------------------
fempic_opp+sycl_mpi:
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c fempic_opp.cpp -o $(OBJ)/fempic_opp_sycl_mpi.o $(ALL_INC)
opp_kernels+sycl_mpi:
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c sycl/opp_kernels.cpp -o $(OBJ)/opp_kernels_sycl_mpi.o $(ALL_INC) 
field_solver+sycl_mpi:
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c field_solver/field_solver.cpp -o $(OBJ)/field_solver_sycl_mpi.o $(ALL_INC) 
fempic_opp_hdf5+sycl_mpi:
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c fempic_hdf5_opp.cpp -o $(OBJ)/fempic_hdf5_opp_sycl_mpi.o $(ALL_INC)
opp_mpi_hdf5+sycl_mpi:
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c $(OPP_PATH)/src/hdf5/opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5_sycl_mpi_hdf5.o $(ALL_INC)

sycl_mpi: mklib fempic_opp+sycl_mpi opp_kernels+sycl_mpi field_solver+sycl_mpi
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -o $(BIN)/sycl_mpi \
		 	$(OBJ)/fempic_opp_sycl_mpi.o \
			$(OBJ)/opp_kernels_sycl_mpi.o \
			$(OBJ)/field_solver_sycl_mpi.o \
			$(ALL_LIBS) $(SYCL_LIB) -lopp_sycl_mpi_pet -lz

sycl_mpi_hdf5: mklib fempic_opp_hdf5+sycl_mpi opp_kernels+sycl_mpi field_solver+sycl_mpi opp_mpi_hdf5+sycl_mpi
		 $(MPICPP) $(SYCLFLAGS) -DUSE_MPI -o $(BIN)/sycl_mpi_hdf5 \
		 	$(OBJ)/fempic_hdf5_opp_sycl_mpi.o \
			$(OBJ)/opp_kernels_sycl_mpi.o \
			$(OBJ)/opp_mpi_hdf5_sycl_mpi_hdf5.o \
			$(OBJ)/field_solver_sycl_mpi.o \
			$(ALL_LIBS) $(SYCL_LIB) -lopp_sycl_mpi_pet -lz -lhdf5

# ----------------------------------------------------------------------------------------------------------------------
mpi_hdf5_conv: mklib
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c fempic_convert_hdf5.cpp -o $(OBJ)/fempic_convert_hdf5.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c opp_mpi_hdf5.cpp -o $(OBJ)/opp_mpi_hdf5.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c fempic_ori/meshes.cpp -o $(OBJ)/meshes.o $(ALL_INC)
	     $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c fempic_ori/particles.cpp -o $(OBJ)/particles.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -c fempic_ori/maths.cpp -o $(OBJ)/maths.o $(ALL_INC)
		 $(MPICPP) $(CPPFLAGS) -DUSE_MPI -o $(BIN)/mpi_hdf5_conv_2 $(OBJ)/fempic_convert_hdf5.o $(OBJ)/meshes.o $(OBJ)/particles.o $(OBJ)/maths.o $(OBJ)/opp_mpi_hdf5.o -lopp_mpi_pet -lz $(ALL_INC) $(ALL_LIBS) -lhdf5

# ----------------------------------------------------------------------------------------------------------------------
clean: 
	     rm -f *.o *.d *.a
	     rm -f $(OBJ)/*
	     rm -f $(BIN)/*
	     rm -f $(LIB)/*
		 
