
# module load gnu-6.3.0/gnu-6.3.0
# module load gnu-6.3.0/ompi-1.10.6
# module load cuda/toolkit-10.2.89

# module load intel-oneapi-2021.3/intel-classic-2021.3 intel-oneapi-2021.3/impi-classic-2021.3
# module load cuda/toolkit-10.2.89

# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3/lib:$LD_LIBRARY_PATH
# export LD_LIBRARY_PATH=/ext-home/zl/phd/OP-PIC/app_fempic_oppic/lib:/ext-home/zl/lib_install/petsc-3.18.2-oneapi-2021.3-release/lib:$LD_LIBRARY_PATH

# scorep-mpiicpc  -qopt-report=2 -qopt-report-phase=vec -xHOST     --generate-line-info  -fsanitize=address
# -qopt-zmm-usage=high    -march=znver3 -mtune=znver3

SRC = src
INC = include
LIB = lib
OBJ = obj

CPP        = $(CC_COMPILER)
MPICPP	   = $(MPI_COMPILER)
NVCC       = nvcc
HIPCC      = $(HIP_COMPILER)
SYCLCC     = icpx

CPPFLAGS   = -std=c++17 -Wall -fPIC -march=native $(CPPFLAGS_ADD)
NVCCFLAGS  = -std=c++17 --extended-lambda -m64 -Xptxas -dlcm=ca $(NVCCFLAGS_ADD)
HIPCCFLAGS = $(HIPCCFLAGS_ADD) -std=c++17 -isystem -D__HIP_PLATFORM_AMD__ -munsafe-fp-atomics -march=native 
SYCLFLAGS  = -w -fsycl -std=c++17 $(SYCLFLAGS_ADD) -DUSE_SYCL

OPP_INC    = -I$(OPP_PATH)/include

DEBUG ?= 0
ifeq ($(DEBUG), 1)
	CPPFLAGS   += -O0 -g
	NVCCFLAGS  += -O0 -g -lineinfo
	HIPCCFLAGS += -O0 -g -lineinfo
	SYCLFLAGS  += -O0
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

EXT_HALOS ?= 1
ifeq ($(EXT_HALOS), 0)
	CPPFLAGS   += -DNO_EXTENDED_HALOS
	NVCCFLAGS  += -DNO_EXTENDED_HALOS
	HIPCCFLAGS += -DNO_EXTENDED_HALOS
	SYCLFLAGS  += -DNO_EXTENDED_HALOS
endif

PETSC ?= 0
ifeq ($(PETSC), 1)
	CPPFLAGS   += -DUSE_PETSC
	NVCCFLAGS  += -DUSE_PETSC
	HIPCCFLAGS += -DUSE_PETSC
	SYCLFLAGS  += -DUSE_PETSC
	OPP_INC    += -I$(PETSC_INSTALL_PATH)/include
endif

PARMETIS ?= 0
ifeq ($(PARMETIS), 1)
	CPPFLAGS  += -DHAVE_PARMETIS
endif

CUDA_INC   += -I$(CUDA_INSTALL_PATH)/include
HIP_INC   += -I$(ROCM_INSTALL_DIR)/include -I$(ROCM_THRUST_DIR)/include -I$(ROCM_PRIM_DIR)/include

.PHONY: clean mklib

mklib:
	@mkdir -p $(LIB) $(OBJ)

# ------------------------------------------------------------------------------------------
SRC_SEQ_CPP = 	$(SRC)/core/opp_lib_core.cpp \
			  	$(SRC)/core/opp_util.cpp \
			  	$(SRC)/core/opp_params.cpp \
			  	$(SRC)/seq/opp_seq.cpp
OBJ_SEQ_CPP = 	$(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+seq.o,$(SRC_SEQ_CPP)) \
				$(patsubst $(SRC)/seq/%.cpp,$(OBJ)/%+seq.o,$(SRC_SEQ_CPP))
$(OBJ)/%+seq.o: $(SRC)/core/%.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@ $(OPP_INC)
$(OBJ)/%+seq.o: $(SRC)/seq/%.cpp
	$(CPP) $(CPPFLAGS) -c $< -o $@ $(OPP_INC)

seq_compile: $(OBJ_SEQ_CPP)
seq: clean mklib
	$(MAKE) seq_compile
	ar -r $(LIB)/libopp_seq.a $(OBJ_SEQ_CPP)

# ------------------------------------------------------------------------------------------
SRC_OMP_CPP = 	$(SRC)/core/opp_lib_core.cpp \
			  	$(SRC)/core/opp_util.cpp \
			  	$(SRC)/core/opp_params.cpp \
			  	$(SRC)/omp/opp_omp.cpp
OBJ_OMP_CPP = 	$(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+omp.o,$(SRC_OMP_CPP)) \
				$(patsubst $(SRC)/omp/%.cpp,$(OBJ)/%+omp.o,$(SRC_OMP_CPP))
$(OBJ)/%+omp.o: $(SRC)/core/%.cpp
	$(CPP) $(CPPFLAGS) -fopenmp -c $< -o $@ $(OPP_INC)
$(OBJ)/%+omp.o: $(SRC)/omp/%.cpp
	$(CPP) $(CPPFLAGS) -fopenmp -c $< -o $@ $(OPP_INC)

omp_compile: $(OBJ_OMP_CPP)
omp: clean mklib
	$(MAKE) omp_compile
	ar -r $(LIB)/libopp_omp.a $(OBJ_OMP_CPP)

# ------------------------------------------------------------------------------------------
SRC_CUDA_CPP = 	$(SRC)/core/opp_lib_core.cpp \
           		$(SRC)/core/opp_util.cpp \
           		$(SRC)/core/opp_params.cpp
SRC_CUDA_CU = 	$(SRC)/cuda/opp_cuda.cu \
           		$(SRC)/cuda/opp_particle_sorter.cu \
           		$(SRC)/cuda/opp_cuda_halo.cu
OBJ_CUDA_CPP = 	$(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+cuda.o,$(SRC_CUDA_CPP))
OBJ_CUDA_CU = 	$(patsubst $(SRC)/cuda/%.cu,$(OBJ)/%+cuda.o,$(SRC_CUDA_CU))
$(OBJ)/%+cuda.o: $(SRC)/core/%.cpp
	$(CPP) $(CPPFLAGS) -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(CUDA_INC)
$(OBJ)/%+cuda.o: $(SRC)/cuda/%.cu
	$(NVCC) $(NVCCFLAGS) -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(CUDA_INC)

cuda_compile: $(OBJ_CUDA_CPP) $(OBJ_CUDA_CU)
cuda: clean mklib
	$(MAKE) cuda_compile
	ar -r $(LIB)/libopp_cuda.a $(OBJ_CUDA_CPP) $(OBJ_CUDA_CU)

# ------------------------------------------------------------------------------------------
SRC_SYCL_CPP = 	$(SRC)/sycl/opp_sycl.cpp \
				$(SRC)/sycl/opp_particle_organize.cpp \
				$(SRC)/sycl/opp_sycl_halo.cpp \
				$(SRC)/core/opp_lib_core.cpp \
				$(SRC)/core/opp_util.cpp \
				$(SRC)/core/opp_params.cpp
OBJ_SYCL_CPP =  $(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+sycl.o,$(SRC_SYCL_CPP)) \
				$(patsubst $(SRC)/sycl/%.cpp,$(OBJ)/%+sycl.o,$(SRC_SYCL_CPP))
$(OBJ)/%+sycl.o: $(SRC)/core/%.cpp
	$(SYCLCC) $(SYCLFLAGS) -c $< -o $@ $(OPP_INC)
$(OBJ)/%+sycl.o: $(SRC)/sycl/%.cpp
	$(SYCLCC) $(SYCLFLAGS) -c $< -o $@ $(OPP_INC)

sycl_compile: $(OBJ_SYCL_CPP)
sycl: clean mklib
	$(MAKE) sycl_compile
	ar -r $(LIB)/libopp_sycl.a $(OBJ_SYCL_CPP)

# ------------------------------------------------------------------------------------------
SRC_HIP_CPP = 	$(SRC)/core/opp_lib_core.cpp \
				$(SRC)/core/opp_util.cpp \
				$(SRC)/core/opp_params.cpp \
				$(SRC)/hip/opp_hip.cpp \
				$(SRC)/hip/opp_particle_sorter.cpp \
				$(SRC)/hip/opp_hip_halo.cpp 
OBJ_HIP_CPP =  	$(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+hip.o,$(SRC_HIP_CPP)) \
				$(patsubst $(SRC)/hip/%.cpp,$(OBJ)/%+hip.o,$(SRC_HIP_CPP))
$(OBJ)/%+hip.o: $(SRC)/core/%.cpp
	$(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(HIP_INC)
$(OBJ)/%+hip.o: $(SRC)/hip/%.cpp
	$(HIPCC) $(HIPCCFLAGS) -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(HIP_INC)

hip_compile: $(OBJ_HIP_CPP)
hip: clean mklib
	$(MAKE) hip_compile
	ar -r $(LIB)/libopp_hip.a $(OBJ_HIP_CPP)

# ------------------------------------------------------------------------------------------
SRC_MPI_CPP = 	$(SRC)/core/opp_lib_core.cpp \
				$(SRC)/core/opp_util.cpp \
				$(SRC)/core/opp_params.cpp \
				$(SRC)/mpi/opp_mpi_halo.cpp \
				$(SRC)/mpi/opp_mpi_halo_core.cpp \
				$(SRC)/mpi/opp_mpi_partition.cpp \
				$(SRC)/mpi/opp_mpi.cpp \
				$(SRC)/mpi/opp_mpi_core.cpp \
				$(SRC)/mpi/opp_mpi_utils.cpp \
				$(SRC)/mpi/opp_mpi_particle_comm.cpp \
				$(SRC)/mpi/opp_mpi_double_ind_reducs.cpp 
OBJ_MPI_CPP =  	$(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+mpi.o,$(SRC_MPI_CPP)) \
				$(patsubst $(SRC)/mpi/%.cpp,$(OBJ)/%+mpi.o,$(SRC_MPI_CPP))
$(OBJ)/%+mpi.o: $(SRC)/core/%.cpp
	$(MPICPP) $(CPPFLAGS) -DUSE_MPI -c $< -o $@ $(OPP_INC)
$(OBJ)/%+mpi.o: $(SRC)/mpi/%.cpp
	$(MPICPP) $(CPPFLAGS) -DUSE_MPI -c $< -o $@ $(OPP_INC)

mpi_compile: $(OBJ_MPI_CPP)
mpi: clean mklib
	$(MAKE) mpi_compile
	ar -r $(LIB)/libopp_mpi.a $(OBJ_MPI_CPP)

# ------------------------------------------------------------------------------------------
SRC_SYCL_MPI_CPP = 	$(SRC)/sycl/opp_sycl.cpp \
					$(SRC)/sycl/opp_sycl_halo.cpp \
					$(SRC)/sycl/opp_particle_organize.cpp \
					$(SRC)/core/opp_lib_core.cpp \
					$(SRC)/core/opp_util.cpp \
					$(SRC)/core/opp_params.cpp \
					$(SRC)/mpi/opp_mpi_halo_core.cpp \
					$(SRC)/mpi/opp_mpi_partition.cpp \
					$(SRC)/mpi/opp_mpi_core.cpp \
					$(SRC)/mpi/opp_mpi_particle_comm.cpp \
					$(SRC)/mpi/opp_mpi_double_ind_reducs.cpp \
					$(SRC)/mpi/opp_mpi_utils.cpp
OBJ_SYCL_MPI_CPP =  $(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+sycl_mpi.o,$(SRC_SYCL_MPI_CPP)) \
					$(patsubst $(SRC)/sycl/%.cpp,$(OBJ)/%+sycl_mpi.o,$(SRC_SYCL_MPI_CPP)) \
					$(patsubst $(SRC)/mpi/%.cpp,$(OBJ)/%+sycl_mpi.o,$(SRC_SYCL_MPI_CPP))
$(OBJ)/%+sycl_mpi.o: $(SRC)/core/%.cpp
	$(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c $< -o $@ $(OPP_INC)
$(OBJ)/%+sycl_mpi.o: $(SRC)/sycl/%.cpp
	$(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c $< -o $@ $(OPP_INC)
$(OBJ)/%+sycl_mpi.o: $(SRC)/mpi/%.cpp
	$(MPICPP) $(SYCLFLAGS) -DUSE_MPI -c $< -o $@ $(OPP_INC)

sycl_mpi_compile: $(OBJ_SYCL_MPI_CPP)
sycl_mpi: clean mklib
	$(MAKE) sycl_mpi_compile
	ar -r $(LIB)/libopp_sycl_mpi.a $(OBJ_SYCL_MPI_CPP)

# ------------------------------------------------------------------------------------------
SRC_CUDA_MPI_CPP = 	$(SRC)/core/opp_lib_core.cpp \
					$(SRC)/core/opp_util.cpp \
					$(SRC)/core/opp_params.cpp \
					$(SRC)/mpi/opp_mpi_halo_core.cpp \
					$(SRC)/mpi/opp_mpi_partition.cpp \
					$(SRC)/mpi/opp_mpi_core.cpp \
					$(SRC)/mpi/opp_mpi_particle_comm.cpp \
					$(SRC)/mpi/opp_mpi_double_ind_reducs.cpp \
					$(SRC)/mpi/opp_mpi_utils.cpp
SRC_CUDA_MPI_CU = 	$(SRC)/cuda/opp_cuda.cu \
					$(SRC)/cuda/opp_particle_sorter.cu \
					$(SRC)/cuda/opp_cuda_halo.cu
OBJ_CUDA_MPI_CPP = 	$(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+cuda_mpi.o,$(SRC_CUDA_MPI_CPP)) \
					$(patsubst $(SRC)/mpi/%.cpp,$(OBJ)/%+cuda_mpi.o,$(SRC_CUDA_MPI_CPP))
OBJ_CUDA_MPI_CU = 	$(patsubst $(SRC)/cuda/%.cu,$(OBJ)/%+cuda_mpi.o,$(SRC_CUDA_MPI_CU))
$(OBJ)/%+cuda_mpi.o: $(SRC)/core/%.cpp
	$(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(CUDA_INC)
$(OBJ)/%+cuda_mpi.o: $(SRC)/mpi/%.cpp
	$(CPP) $(CPPFLAGS) -DUSE_MPI -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(CUDA_INC)
$(OBJ)/%+cuda_mpi.o: $(SRC)/cuda/%.cu
	$(NVCC) $(NVCCFLAGS) -DUSE_MPI -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(CUDA_INC)

cuda_mpi_compile: $(OBJ_CUDA_MPI_CPP) $(OBJ_CUDA_MPI_CU)
cuda_mpi: clean mklib
	$(MAKE) cuda_mpi_compile
	ar -r $(LIB)/libopp_cuda_mpi.a $(OBJ_CUDA_MPI_CPP) $(OBJ_CUDA_MPI_CU)

# ------------------------------------------------------------------------------------------
SRC_HIP_MPI_CPP = 	$(SRC)/core/opp_lib_core.cpp \
					$(SRC)/core/opp_util.cpp \
					$(SRC)/core/opp_params.cpp \
					$(SRC)/mpi/opp_mpi_halo_core.cpp \
					$(SRC)/mpi/opp_mpi_partition.cpp \
					$(SRC)/mpi/opp_mpi_core.cpp \
					$(SRC)/mpi/opp_mpi_particle_comm.cpp \
					$(SRC)/mpi/opp_mpi_double_ind_reducs.cpp \
					$(SRC)/mpi/opp_mpi_utils.cpp \
					$(SRC)/hip/opp_hip.cpp \
					$(SRC)/hip/opp_hip_halo.cpp \
					$(SRC)/hip/opp_particle_sorter.cpp
OBJ_HIP_MPI_CPP = 	$(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+hip_mpi.o,$(SRC_HIP_MPI_CPP)) \
					$(patsubst $(SRC)/mpi/%.cpp,$(OBJ)/%+hip_mpi.o,$(SRC_HIP_MPI_CPP)) \
					$(patsubst $(SRC)/hip/%.cpp,$(OBJ)/%+hip_mpi.o,$(SRC_HIP_MPI_CPP))
$(OBJ)/%+hip_mpi.o: $(SRC)/core/%.cpp
	$(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(HIP_INC)
$(OBJ)/%+hip_mpi.o: $(SRC)/mpi/%.cpp
	$(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(HIP_INC)
$(OBJ)/%+hip_mpi.o: $(SRC)/hip/%.cpp
	$(HIPCC) $(HIPCCFLAGS) -DUSE_MPI -DUSE_THRUST -c $< -o $@ $(OPP_INC) $(HIP_INC)

hip_mpi_compile: $(OBJ_HIP_MPI_CPP)
hip_mpi: clean mklib
	$(MAKE) hip_mpi_compile
	ar -r $(LIB)/libopp_hip_mpi.a $(OBJ_HIP_MPI_CPP)

# ------------------------------------------------------------------------------------------
SRC_OMP_MPI_CPP = 	$(SRC)/core/opp_lib_core.cpp \
					$(SRC)/core/opp_util.cpp \
					$(SRC)/core/opp_params.cpp \
					$(SRC)/omp/opp_omp.cpp \
					$(SRC)/mpi/opp_mpi_halo.cpp \
					$(SRC)/mpi/opp_mpi_halo_core.cpp \
					$(SRC)/mpi/opp_mpi_partition.cpp \
					$(SRC)/mpi/opp_mpi_core.cpp \
					$(SRC)/mpi/opp_mpi_particle_comm.cpp \
					$(SRC)/mpi/opp_mpi_double_ind_reducs.cpp \
					$(SRC)/mpi/opp_mpi_utils.cpp 
OBJ_OMP_MPI_CPP = 	$(patsubst $(SRC)/core/%.cpp,$(OBJ)/%+omp_mpi.o,$(SRC_OMP_MPI_CPP)) \
					$(patsubst $(SRC)/omp/%.cpp,$(OBJ)/%+omp_mpi.o,$(SRC_OMP_MPI_CPP)) \
					$(patsubst $(SRC)/mpi/%.cpp,$(OBJ)/%+omp_mpi.o,$(SRC_OMP_MPI_CPP))
$(OBJ)/%+omp_mpi.o: $(SRC)/core/%.cpp
	$(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -c $< -o $@ $(OPP_INC)
$(OBJ)/%+omp_mpi.o: $(SRC)/mpi/%.cpp
	$(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -c $< -o $@ $(OPP_INC)
$(OBJ)/%+omp_mpi.o: $(SRC)/omp/%.cpp
	$(MPICPP) $(CPPFLAGS) -DUSE_MPI -fopenmp -c $< -o $@ $(OPP_INC)

omp_mpi_compile: $(OBJ_OMP_MPI_CPP)
omp_mpi: clean mklib
	$(MAKE) omp_mpi_compile
	ar -r $(LIB)/libopp_omp_mpi.a $(OBJ_OMP_MPI_CPP)

# ------------------------------------------------------------------------------------------
clean: 
	     rm -rf *.o *.d *.a
	     rm -rf $(OBJ)/*

#  rm -rf $(LIB)/*

# ------------------------------------------------------------------------------------------