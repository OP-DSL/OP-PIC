# ------------------------------------------------------------------------------

# This dockerfile compiles and runs both Mini-FEM-PIC and CabanaPIC applications
# on CPUs using an Ubuntu Operating System. Further instructions on deploying 
# (e.g. on GPUs) can be found below.

# Please copy and keep this Dockerfile and the extracted artifact in the same 
# folder prior running this script. Please check the below commented commands.

# mkdir op-pic_docker
# cd op-pic_docker
# copy this docker file to op-pic_docker folder
# wget https://zenodo.org/records/11876809/files/OP-PIC_Artifacts_icpp24.tar.gz
# tar -xvf OP-PIC_Artifacts_icpp24.tar.gz

# docker build -t op-pic:latest .
# docker images
# docker run -it <IMAGE ID> /bin/bash

# ------------------------------------------------------------------------------
FROM ubuntu

# Install necessary packages
RUN apt-get update && \
    apt-get install -y \
    gcc \
    g++ \
    openmpi-bin \
    libopenmpi-dev \
    make \
    git \
    wget \
    software-properties-common \
    petsc-dev \
    libhdf5-dev \
    python3 \
    python3-venv \
    vim

# switch user
RUN useradd -ms /bin/bash myuser
USER myuser

# setup compilers
ENV MPICH_CXX=/usr/bin/g++
ENV MPICH_CC=/usr/bin/gcc
ENV CC_COMPILER=mpicxx
ENV MPI_COMPILER=mpicxx

# setup user installed libs
ENV PETSC_INSTALL_PATH=/usr/lib/petsc

ENV HDF5_INSTALL_PATH=/usr/lib/aarch64-linux-gnu/hdf5/openmpi
# 'x86_64-linux-gnu' for x86_64, 'aarch64-linux-gnu' for Mac M1

# clone the OP-PIC library
WORKDIR /home/myuser/
RUN git clone https://github.com/OP-DSL/OP-PIC.git /home/myuser/OP-PIC 
WORKDIR /home/myuser/OP-PIC

# setup OP-PIC related environment variables
ENV OPP=/home/myuser/OP-PIC
ENV OPP_INSTALL_PATH=/home/myuser/OP-PIC/opp_install
ENV OPP_TRANSLATOR=/home/myuser/OP-PIC/opp_translator/translator
ENV OPP_PATH=/home/myuser/OP-PIC/opp_lib

# setup OP-PIC translator (one time process)
WORKDIR /home/myuser/OP-PIC/opp_translator
RUN chmod +x setup_venv.sh
RUN ./setup_venv.sh

# Compile Mini-FEM-PIC --------------------------------------------------------

# activate OP-PIC translator environment and code-generate Mini-FEM-PIC
WORKDIR /home/myuser/OP-PIC/apps/app_fempic
RUN . /home/myuser/OP-PIC/opp_translator/opp_venv/bin/activate && \
  python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths fempic.cpp && \
  python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths fempic_hdf5.cpp

# compile lib for Mini-FEM-PIC
WORKDIR /home/myuser/OP-PIC/opp_lib
RUN make PETSC=1 seq -j4
RUN make PETSC=1 mpi -j4
RUN make PETSC=1 omp -j4

# compile Mini-FEM-PIC application
WORKDIR /home/myuser/OP-PIC/apps/app_fempic
RUN make PETSC=1 seq -j4
RUN make PETSC=1 mpi -j4
RUN make PETSC=1 omp -j4
RUN make PETSC=1 mpi_hdf5 -j4 || \
    echo "Script is written for architecture x86_64, please check the arch using echo $(uname -m) and change HDF5_INSTALL_PATH if required" 

# Copy the unzipped artifacts
ADD OP-PIC_Artifacts /home/myuser/OP-PIC_Artifacts
USER root
RUN chmod -R 777 /home/myuser/OP-PIC_Artifacts/mesh_files
USER myuser

# extract the mesh files
WORKDIR /home/myuser/OP-PIC_Artifacts/mesh_files
RUN python3 extract_tar.py

# update the config file with the mesh directory
WORKDIR /home/myuser/OP-PIC/apps/app_fempic
RUN sed -i 's|STRING global_mesh  = .*|STRING global_mesh  = /home/myuser/OP-PIC_Artifacts/mesh_files/48000/mesh.dat|' /home/myuser/OP-PIC/apps/app_fempic/configs/coarse.param
RUN sed -i 's|STRING inlet_mesh   = .*|STRING inlet_mesh   = /home/myuser/OP-PIC_Artifacts/mesh_files/48000/inlet.dat|' /home/myuser/OP-PIC/apps/app_fempic/configs/coarse.param
RUN sed -i 's|STRING wall_mesh    = .*|STRING wall_mesh    = /home/myuser/OP-PIC_Artifacts/mesh_files/48000/wall.dat|' /home/myuser/OP-PIC/apps/app_fempic/configs/coarse.param
RUN sed -i 's|STRING hdf_filename = .*|STRING hdf_filename = /home/myuser/OP-PIC_Artifacts/mesh_files/box_48000.hdf5|' /home/myuser/OP-PIC/apps/app_fempic/configs/coarse.param
RUN sed -i 's|STRING rand_file    = .*|STRING rand_file    = /home/myuser/OP-PIC_Artifacts/mesh_files/random_100k.dat|' /home/myuser/OP-PIC/apps/app_fempic/configs/coarse.param

# Run Mini-FEM-PIC for 250 iterations with smaller configurations ------------

# sequential run with one CPU core
RUN bin/seq configs/coarse.param > run_fempic_seq.log

# MPI run with 4 CPU cores
RUN mpirun -np 4 bin/mpi configs/coarse.param > run_fempic_mpi.log

# OpenMP run with 4 CPU cores
ENV OMP_NUM_THREADS=4
ENV OMP_PROC_BIND=close
ENV OMP_PLACES=cores
RUN bin/omp configs/coarse.param > run_fempic_omp.log

# MPI run with 4 CPU cores
RUN mpirun -np 4 bin/mpi_hdf5 configs/coarse.param > run_fempic_mpi_hdf5.log

# Compile CabanaPIC ----------------------------------------------------------
  
# activate OP-PIC translator environment and code-generate CabanaPIC
WORKDIR /home/myuser/OP-PIC/apps/app_cabanapic
RUN . /home/myuser/OP-PIC/opp_translator/opp_venv/bin/activate && \
  python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths cabana.cpp

# compile lib for CabanaPIC
WORKDIR /home/myuser/OP-PIC/opp_lib
RUN make PETSC=0 seq -j4
RUN make PETSC=0 mpi -j4
RUN make PETSC=0 omp -j4

# compile CabanaPIC application
WORKDIR /home/myuser/OP-PIC/apps/app_cabanapic
RUN make PETSC=0 seq -j4
RUN make PETSC=0 mpi -j4
RUN make PETSC=0 omp -j4

# Run CabanaPIC for 100 iterations with smaller configurations ----------------

# sequential run with one CPU core
RUN bin/seq configs/cabana.param > run_cabana_seq.log

# MPI run with 4 CPU cores
RUN mpirun -np 4 bin/mpi configs/cabana.param > run_cabana_mpi.log

# OpenMP run with 4 CPU cores
ENV OMP_NUM_THREADS=4
ENV OMP_PROC_BIND=close
ENV OMP_PLACES=cores
RUN bin/omp configs/cabana.param > run_cabana_omp.log

# -----------------------------------------------------------------------------
WORKDIR /home/myuser/

# There will be run log files generated for both applications using sequential, 
# openmp and mpi backends.
# The Mini-FEM-PIC log files can be found in /home/myuser/OP-PIC/app_fempic 
# with names run_fempic_seq.log, run_fempic_omp.log, run_fempic_mpi.log and
# run_fempic_mpi_hdf5.log
# CabanaPIC log files can be found in /home/myuser/OP-PIC/app_cabanapic 
# with names run_cabana_seq.log, run_cabana_omp.log and run_cabana_mpi.log

# Check the created docker image using below.
# docker images
# docker run -it <IMAGE ID> /bin/bash


# INSTALLATION DEPLOYMENT PROCESS ***************************************************************

# Summary ------------------------------------------------------------------------------------------------

# The OP-PIC source code can be found in the GitHub repository: https://github.com/OP-DSL/OP-PIC and the 
# archived artifacts will contain a copy of the code used for the simulations.

# Read the Docs OP-PIC documentation be found in https://op-dsl.github.io/OP-PIC/

# The current OP-PIC eDSL supports generating code targeting multi-core CPUs with OpenMP threading, many-core 
# GPUs with CUDA offloading, and distributed memory cluster variants of these using MPI. Particularly 
# parallelizations are, sequential, OpenMP, MPI, CUDA, HIP and distributed memory CUDA+MPI, HIP+MPI.

# OP-PIC can also run on desktop machines even on a single CPU core.

# Steps A-E below explain the setup, compile and execution of OP-PIC for the two applicatons, Mini-FEM-PIC
# and CabanaPIC.

# Step A : Setup the environment -------------------------------------------------------------------------

#     1. The system should contain a functional C++ compiler [g++, icpc, CC].
#     2. For mpi builds, the system should contain a functional mpi compiler [OpenMPI, oneapi, cray-mpich].
#     3. For nvidia GPU builds, the system should contain a functional CUDA compiler >= 10.2.89 [nvcc].
#     4. For AMD GPU builds, the system should contain a functional HIP compiler >= 5.6.0 [hipcc].
#     5. The system should also contain external Libraries PETSc and HDF5.
#     6. The system should also contain an installation of Python >= 3.8.
    
#     Note : The above can be user installed or module loaded.

#     7. Use the below exports depending on the target platform
        
#         export CC_COMPILER=<c++ compiler, eg. gcc, icpc, CC>
#         export MPI_COMPILER=<mpi compiler wrapper, eg. mpicxx, CC>
#         export HIP_COMPILER=<hip compiler, eg. hipcc, CC>
#         export NVCCFLAGS_ADD='-gencode arch=compute_70,code=sm_70'
#             V100: 'compute_70,code=sm_70'
#             P100: 'compute_60,code=sm_60'
#             H100: 'compute_90,code=sm_90'

#         If CC is used instead of hipcc, use
#             export HIPCCFLAGS_ADD="-x hip"
    
#     8. Export OP-PIC directores

#         export OPP=<path to OP-PIC folder>
#         export OPP_INSTALL_PATH=$OPP/opp_install
#         export OPP_TRANSLATOR=$OPP/opp_translator/translator
#         export OPP_PATH=$OPP/opp_lib

#     9. If external libs are user installed and not loaded by modules, export their paths

#         export PETSC_INSTALL_PATH=<petsc installed directory>
#         export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH

#         export HDF5_INSTALL_PATH=<hdf5 installed directory>
#         export PATH=$HDF5_INSTALL_PATH/bin:$PATH
#         export LD_LIBRARY_PATH=$HDF5_INSTALL_PATH/lib:$LD_LIBRARY_PATH
    
#     10. Activate the opp_translator python environment 

#         Note: If the python environment is not created already, skip to Step B

#         source $OPP/opp_translator/opp_venv/bin/activate

#     Note: Example source files for setting up the enviornment can be found in the 'source_files' folder.

# Step B : Setup OP-PIC translator python environment (One time process) ---------------------------------

#     Note: If step A.10 is done, skip to Step C

#     1. Change directory to OP-PIC/opp_translator.
#     2. Run setup_venv.sh. 
#         This will download the required libraries, including libclang and create the opp_venv python 
#         environment. It will also create a venv folder within 'OP-PIC/opp_translator' folder.
#     3. Activate the opp_translator python environment
#         source $OPP/opp_translator/opp_venv/bin/activate

# Step C : Compile the OP-PIC library --------------------------------------------------------------------

#     1. Change directory to OP-PIC/opp_lib.
#     2. Compile the required platform specific backend library using make commands.
#         make seq      -- CPU single core sequential lib
#         make mpi      -- CPU shared/distributed memory MPI lib
#         make omp      -- CPU shared memory OpenMP lib
#         make omp_mpi  -- CPU shared/distributed memory OpenMP+MPI lib
#         make cuda     -- NVidia single GPU lib 
#         make cuda_mpi -- NVidia multiple GPU distributed memory MPI lib 
#         make hip      -- AMD single GPU lib 
#         make hip_mpi  -- AMD multiple GPU distributed memory MPI lib 
        
#         PS. The default make commands provided above will link to PETSc library which is essential for 
#             Mini-FEM-PIC. However, for CabanaPIC, build the target library using 'make PETSC=0 <target>'

# Step D : Compile and run Mini-FEM-PIC ------------------------------------------------------------------

#     Note: For Mini-FEM-PIC building OP-PIC library with PETSc library is essential.

#     1. Change directory to OP-PIC/app_fempic.
#     2. Code-generate Mini-FEM-PIC

#         Both the regular and the HDF5 applications can be independently code generated using the below.

#         For regular OP-PIC application without HDF5, use 
#             python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths fempic.cpp

#         If HDF5 is required, use 
#             python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths fempic_hdf5.cpp

#         Once the code-generator is invoked, a fempic_opp.cpp or fempic_hdf5_opp.cpp file and 
#         seq, omp, mpi, cuda and hip folders (including opp_kernels.<cpp|cu> and a loop kernel header  
#         file per unique opp_par_loop or opp_particle_move loop) will get generated.

#         Note: This step can be skipped if pre-code-generated files are being used.
#             They are present in 'OP-PIC/app_fempic_cg' folder.

#     3. Compile Mini-FEM-PIC for the required platform using the make command

#         For regular OP-PIC application without HDF5, use
#             make seq
#             make omp
#             make mpi
#             make cuda
#             make cuda_mpi
#             make hip
#             make hip_mpi

#         If HDF5 is required, use
#             make mpi_hdf5
#             make cuda_mpi_hdf5
#             make hip_mpi_hdf5

#         Note: Upon successful compilation, app-binaries can be found in OP-PIC/app_fempic/bin folder.
    
#     4. Setup the configuration file

#         An example configuration file 'coarse.param' can be found in 'OP-PIC/app_fempic/configs' folder.

#         The Mini-FEM-PIC mesh files can be found in the mesh_files folder of the artifact zip.
        
#         Two non-hdf5 mesh file folders (48000, 96000) are provided for non-hdf5 Mini-FEM-PIC runs.
#         To run with 48000 cell mesh, change the coarse.param file by updating,
#             i.   STRING global_mesh  = <mesh_folder_directory>/48000/mesh.dat
#             ii.  STRING inlet_mesh   = <mesh_folder_directory>/48000/inlet.dat
#             iii. STRING wall_mesh    = <mesh_folder_directory>/48000/wall.dat

#         For hdf5 Mini-FEM-PIC binaries, HDF5 (box_<cell_count>.hdf5) mesh files can be used. 
#         To run with 48000 cell mesh using HDF5 files, change the coarse.param file by updating, 
#             STRING hdf_filename = <hdf5_mesh_folder_directory>/box_48000.hdf5

#         Similarly, change for other mesh files.

#         In addition, update the random number file path using,
#             STRING rand_file    = <mesh_folder_directory>/random_100k.dat

#         To run Mini-FEM-PIC with ~70 million particles, change plasma_den to 1.0e18

#     5. Run the application

#         For regular OP-PIC application without HDF5, use
#             bin/seq configs/coarse.param
#             bin/omp configs/coarse.param
#             bin/cuda configs/coarse.param
#             bin/hip configs/coarse.param
#             mpirun -np <num_ranks> bin/mpi configs/coarse.param
#             mpirun -np <num_ranks> bin/cuda_mpi configs/coarse.param
#             mpirun -np <num_ranks> bin/hip_mpi configs/coarse.param

#         If HDF5 is required, use
#             mpirun -np <num_ranks> bin/mpi_hdf5 configs/coarse.param
#             mpirun -np <num_ranks> bin/cuda_mpi_hdf5 configs/coarse.param
#             mpirun -np <num_ranks> bin/hip_mpi_hdf5 configs/coarse.param

#         Note: srun can also be used for execution

# Step E : Compile and run CabanaPIC ---------------------------------------------------------------------

#     Note: For CabanaPIC the OP-PIC lib should be built without PETSc (using PETSC=0 in the make command).

#     1. Change directory to OP-PIC/app_cabanapic.
#     2. For CabanaPIC code-generation use the below command

#             python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths cabana.cpp

#         Once the code-generator is invoked, a cabana_opp.cpp file and 
#         seq, omp, mpi, cuda and hip folders (including opp_kernels.<cpp|cu> and a loop kernel header  
#         file per unique opp_par_loop or opp_particle_move loop) will get generated.

#         Note: This step can be skipped if pre-code-generated files are being used.
#             They are present in 'OP-PIC/app_cabanapic_cg' folder.

#     3. Compile CabanaPIC for the required platform using the make command

#             make seq
#             make omp
#             make mpi
#             make cuda
#             make cuda_mpi
#             make hip
#             make hip_mpi

#         Note: Upon successful compilation, app-binaries can be found in OP-PIC/app_cabanapic/bin folder.
    
#     4. Setup the configuration file

#         An example configuration file 'cabana.param' can be found in 'OP-PIC/app_fempic/configs' folder.

#         CabanaPIC does not require any mesh files and the mesh is created at runtime.
        
#         To run with 96000 cell mesh with 72 million particles, change the cabana.param file by updating,
#             INT nz = 60
#         For 144 million particles change,
#             INT num_part_per_cell = 1500
#         And change the main loop iteration count as required using 'num_steps'.

#     5. Run the application using the below commands,

#             bin/seq configs/cabana.param
#             bin/omp configs/cabana.param
#             bin/cuda configs/cabana.param
#             bin/hip configs/cabana.param
#             mpirun -np <num_ranks> bin/mpi configs/cabana.param
#             mpirun -np <num_ranks> bin/cuda_mpi configs/cabana.param
#             mpirun -np <num_ranks> bin/hip_mpi configs/cabana.param

#         Note: srun can also be used for execution

# ********************************************************************************************************
