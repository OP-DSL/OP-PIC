# OP-PIC CabanaPIC

**This folder contain fully user written code, ready for code generation.**

##
CabanaPIC is a 3D electromagnetic, two stream PIC code, where particles move in a duct (cuboid) with cuboid cells.
It is implemented with periodic boundaries and has 9 DOFs per cell and 7 DOFs per particle. 
It was originally developed as part of the Exascale Computing Project, Co-design center for Particle Applications (CoPA) (https://github.com/ECP-copa/CabanaPIC) using the Cabana library, based on Kokkos as a structured-mesh PIC application.

Here, we have implemented the application with OP-PIC, using unstructured-mesh mappings solving the same physics as the original.

## Structure
 * `cabana.cpp` : The main file containing OP-PIC API calls. 
 * `kernels.h` : The user written elemental kernel functions.
 * `cabana_defs.h` : The defines and the structure used to hold data till OP-PIC is initialized.
 * `cabana_misc.h` : Miscellaneous functions including particle initialization (enrichment prior simulation).
 * `cabana_misc_mesh_color.h` : This include a user defined partitioning scheme (Only used in MPI builds).
 * `cabana_misc_mesh_loader.h` : Mesh creating code, used prior OP-PIC DSL initializing. 

## Code Generation
Code generation can be done by invoking the below command.

`python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths cabana.cpp`

Once the code-generator is invoked, a `cabana_opp.cpp` file and `seq`, `omp`, `mpi`, `cuda`, `hip` and `sycl` folders, including `opp_kernels.<cpp|cu>` and a loop kernel header file per unique `opp_par_loop` or `opp_particle_move` loop will get generated.

## Compile

Once the platform specific target files are generated, use the provided `MakeFile` to compile the application.
 * `make seq`
 * `make mpi`
 * `make omp`
 * `make omp_mpi`
 * `make cuda`
 * `make cuda_mpi`
 * `make hip`
 * `make hip_mpi`
 * `make sycl`
 * `make sycl_mpi`

## Configuration
An example configuration file is provided in `OP-PIC/app_cabanapic/configs` folder.

This file can be used to change the application configurations such as number of steps in the main iterating loop (`num_steps`), domain mesh size (`nx, ny, nz`) and other parameters included in the original CabanaPIC application. 

In addition, the config file include OP-PIC DSL simulation parameters, such as gpu threads per block (`opp_threads_per_block`), particle move finalizing mechanism (`opp_fill`=[`HoleFill_All`, `Sort_All`, `Shuffle_All`, `Sort_Periodic`, `Shuffle_Periodic`]).

# Run
To run the application, below commands can be used.
 * `bin/seq configs/cabana.param`
 * `bin/omp configs/cabana.param`
 * `bin/cuda configs/cabana.param`
 * `bin/hip configs/cabana.param`
 * `mpirun -np <num_ranks> bin/mpi configs/cabana.param`
 * `mpirun -np <num_ranks> bin/cuda_mpi configs/cabana.param`
 * `mpirun -np <num_ranks> bin/hip_mpi configs/cabana.param`

In addition, `srun` can be used for execution.