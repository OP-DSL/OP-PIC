# OP-PIC NESO-Advection (Already code generated for reference)

**Please follow the instructions provided in the main readme file first.**

**This README file contains only additional information regarding the application.**

##
**This folder contain both user written code and OP-PIC code generated code.**

All the user written code will have the below comment;
```cpp
// *********************************************
// USER WRITTEN CODE                            
// *********************************************
```

The code-generator has added the below comment to all the generated code; 
```cpp
// *********************************************
// AUTO GENERATED CODE                          
// *********************************************
```

If code-generator is invoked in this folder, the available generated code will be replaced with the newly generated code (May generate the same if `advec.cpp` or `kernels.h` is not changed)

##
NESO Advection is a simple 2D demonstrator unstructured mesh advection benchmark code that can be used as a particle mover of a PIC application (https://github.com/will-saunders-ukaea/hybrid_move_benchmark/).

It is based on quadrilateral mesh cells and particles are periodically move through the mesh by updating its position using with `s=ut`, without any field influence. 
Overall NESO_Advection has 3 DOFs per particle.

Here, we have implemented the application with OP-PIC, solving the same physics as the original.

The main function is included in `advec.cpp` containing OP-PIC API calls. 
It contains one mesh cell set, and a particle set. 
NESO-Advection also contain two maps. 
They are, `cell to cell map` and `particle to cell map`.

## Structure
 * `advec.cpp` : The main file containing OP-PIC API calls. 
 * `kernels.h` : The user written elemental kernel functions.
 * `advec_defs.h` : The defines used in the simulation
 * `advec_misc.h` : Miscellaneous functions including particle distribution initialization and the structure used to hold mesh data till OP-PIC is initialized..

## Code Generation
Code generation can be done by invoking the below command.
```bash
python3 $OPP_TRANSLATOR -I$OPP_PATH/include/ --file_paths advec.cpp
```

Once the code-generator is invoked, a `advec_opp.cpp` file and `seq`, `omp`, `mpi`, `cuda`, `hip` and `sycl` folders, including `opp_kernels.<cpp|cu>` and a loop kernel header file per unique `opp_par_loop` or `opp_particle_move` loop will get generated.

## Compile
Once the platform specific target files are generated, use the provided `MakeFile` to compile the application.
```bash
 make seq
 make mpi
 make omp
 make omp_mpi
 make cuda
 make cuda_mpi
 make hip
 make hip_mpi
 make sycl
 make sycl_mpi
```

### Configuration
An example configuration file is provided in `configs` folder.

This file can be used to change the application configurations such as number of steps in the main iterating loop (`num_steps`), domain mesh size (`nx, ny`) and other parameters included in the original neso_advection application. 

In addition, the config file include OP-PIC DSL simulation parameters, such as gpu threads per block (`opp_threads_per_block`), particle move finalizing mechanism (`opp_fill`=[`HoleFill_All`, `Sort_All`, `Shuffle_All`, `Sort_Periodic`, `Shuffle_Periodic`]).

Generated direct-hop overlay meshes can be found in [OP-PIC-Meshes](https://github.com/ZamanLantra/OP-PIC_meshes)

# Run
To run the application, below commands can be used.
```bash
    bin/seq configs/advec.param
    bin/omp configs/advec.param
    bin/cuda configs/advec.param
    bin/hip configs/advec.param
    bin/sycl configs/advec.param
    mpirun -np <num_ranks> bin/mpi configs/advec.param
    mpirun -np <num_ranks> bin/cuda_mpi configs/advec.param
    mpirun -np <num_ranks> bin/hip_mpi configs/advec.param
    mpirun -np <num_ranks> bin/sycl_mpi configs/advec.param
```

In addition, `srun` can be used for execution.