# OP-PIC FemPIC (Already code generated for reference)

**This folder contain both user written code and OP-PIC code generated code.**

All the user written code will have the below comment;

`// *********************************************`<br> 
`// USER WRITTEN CODE                            `<br> 
`// *********************************************`

The code-generator has added the below comment to all the generated code; 

`// *********************************************`<br>
`// AUTO GENERATED CODE                          `<br>
`// *********************************************`

If code-generator is invoked in this folder, the available generated code will be replaced with the newly generated code (May generate the same if `fempic.cpp`|`fempic_hdf5.cpp` or `kernels.h` is not changed)

##
FemPIC is a sequential electrostatic 3D unstructured mesh finite element PIC code (https://github.com/ExCALIBUR-NEPTUNE/Documents/blob/main/reports/2057699/TN-03-3.pdf).

It is based on tetrahedral mesh cells, nodes and faces forming a duct. 
Faces on one end of the duct are designated as inlet faces and the outer wall is fixed at a higher potential to retain the ions within the duct. 
Charged particles are injected at a constant rate from the inlet faces of the duct (one-stream) at a fixed velocity, and the particles move through the duct under the influence of the electric field. 
The particles are removed when they leave the boundary face. Overall Mini-FEM-PIC has 1 degree of freedom (DOF) per cell, 2 DOFs per node and 7 DOFs per particle.

Here, we have implemented the application with OP-PIC, using unstructured-mesh mappings solving the same physics as the original.

The main function is included in `fempic.cpp` containing OP-PIC API calls. 
It contains one mesh cell set, a node set, an inlet face set and a particle set. 
FemPIC also contain five maps. 
They are, `cell to cell map`, `cell to node map`, `inlet face to cell map`, a `inlet face to cell map` and `particle to cell map`.

Additionally, FemPIC is also implemented with OP-PIC HDF5 API calls (`opp_decl_set_hdf5`, `opp_decl_particle_set_hdf5`, `opp_decl_map_hdf5` and `opp_decl_dat_hdf5`) in `fempic_hdf5.cpp` file, that can be code generated and compiled as separate binaries. 

## Structure
 * `fempic.cpp` or `fempic_hdf5.cpp` : The main file containing OP-PIC API calls (without or with hdf5). 
 * `kernels.h` : The user written elemental kernel functions.
 * `field_solver.h` and `field_solver/cpu.cpp|cuda.cu|hip.cpp`: PETSc matrix based linear field solver.
 * `fempic_defs.h` : The defines and the structure used to hold data till OP-PIC is initialized.
 * `fempic_misc.h` : Miscellaneous functions including particle distribution initialization.
 * `fempic_misc_mesh_colour.h` : This include a user defined partitioning scheme (Only used in MPI builds).
 * `fempic_misc_mesh_loader.h` : Mesh creating code, used prior OP-PIC DSL initializing. 
 * `minifempic_funcs.h` : Some original MiniFemPIC code used to read the mesh files.

## Code Generation

Both the regular and the HDF5 applications can be independently code generated using the below.

For regular OP-PIC application without HDF5, use
`python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths fempic.cpp`

If HDF5 is required, invoke the below command.
`python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths fempic_hdf5.cpp`

Once the code-generator is invoked, a `fempic_opp.cpp` or `fempic_hdf5_opp.cpp` file and `seq`, `omp`, `mpi`, `cuda`, `hip` and `sycl` folders, including `opp_kernels.<cpp|cu>` and a loop kernel header file per unique `opp_par_loop` or `opp_particle_move` loop will get generated.

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
 * `make mpi_hdf5`
 * `make cuda_mpi_hdf5`
 * `make hip_mpi_hdf5`
 * `make sycl_mpi_hdf5`

### Configuration
An example configuration file is provided in `OP-PIC/app_fempic_cg/configs` folder.

This file can be used to change the application configurations such as number of steps in the main iterating loop (`num_steps`), domain mesh file (`global_mesh, inlet_mesh, wall_mesh` or `hdf_filename`) and other parameters such as `plasma_den`, `dt`. 

In addition, the config file include OP-PIC DSL simulation parameters, such as gpu threads per block (`opp_threads_per_block`), boolean to switch between segmented reductions and atomics (`use_reg_red`), a switch for direct-hop and multi-hop (`opp_global_move`=true is DH), particle move finalizing mechanism (`opp_fill`=[`HoleFill_All`, `Sort_All`, `Shuffle_All`, `Sort_Periodic`, `Shuffle_Periodic`]).

# Run
To run the application, below commands can be used.
 * `bin/seq configs/coarse.param`
 * `bin/omp configs/coarse.param`
 * `bin/cuda configs/coarse.param`
 * `bin/hip configs/coarse.param`
 * `bin/sycl configs/coarse.param`
 * `mpirun -np <num_ranks> bin/mpi configs/coarse.param`
 * `mpirun -np <num_ranks> bin/cuda_mpi configs/coarse.param`
 * `mpirun -np <num_ranks> bin/hip_mpi configs/coarse.param`
 * `mpirun -np <num_ranks> bin/sycl_mpi configs/coarse.param`

 * `mpirun -np <num_ranks> bin/mpi_hdf5 configs/coarse.param`
 * `mpirun -np <num_ranks> bin/cuda_mpi_hdf5 configs/coarse.param`
 * `mpirun -np <num_ranks> bin/hip_mpi_hdf5 configs/coarse.param`
 * `mpirun -np <num_ranks> bin/sycl_mpi_hdf5 configs/coarse.param`

In addition, `srun` can be used for execution.