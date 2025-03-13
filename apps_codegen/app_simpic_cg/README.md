# OP-PIC SimPIC (Already code generated for reference)

**This folder contain both user written code and OP-PIC code generated code.**

All the user written code will have the below comment;

`// *********************************************`<br>
`// USER WRITTEN CODE                            `<br>
`// *********************************************`

The code-generator has added the below comment to all the generated code; 

`// *********************************************`<br>
`// AUTO GENERATED CODE                          `<br>
`// *********************************************`

If code-generator is invoked in this folder, the available generated code will be replaced with the newly generated code (May generate the same if `simpic.cpp` or `kernels.h` is not changed)

##
SimPIC is a 1D electrostatic, PIC code (https://bitbucket.org/lecad-peg/simpic/src/master/).

It has 3 DOFs per cell and 3 DOFs per particle.

Here, we have implemented the application with OP-PIC, using unstructured-mesh mappings solving the same physics as the original.

**The MPI based user code is currently being developed. i.e. Distributed memory mesh and particle distribution, MPI partitioning scheme, make the tri-diagonal solver to use the halo updated data need to be implemented.**

**For now, `seq`, `omp`, `cuda` and `hip` shared memory versions are functional.**

*The original SimPIC MPI code does not seem to give the exact same answer compared to its serial versions.*

* Only the particles within the current MPI rank provide contributions to the mesh. If communicated, the contribution to the received rank is not computed!
https://bitbucket.org/lecad-peg/simpic/src/ab6a92dea645ee39747aff214884a51a14802f76/mpi/simpic.cpp#lines-96
* When communicating fields, left most and right most `phiarray` values are sent to all the MPI ranks (not only to the adjoining neighbour rank) using `MPI_Allgather` and once the communication is done, the edge `phiarray` values of all MPI ranks are added to all the `phiarray` values of the current rank.
https://bitbucket.org/lecad-peg/simpic/src/ab6a92dea645ee39747aff214884a51a14802f76/mpi/fields.cpp#lines-130

*Both these issues will make the shared memory serial and distributed memory runs deviate even if we use the same number of particles per cell and the same total number of particles and mesh elements per simulation*

## TODO
If MPI versions of Simpic needs to be developed with OP-PIC, the above issues should be addressed.

In addition, the mesh/particles should be distributed, an appropriate partitioning routine should be called and the field solver needs to be re-implemented, perhaps consider using PETSc.

## Structure
 * `simpic.cpp` : The main file containing OP-PIC API calls. 
 * `kernels.h` : The user written elemental kernel functions.
 * `simpic_defs.h` : The defines, the structure used to hold data till OP-PIC is initialized and the mesh/particle loading functions.
 * `funcs_from_simpic.cpp` : Original matrix tri-diagonal solver used in the simulation.

## Code Generation
Code generation can be done by invoking the below command.

`python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths simpic.cpp`

Once the code-generator is invoked, a `simpic_opp.cpp` file and `seq`, `omp`, `mpi`, `cuda` and `hip` folders, including `opp_kernels.<cpp|cu>` and a loop kernel header file per unique `opp_par_loop` or `opp_particle_move` loop will get generated.

## Compile
Once the platform specific target files are generated, use the provided `MakeFile` to compile the application.
 * `make seq`
 * `make omp`
 * ~~`make mpi`~~
 * `make cuda`
 * ~~`make cuda_mpi`~~
 * `make hip`
 * ~~`make hip_mpi`~~

## Configuration
An example configuration file is provided in `OP-PIC/app_simpic_cg/configs` folder.

This file can be used to change the application configurations such as number of steps in the main iterating loop (`num_steps`), number of cells per process (`ncpp`), number of particles per process (`ppc`) and other parameters included in the original SimPIC application. 

In addition, the config file include OP-PIC DSL simulation parameters, such as gpu threads per block (`opp_threads_per_block`), particle move finalizing mechanism (`opp_fill`=[`HoleFill_All`, `Sort_All`, `Shuffle_All`, `Sort_Periodic`, `Shuffle_Periodic`]).

# Run
To run the application, below commands can be used.
 * `bin/seq configs/system.param`
 * `bin/omp configs/system.param`
 * `bin/cuda configs/system.param`
 * `bin/hip configs/system.param`
 * ~~`mpirun -np <num_ranks> bin/mpi configs/system.param`~~
 * ~~`mpirun -np <num_ranks> bin/cuda_mpi configs/system.param`~~
 * ~~`mpirun -np <num_ranks> bin/hip_mpi configs/system.param`~~

In addition, `srun` can be used for execution.