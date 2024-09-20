# OP-PIC Library

The library include sequential(`seq`), OpenMP(`omp`), MPI(`mpi`), CUDA(`cuda`), HIP(`hip`), SYCL(`sycl`) and their CUDA+MPI(`cuda_mpi`), HIP+MPI(`hip_mpi`), SYCL+MPI(`sycl_mpi`) distributed memory parallelizations. 

`PETSc` and `ParMETIS` is optional, and can be embedded to the library using flags `PETSC=1` and `PARMETIS=1` during compilation. 

`HDF5` is also optional, but it is not directly embedded to the library archive file, but can be integrated during application compilation. 
For an example, check the `MakeFile` of `OP-PIC/app_fempic` (`mpi_hdf5`, `cuda_mpi_hdf5` and `hip_mpi_hdf5`).

The `CUDA` and `HIP` GPU library will be using `THRUST` library.

## Compiling

Build OP-PIC platform specific OP-PIC library archive using,
 * `make seq`
 * `make omp`
 * `make mpi`
 * `make cuda`
 * `make cuda_mpi`
 * `make hip`
 * `make hip_mpi`
 * `make sycl`
 * `make sycl_mpi`

The above `make` commands will build the release version using `-O3` flags without debug information.

If only debug logs are required, build using `DEBUG_LOG=1` and, to compile with both logs and with `-g -O0`, use `DEBUG=1` in the make command.

The compiled archives will be stored at `OP-PIC/opp_lib/lib` folder.