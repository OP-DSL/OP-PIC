# OP-PIC
OP-PIC is a high-level embedded domain specific language for writing **unstructured mesh particle in cell** simulations with automatic parallelization on multi-core and many-core architectures. The parallelizations include, sequential(`seq`), OpenMP(`omp`), MPI(`mpi`), CUDA(`cuda`), HIP(`hip`), SYCL(`sycl`) and their CUDA+MPI(`cuda_mpi`), HIP+MPI(`hip_mpi`), SYCL+MPI(`sycl_mpi`) distributed memory parallelizations, The API is embedded in C/C++.

This repository contains the implementation of the code translation tools and run-time support libraries, and is structured as follows:
 * `opp_lib`: The C/C++ OP-PIC run-time libraries.
 * `opp_translator`: The Python code translator for C/C++.
 * `app_<app_name>`: Example applications that demonstrate use of the API ready to be code-generated.
 * `app_<app_name>_cg`: Copy of the example applications that are already being code generated. 
 * `scripts/source`: Some example source files used during implementation
 * `scripts/batch`: Some example slurm batch files used during implementation
 * `app_handcoded`: Some hand-coded application code, written prior code generation

## Documentation
Documentation is available on [Read the Docs](https://zamanlantra.github.io/git_docs_test/index.html).

## Dependencies
OP-PIC has a variety of toolchain dependencies that you will likely be able to obtain from your package manager or programming environment:
 * GNU Make > 4.2
 * A C/C++ compiler: Currently supported compilers are Intel, GCC, Clang, Cray, and NVHPC.
 * Python >= 3.8
 * (Optional) An MPI implementation: Any implementation with the `mpicc` and `mpicxx` wrappers is supported.
 * (Optional) NVIDIA CUDA >= 10.2.89
 * (Optional) AMD HIP-ROCM >= 5.6.0
 * (Optional) Intel-oneAPI DPCPP compiler (including oneDPL)

In addition, there are a few optional library dependencies that you will likely have to build manually, although some package managers or programming environments may be able to provide appropriate versions:
 * (Optional) [HDF5](https://www.hdfgroup.org/solutions/hdf5/): Used for HDF5 I/O. You may build with and without `--enable-parallel` (depending on if you need MPI), and then specify via the environment variable `HDF5_INSTALL_PATH`.
* (Optional) [PETSc](https://petsc.org/release/install/download/): Used for linear/nonlinear sparse matrix solvers. You may build with and without `--with-mpi=1` (depending on if you need MPI), and then specify via the environment variable `PETSC_INSTALL_PATH`.
 * (Optional) [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview): Can be used for MPI mesh partitioning. Build *with* 32-bit indicies (`-DIDXSIZE32`) and *without* `-DSCOTCH_PTHREAD`. Specify ParMETIS via the environment variable `PARMETIS_INSTALL_PATH`.

## Compile OP-PIC library bindings
 * Source the compiler and other environment variables (some example source files can be found at `OP-PIC/scripts/source` directory).
 * Change directory to `OP-PIC/opp_lib`.
 * Compile the required platform specific backend library using make commands (e.g. `make cuda_mpi`).

## Setup OP-PIC translator python environment (One time process)
 * Change directory to `OP-PIC/opp_translator`.
 * Run `setup_venv.sh` (this will download the required libraries, including libclang and create the `opp_venv` python environment).

## Compile application
 * Change directory to the required demonstrator application folder (the applications ready to code-generate will have the name in  format `app_<app_name>`).
 * Run `python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths <application_cpp_file>`. 
 * Compile the required application version using the make command. (`make` followed by the required parallelization). 
 
 For example, `cd app_cabanapic; python3 $OPP_TRANSLATOR -v -I$OPP_PATH/include/ --file_paths cabanapic.cpp; make cuda_mpi`.
 A detail explanation can be found in the readme file of `opp_translator` folder and the required application folder.
 
 In addition, `app_<app_name>_cg` will additionally include code-generated files, ready to compile directly using make commands.

## Run the compiled application
 * To run the application, use `bin/<parallelization> <config_file>`. For example, `bin/hip configs/cabana.params`. For distributed memory MPI builds use `mpirun`, or `srun` for slurm runs. Some example slurm scripts can be found at `batch/app_<app_name>` folder.

## License 
OP-PIC is released as an open-source project under the BSD 3-Clause License. See the file called [LICENSE](https://github.com/OP-DSL/OP-PIC/blob/main/LICENSE) for more information.

## Citing
To cite OP-PIC, please reference the following paper:

[Zaman Lantra, Steven A. Wright, and Gihan R. Mudalige. 2024. OP-PIC - an Unstructured-Mesh Particle-in-Cell DSL for Developing Nuclear Fusion Simulations. In Proceedings of the 53rd International Conference on Parallel Processing (ICPP '24). Association for Computing Machinery, New York, NY, USA, 294–304. https://doi.org/10.1145/3673038.3673130](https://dl.acm.org/doi/10.1145/3673038.3673130)

```
@inproceedings{10.1145/3673038.3673130,
author = {Lantra, Zaman and Wright, Steven A. and Mudalige, Gihan R.},
title = {OP-PIC - an Unstructured-Mesh Particle-in-Cell DSL for Developing Nuclear Fusion Simulations},
year = {2024},
isbn = {9798400717932},
publisher = {Association for Computing Machinery},
address = {New York, NY, USA},
url = {https://doi.org/10.1145/3673038.3673130},
doi = {10.1145/3673038.3673130},
pages = {294–304},
numpages = {11},
keywords = {DSL, OP-PIC, PIC, Particle-In-Cell, Unstructured-mesh},
location = {Gotland, Sweden},
series = {ICPP '24}
}
```
