# OP-PIC
OP-PIC is a high-level embedded domain specific language for writing **unstructured mesh particle in cell** simulations with automatic parallelization on multi-core and many-core architectures. The parallelizations include, sequential(`seq`), OpenMP(`omp`), MPI(`mpi`), CUDA(`cuda`), HIP(`hip`), SYCL(`sycl`) and their CUDA+MPI(`cuda_mpi`), HIP+MPI(`hip_mpi`), SYCL+MPI(`sycl_mpi`), OMP+MPI(`omp_mpi`) distributed memory parallelizations, The API is embedded in C/C++.

This repository contains the implementation of the code translation tools and run-time support libraries, and is structured as follows:
 * `opp_lib`: The C/C++ OP-PIC run-time libraries.
 * `opp_translator`: The Python code translator for C/C++.
 * `apps`: Example applications that demonstrate use of the API ready to be code-generated.
 * `apps_codegen`: Copy of the example applications that are already being code generated. 
 * `apps_handcoded`: Some hand-coded application code, written prior code generation
 * `scripts/source`: Some example source files used during implementation
 * `scripts/batch`: Some example slurm batch files used during implementation

## Documentation
Documentation is available on [Read the Docs](https://op-dsl.github.io/OP-PIC/).

## Dependencies
OP-PIC has a variety of toolchain dependencies that you will likely be able to obtain from your package manager or programming environment:
 * GNU Make > 4.2
 * A C/C++ compiler: Tested compilers are Intel, GCC, Clang, Cray, and NVHPC.
 * Python >= 3.8
 * (Optional) An MPI implementation: Any implementation with the `mpicc` and `mpicxx` wrappers are supported.
 * (Optional) NVIDIA CUDA >= 10.2.89
 * (Optional) AMD HIP-ROCM >= 5.6.0
 * (Optional) Intel-oneAPI DPCPP compiler (including oneDPL and DPCT) >= 2024.2.0

In addition, there are a few optional library dependencies that you will likely have to build manually, although some package managers or programming environments may be able to provide appropriate versions:
 * (Optional) [HDF5](https://www.hdfgroup.org/solutions/hdf5/): Used for HDF5 I/O. You may build with and without `--enable-parallel` (depending on if you need MPI), and then specify via the environment variable `HDF5_INSTALL_PATH`.
* (Optional) [PETSc](https://petsc.org/release/install/download/): Used for linear/nonlinear sparse matrix solvers. You may build with and without `--with-mpi=1` (depending on if you need MPI), and then specify via the environment variable `PETSC_INSTALL_PATH`.
 * (Optional) [ParMETIS](http://glaros.dtc.umn.edu/gkhome/metis/parmetis/overview): Can be used for MPI mesh partitioning. Build *with* 32-bit indicies (`-DIDXSIZE32`) and *without* `-DSCOTCH_PTHREAD`. Specify ParMETIS via the environment variable `PARMETIS_INSTALL_PATH`.

## Environment Variables Required
```bash
export OPP=<path_to_OP-PIC_repo>
export OPP_INSTALL_PATH=$OPP/opp_install
export OPP_TRANSLATOR=$OPP/opp_translator/translator
export OPP_PATH=$OPP/opp_lib

export CC_COMPILER=<C++_Compiler>
export MPI_COMPILER=<MPI_Compiler>
export HIP_COMPILER=<HIP_Compiler>
```

Optionally, you may require,
```bash
export MPI_INSTALL_PATH=<MPI_install_path>
export CUDA_INSTALL_PATH=<CUDA_toolkit_installed_path>
export ROCM_INSTALL_DIR=<ROCm_installed_path>
export PETSC_INSTALL_PATH=<PETSc_install_path>
export HDF5_INSTALL_PATH=<HDF5_install_path>

export CPPFLAGS_ADD=<Additional_CPP_flags>
export NVCCFLAGS_ADD=<Additional_CUDA_flags>
export HIPCCFLAGS_ADD=<Additional_HIP_flags>
export SYCLFLAGS_ADD=<Additional_SYCL_flags>
```

## OP-PIC parallelization names
 `seq`, `omp`, `mpi`, `cuda`, `hip`, `sycl`, `omp_mpi`, `cuda_mpi`, `hip_mpi`, `sycl_mpi`

## Compile OP-PIC library bindings
 * Source the compiler and other environment variables (some example source files can be found at `OP-PIC/scripts/source` directory).
 * Change directory to `OP-PIC/opp_lib`.
 * Compile the required platform specific backend library using make commands (e.g. `make cuda_mpi`).

## Setup OP-PIC translator python environment (One time process)
 * Change directory to `OP-PIC/opp_translator`.
 * Run `setup_venv.sh` (this will download the required libraries, including libclang and create the `opp_venv` python environment).

## Compile application
 * Change directory to the required demonstrator application folder (the applications ready to code-generate will have the name in  format `app_<app_name>`).
 * Run `python3 $OPP_TRANSLATOR -I$OPP_PATH/include/ --file_paths <application_cpp_file>`. 
 * Compile the required application version using the make command. (`make` followed by the required parallelization). 
 
 For example, `cd app_cabanapic; python3 $OPP_TRANSLATOR -I$OPP_PATH/include/ --file_paths cabanapic.cpp; make cuda_mpi`.
 A detail explanation can be found in the readme file of `opp_translator` folder and the required application folder.
 
 In addition, `apps_codegen` will additionally include code-generated applications, ready to compile directly using make commands.

## Run the compiled application
 * To run the application, use `bin/<parallelization> <config_file>`. For example, `bin/hip configs/cabana.params`. For distributed memory MPI builds use `mpirun`, or `srun` for slurm runs. Some example slurm scripts can be found at `batch/<app_name>` folder.

## License 
OP-PIC is released as an open-source project under the BSD 3-Clause License. See the file called [LICENSE](https://github.com/OP-DSL/OP-PIC/blob/main/LICENSE) for more information.