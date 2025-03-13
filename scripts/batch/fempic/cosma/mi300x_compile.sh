#!/bin/bash

#SBATCH --job-name=comp_fem   
#SBATCH --partition=mi300x 
#SBATCH --nodes=1          
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0-00:03:00  
#SBATCH --account=do018

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

cd $HOME/phd/OP-PIC/opp_lib;
make PETSC=1 hip_mpi -j

cd $HOME/phd/OP-PIC/app_fempic_cg;
make PETSC=1 HDF5=1 DEBUG_LOG=0 hip_mpi_hdf5 -j

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"

# srun --nodes=1 --ntasks=1 --cpus-per-task=1 --partition=mi300x --job-name=comp_fem --time=1:00:00 --account=do018 --pty bash
# srun --nodes=1 --ntasks=96 --cpus-per-task=1 --partition=mi300x --job-name=comp_fem --time=2:00:00 --account=do018 --exclusive --pty bash


# lscpu

# rocm-smi

# rocminfo



# lspci -vvv | grep -i 'MI300'