#!/bin/bash

#SBATCH --job-name=comp_cab   
#SBATCH --partition=mi300x 
#SBATCH --nodes=1          
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=12
#SBATCH --time=0-00:05:00  
#SBATCH --account=do018
##SBATCH --exclusive

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

cd /cosma/home/do018/dc-lant1/phd/OP-PIC/opp_lib;
make hip_mpi -j

cd /cosma/home/do018/dc-lant1/phd/OP-PIC/app_cabanapic_cg;
make DEBUG_LOG=0 hip_mpi -j

# echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"

# srun --nodes=1 --ntasks=1 --cpus-per-task=1 --partition=mi300x --job-name=comp_fem --time=1:00:00 --account=do018 --pty bash
# srun --nodes=1 --ntasks=12 --cpus-per-task=1 --partition=mi300x --job-name=comp_fem --time=1:00:00 --account=do018 --pty bash


# srun /cosma/home/do018/dc-lant1/phd/OP-PIC/app_cabanapic_cg/bin/hip_mpi /cosma/home/do018/dc-lant1/phd/OP-PIC/app_cabanapic_cg/configs/cabana.param