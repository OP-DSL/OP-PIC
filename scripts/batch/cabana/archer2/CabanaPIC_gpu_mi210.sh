#!/bin/bash

#SBATCH --job-name=cabG_mi210
#SBATCH --gpus=1
#SBATCH --time=00:10:00
#SBATCH --account=e723-neptune 
#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc
#SBATCH --nodes=1
#SBATCH --exclusive

module purge
module load PrgEnv-amd
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan

export OMP_NUM_THREADS=1
export OMP_PROC_BIND=close
export OMP_PLACES=cores

export MPICH_GPU_SUPPORT_ENABLED=1

srun --ntasks=1 rocm-smi

srun --ntasks=1 --cpus-per-task=1 /work/e723/e723/csrcnj/phd/OP-PIC/app_cabanapic_cg/bin/hip_mpi cabana.param

srun --ntasks=1 rocm-smi
