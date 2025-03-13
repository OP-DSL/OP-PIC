#!/bin/bash

# Comma separated GPUs that can be accessed, first GPU should be 0
export CUDA_VISIBLE_DEVICES=0,1
# export HIP_VISIBLE_DEVICES=3,4,5,6

# For SYCL 
# local_rank=${PMI_RANK:-0}                  # For intel-mpi
# local_rank=${OMPI_COMM_WORLD_LOCAL_RANK:-0 # For openmpi
# gpu_index=$((local_rank % 2))
# export ONEAPI_DEVICE_SELECTOR="level_zero:${gpu_index}"

# Run example: mpirun -np 2 ./gpu_wrapper.sh bin/cuda_mpi configs/coarse.param
exec "$@"