#!/bin/bash

# Get the local MPI rank from environment variable
local_rank=${PMI_RANK:-0}
# local_rank=${OMPI_COMM_WORLD_LOCAL_RANK:-0}

# Calculate GPU index by taking modulo of available GPUs on the node
gpu_index=$((local_rank % 2))

# Set SYCL device based on calculated index
export ONEAPI_DEVICE_SELECTOR="level_zero:${gpu_index}"

echo "intra-rank:"$local_rank " ONEAPI_DEVICE_SELECTOR:"$ONEAPI_DEVICE_SELECTOR 

# Execute the original command with all arguments
exec "$@"