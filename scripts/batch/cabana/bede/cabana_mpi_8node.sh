#!/bin/bash

#SBATCH --job-name=cab8N
#SBATCH --time=01:00:00
#SBATCH --nodes=8

#SBATCH --account=bdyrk17            
#SBATCH --partition=gpu
#SBATCH --gres=gpu:4

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI8_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder" $runFolder

binary="/users/csrcnl/phd/OP-PIC/cabana_mpi/bin/cuda_mpi"
echo "Using Binary" $binary

module load gcc/12.2 openmpi/4.0.5 cuda/12.0.1

# module purge
# module load intel/2022a

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close
# export OMP_NUM_THREADS=1

num_nodes=$SLURM_JOB_NUM_NODES

configFile="cabana.param"
file=$PWD/$configFile

for config in 750 1500 2000; do
    for run in 1 2 3; do

        echo $file $config

        # ****************************************
        echo "Running CUDA-MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nz = 30/INT nz = 1920/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        sed -i "s/STRING cluster = pencil/STRING cluster = block/" ${currentfilename}
        bede-mpirun --bede-par 1ppg $binary $currentfilename  | tee $folder/log_60_N${num_nodes}_Block_D${config}_R${run}.log;
    done
done

for config in 1500 2000 3000; do
    for run in 1 2 3; do

        echo $file $config

        # ****************************************
        echo "Running CUDA-MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nz = 30/INT nz = 960/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        sed -i "s/STRING cluster = pencil/STRING cluster = block/" ${currentfilename}
        bede-mpirun --bede-par 1ppg $binary $currentfilename  | tee $folder/log_30_N${num_nodes}_Block_D${config}_R${run}.log;
    done
done

for config in 750 1500 2000; do
    for run in 1 2 3; do

        echo $file $config

        # ****************************************
        echo "Running CUDA-MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nz = 30/INT nz = 1920/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        sed -i "s/STRING cluster = pencil/STRING cluster = cart/" ${currentfilename}
        bede-mpirun --bede-par 1ppg $binary $currentfilename  | tee $folder/log_N${num_nodes}_Cart_D${config}_R${run}.log;
    done
done

for config in 750 1500 2000; do
    for run in 1 2 3; do

        echo $file $config

        # ****************************************
        echo "Running CUDA-MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nz = 30/INT nz = 1920/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        bede-mpirun --bede-par 1ppg $binary $currentfilename  | tee $folder/log_N${num_nodes}_Pencil_D${config}_R${run}.log;
    done
done


echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
