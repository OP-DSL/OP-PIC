#!/bin/bash

#SBATCH --job-name=fempic1
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=dp360
#SBATCH --qos=standard

#SBATCH --nodes=1
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

# partitions: gpu, gpu-a100-80, gpu-a100-40 | qos: standard, dev

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_SR0_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=/home/dp360/dp360/dc-lant1/phd/OP-PIC/app_fempic_cg/bin
binary=$binpath"/cuda_mpi_hdf5"

configFile="coarse.param"
file=/home/dp360/dp360/dc-lant1/phd/OP-PIC/scripts/batch/fempic/tursa/$configFile

echo "Creating running folder -> " $runFolder
echo "Using Binary -> " $binary
echo "Config file -> " $file
echo "********************************************************"
cd $binpath
gitbranch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
gitcommit=$(git log -n 1 $gitbranch)
echo "Git branch " $gitbranch
echo "Git commit " $gitcommit
cd -
echo "********************************************************"

module load gcc/12.2.0
module load openmpi/4.1.5-cuda12.3
module load cuda/12.3 

nvidia-smi

cat << EOF > gpu_launch.sh
#!/bin/bash

# Compute the raw process ID for binding to GPU and NIC
lrank=$((SLURM_PROCID % SLURM_NTASKS_PER_NODE))

# Bind the process to the correct GPU and NIC
export CUDA_VISIBLE_DEVICES=${lrank}
export UCX_NET_DEVICES=mlx5_${lrank}:1

$@
EOF

chmod +x ./gpu_launch.sh

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

plasma_den=1e18
use_hole_fill=1

for gpus in 4 2 1; do
    for config in 48000 96000; do
        for gpu_red_arrays in 1024 512 256 128 64 32 16 8 4 2 1; do
            for run in 1 2; do
                (( actual_config=config*gpus ))
                echo "Running CUDA BLOCK " $actual_config $run $gpus $gpu_red_arrays

                folder=$runFolder/$config"_mpi"

                mkdir -p $folder
                cp $file $folder
                currentfilename=$folder/$configFile
        
                sed -i "s|STRING hdf_filename = <path_to_hdf5_mesh_files>/box_48000.hdf5|STRING hdf_filename = /home/dp360/dp360/dc-lant1/phd/box_mesh_gen/hdf5/box_${actual_config}.hdf5|" ${currentfilename}

                if [ "$use_hole_fill" -eq 1 ]; then
                    sed -i "s/STRING opp_fill = Shuffle_Periodic/STRING opp_fill = HoleFill_All/" ${currentfilename}
                fi
                sed -i "s/REAL plasma_den           = 1.0e18/REAL plasma_den     = ${plasma_den}/" ${currentfilename}
                sed -i "s/INT gpu_reduction_arrays = 16/INT gpu_reduction_arrays = ${gpu_red_arrays}/" ${currentfilename}

                srun --ntasks=$gpus --gres=gpu:$gpus --hint=nomultithread --distribution=block:block gpu_launch.sh $binary $currentfilename > $folder/log_G${gpus}_M${actual_config}_D${plasma_den}_ARR${gpu_red_arrays}_R${run}.log;
            done
        done
    done
done

rm -rf ./gpu_launch.sh

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
