#!/bin/bash

#SBATCH --job-name=fem_2
#SBATCH --time=02:00:00
#SBATCH --partition=gpu-a100-80
#SBATCH --account=dp360
#SBATCH --qos=standard

#SBATCH --nodes=2
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=32
#SBATCH --gres=gpu:4

#SBATCH --qos=reservation
#SBATCH --reservation=Q187755

# partitions: gpu, gpu-a100-80, gpu-a100-40 | qos: standard, dev

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/CUDA_NEW"$SLURM_JOB_NUM_NODES"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

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

lrank=$((SLURM_PROCID % SLURM_NTASKS_PER_NODE))
echo "lrank: " $lrank
echo "SLURM_PROCID: " $SLURM_PROCID
echo "SLURM_NTASKS_PER_NODE: " $SLURM_NTASKS_PER_NODE

ibv_devices

cat << EOF > gpu_launch.sh
#!/bin/bash

# Compute the raw process ID for binding to GPU and NIC
lrank=$((SLURM_PROCID % SLURM_NTASKS_PER_NODE))
gpu_num=$((lrank / 8))  # 8 ranks per GPU

# Bind the process to the correct GPU and NIC
export CUDA_VISIBLE_DEVICES=${gpu_num}
export UCX_NET_DEVICES=mlx5_${gpu_num}:1

$@
EOF

chmod +x ./gpu_launch.sh

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

plasma_den=0.9625e18
use_hole_fill=1

for gpus in 4; do
    for config in 48000 96000; do
        # for gpu_red_arrays in 1024 512 256 128 64 32 16 8 4 2 1; do
        for gpu_red_arrays in 256 512; do
            for run in 1 2; do
                (( actual_config=$config*$gpus*$SLURM_JOB_NUM_NODES ))
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
                
                (( partitions=$gpus*$SLURM_JOB_NUM_NODES ))
                sed -i "s/INT opp_partitions = 1/INT opp_partitions = ${partitions}/" ${currentfilename}
                (( tasks=$gpus*8 ))
                # gpu_launch.sh 
                srun --hint=nomultithread --distribution=block:block $binary $currentfilename > $folder/log_G${gpus}_M${actual_config}_D${plasma_den}_ARR${gpu_red_arrays}_R${run}.log;
            done
        done
    done
done

rm -rf ./gpu_launch.sh

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
