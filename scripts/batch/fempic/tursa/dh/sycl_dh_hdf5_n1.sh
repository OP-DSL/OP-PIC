#!/bin/bash

#SBATCH --job-name=f-SY-DH1
#SBATCH --time=02:00:00
#SBATCH --partition=gpu-a100-80
#SBATCH --account=dp360
#SBATCH --qos=dev

#SBATCH --nodes=1
#SBATCH --cpus-per-task=1
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

# partitions: gpu, gpu-a100-80, gpu-a100-40 | qos: standard, dev

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/SYCL_DH_N${SLURM_JOB_NUM_NODES}_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=/home/dp360/dp360/dc-lant1/phd/OP-PIC/app_fempic_cg/bin
binary=$binpath"/sycl_mpi_hdf5"

configFile="dh_tursa.param"
file=/home/dp360/dp360/dc-lant1/phd/OP-PIC/scripts/batch/fempic/tursa/dh/$configFile

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
module load cuda/12.3 
. /home/dp360/dp360/dc-lant1/.spack/opt/spack/linux-rhel8-zen2/gcc-12.2.0/intel-oneapi-compilers-2024.2.1-p25odcwhcvz5ucvps6ort6hfx2kxwioo/setvars.sh --include-intel-llvm
export I_MPI_PMI_LIBRARY=/lib64/libpmi.so

nvidia-smi
sycl-ls

cat << EOF > gpu_launch.sh
#!/bin/bash
export ONEAPI_DEVICE_SELECTOR="cuda:$((SLURM_PROCID % 4))"
export UCX_NET_DEVICES=mlx5_$((SLURM_PROCID % 4)):1
$@
EOF

chmod +x ./gpu_launch.sh

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

plasma_den=0.9625e18
use_hole_fill=0

tasks=1
gpus=1
mesh=48000
hdfOriginalFolder=/home/dp360/dp360/dc-lant1/phd/box_mesh_gen/hdf5
(( gpus_full=$gpus*$SLURM_JOB_NUM_NODES ))
(( mesh_full=$mesh*$gpus_full ))

folder=$runFolder/$mesh"_mpi"
currentfilename=$folder/$configFile
escaped_folder="${hdfOriginalFolder//\//\\/}"

mkdir -p $folder

# export OMPI_MCA_errmgr_base_verbose=100
# export OMPI_MCA_mpi_show_handle_leaks=1
# export OMPI_MCA_mpi_abort_print_stack=1
# export OMPI_MCA_btl_base_verbose=100  # Verbose logging for shared memory
# export OMPI_MCA_pml_base_verbose=100  # Verbose logging for point-to-point comm
# export OMPI_MCA_osc_base_verbose=100  # Verbose logging for one-sided comm (RMA)
# export OMPI_MCA_shmem_base_verbose=100  # Logs for shared memory segment issues

# unset  OMPI_MCA_pml
# unset  OMPI_MCA_osc
# unset  OMPI_MCA_btl

cp $file $folder
sed -i "s/STRING hdf_filename = <path_to_hdf5_mesh_file>/STRING hdf_filename = ${escaped_folder}\/box_${mesh_full}.hdf5/" ${currentfilename}
sed -i "s/INT opp_partitions = 1/INT opp_partitions = ${gpus_full}/" ${currentfilename}
sed -i "s/INT gpu_reduction_arrays = 512/INT gpu_reduction_arrays = 32/" ${currentfilename}

srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D0.9625e18_ARR0_r1.log
srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D0.9625e18_ARR0_r2.log

sed -i "s/BOOL opp_global_move = false/BOOL opp_global_move = true/" ${currentfilename}

srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D0.9625e18_ARR1_r1.log
srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D0.9625e18_ARR1_r2.log

sed -i "s/REAL plasma_den           = 0.9625e18/REAL plasma_den           = 3.4e18/" ${currentfilename}
sed -i "s/INT charge_multiple = 1/INT charge_multiple = 200/" ${currentfilename}
sed -i "s/REAL opp_allocation_multiple = 11/REAL opp_allocation_multiple = 5/" ${currentfilename}
sed -i "s/BOOL opp_global_move = true/BOOL opp_global_move = false/" ${currentfilename}

srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D3.4e18_ARR0_r1.log
srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D3.4e18_ARR0_r2.log

sed -i "s/BOOL opp_global_move = false/BOOL opp_global_move = true/" ${currentfilename}

srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D3.4e18_ARR1_r1.log
srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D3.4e18_ARR1_r2.log

sed -i "s/REAL plasma_den           = 3.4e18/REAL plasma_den           = 5.75e18/" ${currentfilename}
sed -i "s/INT charge_multiple = 200/INT charge_multiple = 2000/" ${currentfilename}
sed -i "s/REAL opp_allocation_multiple = 5/REAL opp_allocation_multiple = 3.5/" ${currentfilename}
sed -i "s/BOOL opp_global_move = true/BOOL opp_global_move = false/" ${currentfilename}

srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D5.75e18_ARR0_r1.log
srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D5.75e18_ARR0_r2.log

sed -i "s/BOOL opp_global_move = false/BOOL opp_global_move = true/" ${currentfilename}

srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D5.75e18_ARR1_r1.log
srun --ntasks=$tasks gpu_launch.sh $binary $currentfilename > ${folder}/run_G${gpus_full}_M${mesh}_D5.75e18_ARR1_r2.log

rm -rf ./gpu_launch.sh

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
