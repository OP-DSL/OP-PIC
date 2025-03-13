#!/bin/bash

#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=00:15:00
#SBATCH --account=e723-neptune 
#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc       # gpu-exc, gpu-shd
#SBATCH --exclusive

module purge
module load PrgEnv-amd
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan
export PETSC_INSTALL_PATH=/work/e723/e723/csrcnj/lib_install/petsc-3.20.0_Milan
export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/HIP_"$SLURM_JOB_NUM_NODES"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath='/work/e723/e723/csrcnj/phd/OP-PIC/app_fempic_cg/bin/'
binary=$binpath'hip_mpi_hdf5'
echo "Using Binary -> " $binary

echo "Creating running folder -> " $runFolder
echo "Using Binary -> " $binpath$binary
echo "********************************************************"
cd $binpath
gitbranch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
gitcommit=$(git log -n 1 $gitbranch)
echo "Git branch " $gitbranch
echo "Git commit " $gitcommit
cd -
echo "********************************************************"


export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

hdfOriginalFolder=/work/e723/e723/csrcnj/phd/box_mesh_gen/hdf5
num_nodes=$SLURM_JOB_NUM_NODES

configFile="dh_archer.param"
file=$PWD'/'$configFile


mesh=48000
(( mesh_full=$mesh*$SLURM_JOB_NUM_NODES ))

folder=$runFolder/$mesh"_mpi"
currentfilename=$folder/$configFile
escaped_folder="${hdfOriginalFolder//\//\\/}"

mkdir -p $folder


cp $file $folder
sed -i "s/STRING hdf_filename = <path_to_hdf5_mesh_file>/STRING hdf_filename = ${escaped_folder}\/box_${mesh_full}.hdf5/" ${currentfilename}

srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D0.9625e18_ARR0_r1.log
srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D0.9625e18_ARR0_r2.log

sed -i "s/BOOL opp_global_move = false/BOOL opp_global_move = true/" ${currentfilename}

srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D0.9625e18_ARR1_r1.log
srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D0.9625e18_ARR1_r2.log

sed -i "s/REAL plasma_den           = 0.9625e18/REAL plasma_den           = 3.4e18/" ${currentfilename}
sed -i "s/INT charge_multiple = 1/INT charge_multiple = 200/" ${currentfilename}
sed -i "s/REAL opp_allocation_multiple = 11/REAL opp_allocation_multiple = 5/" ${currentfilename}
sed -i "s/BOOL opp_global_move = true/BOOL opp_global_move = false/" ${currentfilename}

srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D3.4e18_ARR0_r1.log
srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D3.4e18_ARR0_r2.log

sed -i "s/BOOL opp_global_move = false/BOOL opp_global_move = true/" ${currentfilename}

srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D3.4e18_ARR1_r1.log
srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D3.4e18_ARR1_r2.log

sed -i "s/REAL plasma_den           = 3.4e18/REAL plasma_den           = 5.75e18/" ${currentfilename}
sed -i "s/INT charge_multiple = 200/INT charge_multiple = 2000/" ${currentfilename}
sed -i "s/REAL opp_allocation_multiple = 5/REAL opp_allocation_multiple = 3.5/" ${currentfilename}
sed -i "s/BOOL opp_global_move = true/BOOL opp_global_move = false/" ${currentfilename}

srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D5.75e18_ARR0_r1.log
srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D5.75e18_ARR0_r2.log

sed -i "s/BOOL opp_global_move = false/BOOL opp_global_move = true/" ${currentfilename}

srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D5.75e18_ARR1_r1.log
srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D5.75e18_ARR1_r2.log

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"

