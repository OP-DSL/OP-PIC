#!/bin/bash --login
#SBATCH --job-name=PIC_DH1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=48
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --partition=compute


echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_"$SLURM_JOB_NUM_NODES"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath='/home/dcs/csrcnj/phd/neso_test2/OP-PIC/app_fempic_cg/bin/'
binary=$binpath'mpi_hdf5'
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

module load intel/2022a
module load HDF5/1.13.1

# unset I_MPI_PMI_LIBRARY

export PETSC_INSTALL_PATH=/home/dcs/csrcnj/lib_install/petsc-3.21.0_intel2022a
export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH


hdfOriginalFolder=/home/dcs/csrcnj/phd/box_mesh_gen/hdf5N
num_nodes=$SLURM_JOB_NUM_NODES

configFile="dh_box_avon.param"
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
sed -i "s/REAL opp_allocation_multiple = 5/REAL opp_allocation_multiple = 3/" ${currentfilename}
sed -i "s/BOOL opp_global_move = true/BOOL opp_global_move = false/" ${currentfilename}

srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D5.75e18_ARR0_r1.log
srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D5.75e18_ARR0_r2.log

sed -i "s/BOOL opp_global_move = false/BOOL opp_global_move = true/" ${currentfilename}

srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D5.75e18_ARR1_r1.log
srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > ${folder}/run_G${SLURM_JOB_NUM_NODES}_M${mesh}_D5.75e18_ARR1_r2.log

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
