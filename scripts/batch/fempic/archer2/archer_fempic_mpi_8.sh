#!/bin/bash

#SBATCH --job-name=femN8
#SBATCH --time=00:20:00
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --account=e723-neptune             
#SBATCH --partition=standard
#SBATCH --qos=short

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_"$SLURM_JOB_NUM_NODES"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath='/work/e723/e723/csrcnj/phd/OP-PIC/app_fempic_cg/bin/'
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

module load PrgEnv-gnu
module load cray-hdf5-parallel
export PETSC_INSTALL_PATH=/work/e723/e723/csrcnj/lib_install/petsc-3.20.0
export OPPIC_PATH=/work/e723/e723/csrcnj/phd/OP-PIC/lib_oppic
export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH

hdfOriginalFolder=/work/e723/e723/csrcnj/phd/box_mesh_gen/hdf5
num_nodes=$SLURM_JOB_NUM_NODES

configFile="box_archer.param"
file=$PWD'/'$configFile

# monitor_memory() {
#     while true; do
#         free -h >> memory_usage.log
#         sleep 1
#     done
# }
# monitor_memory &

for run in 1 2; do
    for mesh in 48000 96000; do
           
        folder=$runFolder/$mesh"_mpi"

        echo "Running MPI" $mesh $folder $(date +"%Y-%m-%d %H:%M:%S")

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        (( mesh_full=$mesh*$SLURM_JOB_NUM_NODES ))

        copyCommand=$hdfOriginalFolder'/box_'$mesh_full'.hdf5 '$folder'/'
        echo "COPY ->" $copyCommand
        cp $copyCommand

        escaped_folder="${folder//\//\\/}"
        sed -i "s/STRING hdf_filename = <path_to_hdf5_mesh_file>/STRING hdf_filename = ${escaped_folder}\/box_${mesh_full}.hdf5/" ${currentfilename}

        srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D10_ARR1_R${run}.log;

        sed -i "s/REAL plasma_den           = 0.9625e18/REAL plasma_den     = 1.3e18/" ${currentfilename}

        srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename  > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D13_ARR1_R${run}.log;

        rmCommand=$folder'/box_'$mesh_full'.hdf5'
        rm $rmCommand
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
