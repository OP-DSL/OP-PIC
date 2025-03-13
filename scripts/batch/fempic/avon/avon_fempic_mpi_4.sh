#!/bin/bash --login
#SBATCH --job-name=PIC_4
#SBATCH --nodes=4
#SBATCH --ntasks-per-node=48
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --partition=compute

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_"$SLURM_JOB_NUM_NODES"_ADV_S_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder" $runFolder

binpath='/home/dcs/csrcnj/phd/neso_test2/OP-PIC/app_fempic_cg/bin/'
binary=$binpath"mpi_hdf5"
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

module purge

# module load GCC/10.3.0  OpenMPI/4.1.1
# module load PETSc/3.15.1
# module load HDF5/1.12.1

module load intel/2022a
module load HDF5/1.13.1

# unset I_MPI_PMI_LIBRARY

export PETSC_INSTALL_PATH=/home/dcs/csrcnj/lib_install/petsc-3.21.0_intel2022a
export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH

# export I_MPI_PMI_LIBRARY=/usr/lib/x86_64-linux-gnu/libpmi.so
# export I_MPI_PMI_LIBRARY=/usr/lib64/pmix/lib/libpmi.so

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close
# export OMP_NUM_THREADS=1

num_nodes=$SLURM_JOB_NUM_NODES

configFile="box_avon.param"
file=$PWD/$configFile

for run in 1 2; do
    for mesh in 48000 96000; do
    
        echo $file $mesh  $(date +"%Y-%m-%d %H:%M:%S")

        folder=$runFolder/$mesh"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        (( mesh_full=$mesh*$SLURM_JOB_NUM_NODES ))

        sed -i "s|STRING hdf_filename = <path_to_hdf5_mesh_file>|STRING hdf_filename = /home/dcs/csrcnj/phd/box_mesh_gen/hdf5_8/box_${mesh_full}.hdf5|" ${currentfilename}

        srun $binary $currentfilename > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D10_ARR1_R${run}.log;

        sed -i "s/REAL plasma_den           = 0.9625e18/REAL plasma_den     = 1.3e18/" ${currentfilename}

        srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename  > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D13_ARR1_R${run}.log;
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
