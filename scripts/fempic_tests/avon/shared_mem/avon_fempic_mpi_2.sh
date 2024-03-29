#!/bin/bash

#SBATCH --job-name=PIC_2
#SBATCH --time=05:00:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=2
#SBATCH --cpus-per-task=1
#SBATCH --partition=compute

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=spread

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module load GCC/10.3.0  OpenMPI/4.1.1
module load PETSc/3.15.1
module load HDF5/1.12.1

runFolder=$PWD"/MPI_"$SLURM_TASKS_PER_NODE"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binPath='/home/dcs/csrcnj/phd/OP-PIC/fempic_mpi/bin/'
binName="mpi_hdf5"
echo "Using Binary -> " $binName

hdfOriginalFolder="/home/dcs/csrcnj/phd/box_mesh_gen/hdf5"
tasks=$SLURM_TASKS_PER_NODE

configFile="/box_avon.param"
file=$PWD'/..'$configFile

for run in 1 2; do
    for config in 12000 24000 48000 96000 192000; do
           
        folder=$runFolder/$config"_mpi"

        echo "Running MPI" $file $config $folder

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        copyCommand=$hdfOriginalFolder'/box_'$config'.hdf5 '$folder'/'
        echo $copyCommand
        cp $copyCommand
        escaped_folder="${folder//\//\\/}"
        sed -i "s/STRING hdf_filename = \/home\/dcs\/csrcnj\/phd\/box_mesh_gen\/hdf5\/box_96000.hdf5/STRING hdf_filename = ${escaped_folder}\/box_${config}.hdf5/" ${currentfilename}

        # ---------------------
        command=${binPath}${binName}' '${currentfilename}' >'$folder'/log_N'${tasks}'_C'${config}'_D10_R'${run}'.log 2>&1'
        echo "RUNNING ->"$command
        srun $command | tee $folder/log_N${tasks}_C${config}_D10_R${run}.log;
        
        sed -i "s/REAL plasma_den     = 1e18/REAL plasma_den     = 0.5e18/" ${currentfilename}
        command=${binPath}${binName}' '${currentfilename}' | tee '$folder'/log_N'${tasks}'_C'${config}'_D5_R'${run}'.log;'
        echo "RUNNING ->"$command
        srun $command | tee $folder/log_N${tasks}_C${config}_D5_R${run}.log;

        sed -i "s/REAL plasma_den     = 0.5e18/REAL plasma_den     = 1e17/" ${currentfilename}
        command=${binPath}${binName}' '${currentfilename}' | tee '$folder'/log_N'${tasks}'_C'${config}'_D1_R'${run}'.log;'
        echo "RUNNING ->"$command
        srun $command | tee $folder/log_N${tasks}_C${config}_D1_R${run}.log;
        # ---------------------

        rmCommand=$folder'/box_'$config'.hdf5'
        echo $rmCommand
        rm $rmCommand
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"