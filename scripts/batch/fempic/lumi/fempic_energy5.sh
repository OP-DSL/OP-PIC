#!/bin/bash -l
#SBATCH --job-name=fem_5             # Job name
#SBATCH --output=fem_energy5.o%j           # Name of stdout output file
#SBATCH --error=fem_energy5.e%j            # Name of stderr error file
#SBATCH --partition=standard-g    # Partition (queue) name
#SBATCH --nodes=5                    # Total number of nodes 
#SBATCH --ntasks-per-node=8          # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=8            # Allocate one gpu per MPI rank
#SBATCH --time=0-00:30:00            # Run time (d-hh:mm:ss)
##SBATCH --mail-type=all             # Send email at begin and end of job
#SBATCH --account=project_465001068  # Project for billing
##SBATCH --mail-user=username@domain.com

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
export MPICH_GPU_SUPPORT_ENABLED=1

module load LUMI/23.09
module load partition/G
module load cray-hdf5-parallel
module load SCOTCH/6.1.3-cpeCray-23.09

export LD_LIBRARY_PATH=/users/lantraza/phd/lib_install/petsc-3.20.5/lib:$LD_LIBRARY_PATH

runFolder=$PWD"/MPI_N"${SLURM_JOB_NUM_NODES}"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath=/users/lantraza/phd/OP-PIC/app_fempic/bin
binary=$binpath"/hip_mpi_hdf5"
echo "Using Binary -> " $binary

cd $binpath
echo "********************************************************"
gitbranch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
gitcommit=$(git log -n 1 $gitbranch)
echo "Git branch " $gitbranch
echo "Git commit " $gitcommit
echo "********************************************************"
cd -

hdfOriginalFolder=/project/project_465001068/box_mesh_gen/hdf5
num_nodes=$SLURM_JOB_NUM_NODES

configFile="box_fempic.param"
file=$PWD'/'$configFile

for run in 1 2 3; do
    for config in 1536000; do
            
        folder=$runFolder/$config"_mpi"
        (( totalGPUs=8*$SLURM_JOB_NUM_NODES ))

        echo "Running MPI " $config $folder

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        copyCommand=$hdfOriginalFolder'/box_'$config'.hdf5 '$folder'/'
        echo 'Copy box_'$config'.hdf5'
        cp $copyCommand
        copyCommandA=$hdfOriginalFolder'/random_100k.dat '$folder'/'
        echo 'Copy random_100k.dat'
        cp $copyCommandA

        escaped_folder="${folder//\//\\/}"
        sed -i "s/STRING hdf_filename = <path_to_hdf5_mesh_file>/STRING hdf_filename = ${escaped_folder}\/box_${config}.hdf5/" ${currentfilename}
        sed -i "s/STRING rand_file    = <path_to_hdf5_mesh_file>/STRING rand_file    = ${escaped_folder}\/random_100k.dat/" ${currentfilename}

        # ---------------------       
        echo "RUNNING -> 1e18 On "$totalGPUs" GPUs with unsafe atomics"
        srun --cpu-bind=${CPU_BIND} ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_G${totalGPUs}_C${config}_D10_UA_R${run}.log;

        sed -i "s/BOOL use_reg_red = false/BOOL use_reg_red = true/" ${currentfilename}
        echo "RUNNING -> 1e18 On "$totalGPUs" GPUs with seg red"
        srun --cpu-bind=${CPU_BIND} ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_G${totalGPUs}_C${config}_D10_SR_R${run}.log;
        # ---------------------

        echo 'Remove /box_'$config'.hdf5 and random_100k.dat'
        rm $folder'/box_'$config'.hdf5'
        rm $folder'/random_100k.dat'
    done
done

for run in 1 2 3; do
    for config in 1536000; do
            
        folder=$runFolder/$config"_mpi"
        (( totalGPUs=8*$SLURM_JOB_NUM_NODES ))

        echo "Running MPI " $config $folder

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        copyCommand=$hdfOriginalFolder'/box_'$config'.hdf5 '$folder'/'
        echo 'Copy box_'$config'.hdf5'
        cp $copyCommand
        copyCommandA=$hdfOriginalFolder'/random_100k.dat '$folder'/'
        echo 'Copy random_100k.dat'
        cp $copyCommandA

        sed -i "s/REAL plasma_den     = 1e18/REAL plasma_den     = 1.3e18/" ${currentfilename}

        escaped_folder="${folder//\//\\/}"
        sed -i "s/STRING hdf_filename = <path_to_hdf5_mesh_file>/STRING hdf_filename = ${escaped_folder}\/box_${config}.hdf5/" ${currentfilename}
        sed -i "s/STRING rand_file    = <path_to_hdf5_mesh_file>/STRING rand_file    = ${escaped_folder}\/random_100k.dat/" ${currentfilename}

        # ---------------------       
        echo "RUNNING -> 1e18 On "$totalGPUs" GPUs with unsafe atomics"
        srun --cpu-bind=${CPU_BIND} ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_G${totalGPUs}_C${config}_D13_UA_R${run}.log;

        sed -i "s/BOOL use_reg_red = false/BOOL use_reg_red = true/" ${currentfilename}
        echo "RUNNING -> 1e18 On "$totalGPUs" GPUs with seg red"
        srun --cpu-bind=${CPU_BIND} ${binary} ${currentfilename} | tee $folder/log_N${num_nodes}_G${totalGPUs}_C${config}_D13_SR_R${run}.log;
        # ---------------------

        echo 'Remove /box_'$config'.hdf5 and random_100k.dat'
        rm $folder'/box_'$config'.hdf5'
        rm $folder'/random_100k.dat'
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"