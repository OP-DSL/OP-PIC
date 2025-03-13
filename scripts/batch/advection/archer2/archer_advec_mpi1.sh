#!/bin/bash

#SBATCH --job-name=adv_sh1
#SBATCH --time=01:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --account=e723-neptune             
#SBATCH --partition=standard
#SBATCH --qos=short
#SBATCH --exclusive

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module load PrgEnv-gnu
module load cray-hdf5-parallel

runFolder=$PWD"/MPI_"$SLURM_JOB_NUM_NODES"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binary='/work/e723/e723/csrcnj/phd/OP-PIC/app_neso_advection_mdir_cg/bin/mpi'

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

num_nodes=$SLURM_JOB_NUM_NODES

configFile="advec.param"
file=$PWD'/'$configFile

# mesh=512
# for ppc in 425 850 1700; do
#     for run in 1 2; do
        
#         (( actual_ny=mesh*SLURM_JOB_NUM_NODES ))
#         folder=$runFolder/$ppc"_mpi"

#         echo "Running MPI" $file $ppc $folder $nx $actual_ny $ppc $(date +"%Y-%m-%d %H:%M:%S")

#         mkdir -p $folder
#         cp $file $folder
#         currentfilename=$folder/$configFile

#         sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
#         sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
#         sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

#         srun --distribution=block:block --hint=nomultithread ${binary} ${currentfilename} | tee $folder/log_G${num_nodes}_M${mesh}_D${ppc}_ARR1_R${run}.log;
#     done
# done

mesh=256
for ppc in 1700 3400 6800; do
    for run in 1 2; do  
        (( actual_ny=mesh*SLURM_JOB_NUM_NODES ))
        folder=$runFolder/$ppc"_mpi"

        echo "Running MPI" $file $ppc $folder $nx $actual_ny $ppc $(date +"%Y-%m-%d %H:%M:%S")

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
        sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
        sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

        srun --distribution=block:block --hint=nomultithread ${binary} ${currentfilename} | tee $folder/log_G${num_nodes}_M${mesh}_D${ppc}_ARR1_R${run}.log;
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"