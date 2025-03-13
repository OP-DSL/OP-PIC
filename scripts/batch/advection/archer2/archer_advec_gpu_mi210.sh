#!/bin/bash

#SBATCH --job-name=test
#SBATCH --gpus=1
#SBATCH --time=00:30:00
#SBATCH --account=e723-neptune 
#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc
#SBATCH --exclusive

module purge
module load PrgEnv-amd
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan

srun --ntasks=1 rocm-smi

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_SR"${use_seg_red}"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=/work/e723/e723/csrcnj/phd/OP-PIC/app_neso_advection_mdir_cg/bin
binary=$binpath"/hip_mpi"

configFile="advec.param"
file=/work/e723/e723/csrcnj/phd/OP-PIC/scripts/batch/advection/archer2/$configFile

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

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

mesh=512
for gpus in 1; do
    for run in 1 2; do
        for ppc in 425 850 1700; do
            
            (( actual_ny=mesh*gpus ))
            folder=$runFolder/$ppc"_mpi"

            echo "Running MPI" $gpus $ppc $folder $nx $actual_ny $ppc $(date +"%Y-%m-%d %H:%M:%S")

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile

            sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
            sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
            sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

            srun --distribution=block:block --hint=nomultithread --unbuffered ${binary} ${currentfilename} | tee $folder/log_G${gpus}_M${mesh}_D${ppc}_ARR1_R${run}.log;
        done
    done
done

mesh=256
for gpus in 1; do
    for ppc in 1700 3400 6800; do
        for run in 1 2; do  
            (( actual_ny=mesh*gpus ))
            folder=$runFolder/$ppc"_mpi"

            echo "Running MPI" $gpus $ppc $folder $nx $actual_ny $ppc $(date +"%Y-%m-%d %H:%M:%S")

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile

            sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
            sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
            sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

            srun --distribution=block:block --hint=nomultithread --unbuffered ${binary} ${currentfilename} | tee $folder/log_G${gpus}_M${mesh}_D${ppc}_ARR1_R${run}.log;
        done
    done
done

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"

