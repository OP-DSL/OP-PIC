#!/bin/bash

#SBATCH --job-name=cabN1
#SBATCH --time=00:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --account=e723-neptune             
#SBATCH --partition=standard
#SBATCH --qos=short

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_"$SLURM_JOB_NUM_NODES"_ADV_S_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath='/work/e723/e723/csrcnj/phd/OP-PIC/app_cabanapic_cg/bin/'
binary=$binpath'mpi'
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


num_nodes=$SLURM_JOB_NUM_NODES

configFile="cabana.param"
file=$PWD/$configFile

base_nz=60
for ppc in 750 1500 3000; do
    for run in 1 2; do

        echo $file $ppc $(date +"%Y-%m-%d %H:%M:%S")

        folder=$runFolder/$ppc"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        (( nz=$base_nz*$SLURM_JOB_NUM_NODES ))

        sed -i "s/INT nz = 60/INT nz = ${nz}/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 750/INT num_part_per_cell = ${ppc}/" ${currentfilename}
        sed -i "s/INT domain_expansion = 1/INT domain_expansion = ${SLURM_JOB_NUM_NODES}/" ${currentfilename}

        srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores $binary $currentfilename  > $folder/log_G${SLURM_JOB_NUM_NODES}_M${base_nz}_D${ppc}_ARR1_R${run}.log;
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
