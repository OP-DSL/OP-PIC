#!/bin/bash

#SBATCH --job-name=cabN64
#SBATCH --time=01:00:00
#SBATCH --nodes=64
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --account=e723-neptune             
#SBATCH --partition=standard
#SBATCH --qos=standard

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

export SRUN_CPUS_PER_TASK=$SLURM_CPUS_PER_TASK

module load PrgEnv-gnu

runFolder=$PWD"/MPI_"$SLURM_JOB_NUM_NODES"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath='/work/e723/e723/csrcnj/phd/OP-PIC/app_cabanapic_cg/bin/'
binary=$binpath'mpi'
echo "Using Binary -> " $binary

cd $binpath
echo "********************************************************"
gitbranch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
gitcommit=$(git log -n 1 $gitbranch)
echo "Git branch " $gitbranch
echo "Git commit " $gitcommit
echo "********************************************************"
cd -

num_nodes=$SLURM_JOB_NUM_NODES

configFile="cabana.param"
file=$PWD/$configFile

for config in 750 1500 3000; do
    for run in 1 2; do

        echo $file $config

        # ****************************************
        echo "Running MPI"

        folder=$runFolder/$config"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        expand=1
        (( rnz=$num_nodes*60*$expand ))
        (( dexp=$num_nodes*$expand ))
        sed -i "s/INT nz = 30/INT nz = ${rnz}/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
        sed -i "s/INT domain_expansion = 1/INT domain_expansion = ${dexp}/" ${currentfilename}
        
        srun --distribution=block:block --hint=nomultithread --unbuffered $binary $currentfilename  | tee $folder/log_N${num_nodes}_Block_D${config}_R${run}.log;
    done
done

# for config in 750 1500 3000; do
#     for run in 1 2; do

#         echo $file $config

#         # ****************************************
#         echo "Running MPI"

#         folder=$runFolder/$config"_mpi"

#         mkdir -p $folder
#         cp $file $folder
#         currentfilename=$folder/$configFile

#         (( rnz=$num_nodes*60 ))
#         sed -i "s/INT nz = 30/INT nz = ${rnz}/" ${currentfilename}
#         sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
#         srun $binary $currentfilename  | tee $folder/log_N${num_nodes}_Pencil_D${config}_R${run}.log;
#     done
# done

# for config in 750 1500 3000; do
#     for run in 1 2; do

#         echo $file $config

#         # ****************************************
#         echo "Running MPI"

#         folder=$runFolder/$config"_mpi"

#         mkdir -p $folder
#         cp $file $folder
#         currentfilename=$folder/$configFile

#         (( rnz=$num_nodes*60 ))
#         sed -i "s/INT nz = 30/INT nz = ${rnz}/" ${currentfilename}
#         sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
#         sed -i "s/STRING cluster = pencil/STRING cluster = cart/" ${currentfilename}
#         srun $binary $currentfilename  | tee $folder/log_N${num_nodes}_Cart_D${config}_R${run}.log;
#     done
# done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"