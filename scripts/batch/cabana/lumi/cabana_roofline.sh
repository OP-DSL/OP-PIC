#!/bin/bash

#SBATCH --job-name=roof_1             # Job name
#SBATCH --output=roof_1.o%j           # Name of stdout output file
#SBATCH --error=roof_1.e%j            # Name of stderr error file
#SBATCH --partition=standard-g    # Partition (queue) name
#SBATCH --nodes=1                    # Total number of nodes 
#SBATCH --ntasks-per-node=1          # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=1            # Allocate one gpu per MPI rank
#SBATCH --time=0-00:30:00            # Run time (d-hh:mm:ss)
#SBATCH --mail-type=all              # Send email at begin and end of job
#SBATCH --account=project_465001068  # Project for billing
#SBATCH --mail-user=username@domain.com

num_gpus=$SLURM_GPUS_ON_NODE
echo "Number of GPUs allocated: $num_gpus"

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_ROOF_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath=/users/lantraza/phd/OP-PIC/app_cabanapic/bin
binary=$binpath"/hip_mpi"
echo "Using Binary -> " $binary

configFile="cabana.param"
file=/users/lantraza/phd/OP-PIC/scripts/batch/cabana/lumi/$configFile
echo "Config file -> " $file

cd $binpath
echo "********************************************************"
gitbranch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
gitcommit=$(git log -n 1 $gitbranch)
echo "Git branch " $gitbranch
echo "Git commit " $gitcommit
echo "********************************************************"
cd -

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

CPU_BIND="map_cpu:49,57,17,25,1,9,33,41"
export MPICH_GPU_SUPPORT_ENABLED=1

module load CrayEnv
module load buildtools/23.09
module load PrgEnv-cray/8.4.0
module load cce/16.0.1
module load craype-accel-amd-gfx90a
module load craype-x86-trento
module load cray-hdf5-parallel
module use /pfs/lustrep2/projappl/project_462000125/samantao-public/mymodules
module load rocm/5.4.3 omniperf/2.0.1-rocm-5.4.x

num_nodes=$SLURM_JOB_NUM_NODES

monitor_gpu() {
while true; do
    rocm-smi
    sleep 1
done
}

# monitor_gpu &
# monitor_pid=$!

run=1
gpus=1

# ---------------------------------------------------------
config=750

(( totalGPUs=$gpus*$SLURM_JOB_NUM_NODES ))
echo $file $config

echo "Running MPI BLOCK " $config $run $totalGPUs

folder=$runFolder/$config"_mpi"

mkdir -p $folder
cp $file $folder
currentfilename=$folder/$configFile

(( rnz=$totalGPUs*60 ))
sed -i "s/INT nz = 30/INT nz = ${rnz}/" ${currentfilename}
sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
sed -i "s/STRING cluster = pencil/STRING cluster = block/" ${currentfilename}

srun omniperf profile --name all750 --roof-only --kernel opp_dev_interpolate_mesh_fields_kernel,opp_dev_move_deposit_kernel,opp_dev_accumulate_current_to_cells_kernel,opp_dev_half_advance_b_kernel,opp_dev_update_ghosts_B_kernel,opp_dev_update_ghosts_kernel,opp_dev_advance_e_kernel,opp_dev_compute_energy_kernel --kernel-names -- $binary $currentfilename  | tee $folder/log_all750.log;

# ---------------------------------------------------------
config=1500

(( totalGPUs=$gpus*$SLURM_JOB_NUM_NODES ))
echo $file $config

echo "Running MPI BLOCK " $config $run $totalGPUs

folder=$runFolder/$config"_mpi"

mkdir -p $folder
cp $file $folder
currentfilename=$folder/$configFile

(( rnz=$totalGPUs*60 ))
sed -i "s/INT nz = 30/INT nz = ${rnz}/" ${currentfilename}
sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
sed -i "s/STRING cluster = pencil/STRING cluster = block/" ${currentfilename}

srun omniperf profile --name all1500 --roof-only --kernel opp_dev_interpolate_mesh_fields_kernel,opp_dev_move_deposit_kernel,opp_dev_accumulate_current_to_cells_kernel,opp_dev_half_advance_b_kernel,opp_dev_update_ghosts_B_kernel,opp_dev_update_ghosts_kernel,opp_dev_advance_e_kernel,opp_dev_compute_energy_kernel --kernel-names -- $binary $currentfilename  | tee $folder/log_all1500.log;

# ---------------------------------------------------------
config=3000

(( totalGPUs=$gpus*$SLURM_JOB_NUM_NODES ))
echo $file $config

echo "Running MPI BLOCK " $config $run $totalGPUs

folder=$runFolder/$config"_mpi"

mkdir -p $folder
cp $file $folder
currentfilename=$folder/$configFile

(( rnz=$totalGPUs*60 ))
sed -i "s/INT nz = 30/INT nz = ${rnz}/" ${currentfilename}
sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
sed -i "s/STRING cluster = pencil/STRING cluster = block/" ${currentfilename}

srun omniperf profile --name all3000 --roof-only --kernel opp_dev_interpolate_mesh_fields_kernel,opp_dev_move_deposit_kernel,opp_dev_accumulate_current_to_cells_kernel,opp_dev_half_advance_b_kernel,opp_dev_update_ghosts_B_kernel,opp_dev_update_ghosts_kernel,opp_dev_advance_e_kernel,opp_dev_compute_energy_kernel --kernel-names -- $binary $currentfilename  | tee $folder/log_all3000.log;

# ---------------------------------------------------------


echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
