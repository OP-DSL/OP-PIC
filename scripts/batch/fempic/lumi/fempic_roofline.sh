#!/bin/bash -l

#SBATCH --job-name=roof_1             # Job name
#SBATCH --output=roof_1.o%j           # Name of stdout output file
#SBATCH --error=roof_1.e%j            # Name of stderr error file
#SBATCH --partition=standard-g    # Partition (queue) name
#SBATCH --nodes=1                    # Total number of nodes 
#SBATCH --ntasks-per-node=1          # 8 MPI ranks per node, 16 total (2x8)
#SBATCH --gpus-per-node=1           # Allocate one gpu per MPI rank
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

# module load LUMI/23.09
# module load partition/G
# module load cray-hdf5-parallel
# module load SCOTCH/6.1.3-cpeCray-23.09

module load CrayEnv
module load buildtools/23.09
module load PrgEnv-cray/8.4.0
module load cce/16.0.1
module load craype-accel-amd-gfx90a
module load craype-x86-trento
module load cray-hdf5-parallel
module use /pfs/lustrep2/projappl/project_462000125/samantao-public/mymodules
module load rocm/5.4.3 omniperf/2.0.1-rocm-5.4.x

export LD_LIBRARY_PATH=/users/lantraza/phd/lib_install/petsc-3.20.5/lib:$LD_LIBRARY_PATH

runFolder=$PWD"/MPI_ROOF_SR0_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
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

hdfOriginalFolder=/users/lantraza/phd/box_mesh_gen/hdf5
num_nodes=$SLURM_JOB_NUM_NODES

configFile="box_fempic.param"
file=$PWD'/'$configFile

monitor_gpu() {
while true; do
    rocm-smi
    sleep 1
done
}

# monitor_gpu &
# monitor_pid=$!

run=1
# gpus=8
# config=384000
gpus=1
config=48000

folder=$runFolder/$config"_mpi"
totalGPUs=$((gpus*SLURM_JOB_NUM_NODES))

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
# sed -i "s/INT num_steps = 250/INT num_steps = 50/" ${currentfilename}

# ---------------------       
echo "RUNNING -> 1e18 On "$totalGPUs" GPUs"

# srun $binary $currentfilename  | tee $folder/log_all1.log;
# echo "1650Mhz DONE 1"
# srun $binary $currentfilename  | tee $folder/log_all2.log;
# echo "1650Mhz DONE 2"
# srun $binary $currentfilename  | tee $folder/log_all3.log;
# echo "1650Mhz DONE 3"



# srun ${binary} ${currentfilename} | tee $folder/log_opp_dev_all1.log;
# srun ${binary} ${currentfilename} | tee $folder/log_opp_dev_all2.log;
# srun ${binary} ${currentfilename} | tee $folder/log_opp_dev_all3.log;

srun omniperf profile --name allAll --kernel opp_dev_move_kernel,opp_dev_calculate_new_pos_vel_kernel,opp_dev_inject_ions_kernel,opp_dev_deposit_charge_on_nodes_kernel,opp_dev_compute_node_charge_density_kernel,opp_dev_compute_electric_field_kernel,computeF1VectorValuesKernel,computeJmatrixValuesKernel --roof-only --kernel-names -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_all1.log;

# srun omniperf profile --name move_kernel -k opp_dev_move_kernel --roof-only --kernel-names -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_move_kernel.log;
# srun omniperf profile --name calc_pos_vel_kernel -k opp_dev_calculate_new_pos_vel_kernel --roof-only --kernel-names -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_calculate_new_pos_vel_kernel.log;
# srun omniperf profile --name inject_ions_kernel -k opp_dev_inject_ions_kernel --roof-only --kernel-names -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_inject_ions_kernel.log;
# srun omniperf profile --name deposit_charge_kernel -k opp_dev_deposit_charge_on_nodes_kernel --kernel-names --roof-only -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_deposit_charge_on_nodes_kernel.log;
# srun omniperf profile --name compute_ncd_kernel -k opp_dev_compute_node_charge_density_kernel --kernel-names --roof-only -- ${binary} ${currentfilename} | tee $folder/opp_dev_compute_node_charge_density_kernel.log;
# srun omniperf profile --name compute_ef_kernel -k opp_dev_compute_electric_field_kernel --kernel-names --roof-only -- ${binary} ${currentfilename} | tee $folder/opp_dev_compute_electric_field_kernel.log;
# srun omniperf profile --name compute_F1vec -k computeF1VectorValuesKernel --roof-only --kernel-names -- ${binary} ${currentfilename} | tee $folder/opp_dev_computeF1VectorValuesKernel.log;
# srun omniperf profile --name compute_Jmat -k computeJmatrixValuesKernel --roof-only --kernel-names -- ${binary} ${currentfilename} | tee $folder/opp_dev_computeJmatrixValuesKernel.log;

# ---------------------

echo 'Remove /box_'$config'.hdf5 and random_100k.dat'
rm $folder'/box_'$config'.hdf5'
rm $folder'/random_100k.dat'

# kill $monitor_pid
# wait $monitor_pid

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
