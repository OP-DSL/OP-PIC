#!/bin/bash

#SBATCH --job-name=test
#SBATCH --nodes=1
#SBATCH --gpus=1
#SBATCH --time=00:15:00
#SBATCH --account=e723-neptune 
#SBATCH --partition=gpu
#SBATCH --qos=gpu-exc       # gpu-exc, gpu-shd
#SBATCH --exclusive

module purge
module load PrgEnv-amd
module load rocm
module load craype-accel-amd-gfx90a
module load craype-x86-milan

export PETSC_INSTALL_PATH=/work/e723/e723/csrcnj/lib_install/petsc-3.20.0_Milan
export LD_LIBRARY_PATH=$PETSC_INSTALL_PATH/lib:$LD_LIBRARY_PATH

srun --ntasks=1 rocm-smi

currentfilename=box_arche_GPU.param

srun --ntasks=1 --cpus-per-task=1 /work/e723/e723/csrcnj/phd/OP-PIC/app_fempic_cg/bin/hip_mpi_hdf5 ${currentfilename}

# sed -i "s/STRING opp_fill = Shuffle_Periodic/STRING opp_fill = HoleFill_All/" ${currentfilename}

# srun --ntasks=1 --cpus-per-task=1 /work/e723/e723/csrcnj/phd/OP-PIC/app_fempic_cg/bin/hip_mpi_hdf5 ${currentfilename}

# sed -i "s/BOOL opp_segmented_red = false/BOOL opp_segmented_red = true/" ${currentfilename}

# srun --ntasks=1 --cpus-per-task=1 /work/e723/e723/csrcnj/phd/OP-PIC/app_fempic_cg/bin/hip_mpi_hdf5 ${currentfilename}

# sed -i "s/STRING opp_fill = Shuffle_Periodic/STRING opp_fill = HoleFill_All/" ${currentfilename}

# srun --ntasks=1 --cpus-per-task=1 /work/e723/e723/csrcnj/phd/OP-PIC/app_fempic_cg/bin/hip_mpi_hdf5 ${currentfilename}

# sed -i "s/BOOL opp_segmented_red = true/BOOL opp_segmented_red = false/" ${currentfilename}


# srun --ntasks=2 --cpus-per-task=1 /work/e723/e723/csrcnj/phd/OP-PIC/app_fempic_cg/bin/hip_mpi_hdf5 /work/e723/e723/csrcnj/phd/OP-PIC/app_fempic_cg/configs/coarse.param

# export PATH=/work/e723/e723/csrcnj/phd/omniperf/2.0.1/bin:$PATH
# export PYTHONPATH=/work/e723/e723/csrcnj/phd/omniperf/python-libs:$PYTHONPATH

# binary=/work/e723/e723/csrcnj/phd/OP-PIC/app_fempic/bin/hip_mpi
# currentfilename=/work/e723/e723/csrcnj/phd/OP-PIC/app_fempic/configs/coarse.param
# folder=$PWD"/ROOF1"

# srun omniperf profile --name move_kernel -k opp_dev_move_kernel --roof-only -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_move_kernel.log;
# srun omniperf profile --name calc_pos_vel_kernel -k opp_dev_calculate_new_pos_vel_kernel --roof-only -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_calculate_new_pos_vel_kernel.log;
# srun omniperf profile --name inject_ions_kernel -k opp_dev_inject_ions_kernel --roof-only -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_inject_ions_kernel.log;
# srun omniperf profile --name deposit_charge_kernel -k opp_dev_deposit_charge_on_nodes_kernel --roof-only -- ${binary} ${currentfilename} | tee $folder/log_opp_dev_deposit_charge_on_nodes_kernel.log;
# srun omniperf profile --name compute_ncd_kernel -k opp_dev_compute_node_charge_density_kernel --roof-only -- ${binary} ${currentfilename} | tee $folder/opp_dev_compute_node_charge_density_kernel.log;
# srun omniperf profile --name compute_ef_kernel -k opp_dev_compute_electric_field_kernel --roof-only -- ${binary} ${currentfilename} | tee $folder/opp_dev_compute_electric_field_kernel.log;
# srun omniperf profile --name compute_F1vec -k computeF1VectorValuesKernel --roof-only -- ${binary} ${currentfilename} | tee $folder/opp_dev_computeF1VectorValuesKernel.log;
# srun omniperf profile --name compute_Jmat -k computeJmatrixValuesKernel --roof-only -- ${binary} ${currentfilename} | tee $folder/opp_dev_computeJmatrixValuesKernel.log;


# srun --ntasks=1 rocm-smi
