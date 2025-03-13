#!/bin/bash

#SBATCH --job-name=test
#SBATCH --gpus=1
#SBATCH --time=01:00:00
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

runFolder=$PWD"/NEW_"$SLURM_JOB_NUM_NODES"_ADV_DH_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=/work/e723/e723/csrcnj/phd/OP-PIC/app_neso_advection_mdir_cg/bin
binary=$binpath"/hip"

configFile="advec_dh.param"
file=/work/e723/e723/csrcnj/phd/OP-PIC/scripts/batch/advection/archer2/dh/$configFile

echo "Creating running folder -> " $runFolder
echo "Using Binary -> " $binary

dh_mesh_folder=/work/e723/e723/csrcnj/phd/box_mesh_gen/dh_meshes/advection/

command=srun --distribution=block:block --hint=nomultithread --unbuffered

# SLURM_JOB_NUM_NODES=1
# mesh=512
# for run in 1 2; do
#     for ppc in 425 850; do # 105 210
#         for dt in 0.5 2.5 4.5 7.0; do # 0.5 2.5 4.5 7.0
#             folder=$runFolder/$ppc"_mpi"

#             echo "Running MPI" $file $ppc $folder $(date +"%Y-%m-%d %H:%M:%S")

#             mkdir -p $folder
#             cp $file $folder
#             currentfilename=$folder/$configFile

#             (( actual_ny=$mesh*$SLURM_JOB_NUM_NODES ))

#             sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
#             sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
#             sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

#             sed -i "s|REAL dt = 0.5|REAL dt = ${dt}|" ${currentfilename}
#             sed -i "s|STRING opp_dh_mesh_folder  = <dh_folder>|STRING opp_dh_mesh_folder  = ${dh_mesh_folder}dh_${mesh}/|" ${currentfilename}
            
#             ${command} ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H0_R${run}.log;

#             sed -i "s|BOOL opp_global_move = false|BOOL opp_global_move = true|" ${currentfilename}

#             ${command} ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H1_R${run}.log;
#         done
#     done
# done

SLURM_JOB_NUM_NODES=1
mesh=1024
for run in 1 2; do
    for ppc in 105 210; do # 
        for dt in 0.5 2.5 4.5 7.0; do # 0.5 2.5 4.5 7.0
            folder=$runFolder/$ppc"_mpi"

            echo "Running MPI" $file $ppc $folder $(date +"%Y-%m-%d %H:%M:%S")

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile

            (( actual_ny=$mesh*$SLURM_JOB_NUM_NODES ))

            sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
            sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
            sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

            sed -i "s|REAL dt = 0.5|REAL dt = ${dt}|" ${currentfilename}
            sed -i "s|STRING opp_dh_mesh_folder  = <dh_folder>|STRING opp_dh_mesh_folder  = ${dh_mesh_folder}dh_${mesh}/|" ${currentfilename}
            
            ${command} ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H0_R${run}.log;

            sed -i "s|BOOL opp_global_move = false|BOOL opp_global_move = true|" ${currentfilename}

            ${command} ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H1_R${run}.log;
        done
    done
done