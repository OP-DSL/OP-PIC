#!/bin/bash

#SBATCH --job-name=advec_cuda1
#SBATCH --time=01:00:00
#SBATCH --partition=gpu
#SBATCH --account=dp360
#SBATCH --qos=dev

#SBATCH --nodes=2
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4
##SBATCH --exclusive

##SBATCH --qos=reservation
##SBATCH --reservation=Q187755

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_N"${SLURM_JOB_NUM_NODES}"SR"${use_seg_red}"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=$HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin
binary=$binpath"/cuda_mpi"

configFile="advec.param"
file=$HOME/phd/OP-PIC/scripts/batch/advection/tursa/$configFile

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

export I_MPI_PMI_LIBRARY=/lib64/libpmi.so

module load gcc/12.2.0
module load openmpi/4.1.5-cuda12.3
module load cuda/12.3 

unset  OMPI_MCA_pml
unset  OMPI_MCA_osc
unset  OMPI_MCA_btl


nvidia-smi

cat << EOF > gpu_launch.sh
#!/bin/bash

# Compute the raw process ID for binding to GPU and NIC
lrank=$((SLURM_PROCID % SLURM_NTASKS_PER_NODE))

# Bind the process to the correct GPU and NIC
export CUDA_VISIBLE_DEVICES=${lrank}
export UCX_NET_DEVICES=mlx5_${lrank}:1

$@
EOF

chmod +x ./gpu_launch.sh

# mesh=512
# for gpus in 1; do
#     for run in 1 2; do
#         for ppc in 425 850 1700; do
            
#             (( actual_ny=mesh*gpus ))
#             folder=$runFolder/$ppc"_mpi"

#             echo "Running MPI" $gpus $ppc $folder $nx $actual_ny $ppc $(date +"%Y-%m-%d %H:%M:%S")

#             mkdir -p $folder
#             cp $file $folder
#             currentfilename=$folder/$configFile

#             sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
#             sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
#             sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

#             srun --distribution=block:block --hint=nomultithread --unbuffered gpu_launch.sh ${binary} ${currentfilename} | tee $folder/log_G${gpus}_M${mesh}_D${ppc}_ARR1_R${run}.log;
#         done
#     done
# done

# mesh=256
# for gpus in 1; do
#     for ppc in 1700 3400 6800; do
#         for run in 1 2; do  
#             (( actual_ny=mesh*gpus ))
#             folder=$runFolder/$ppc"_mpi"

#             echo "Running MPI" $gpus $ppc $folder $nx $actual_ny $ppc $(date +"%Y-%m-%d %H:%M:%S")

#             mkdir -p $folder
#             cp $file $folder
#             currentfilename=$folder/$configFile

#             sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
#             sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
#             sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

#             srun --distribution=block:block --hint=nomultithread --unbuffered gpu_launch.sh ${binary} ${currentfilename} | tee $folder/log_G${gpus}_M${mesh}_D${ppc}_ARR1_R${run}.log;
#         done
#     done
# done


mesh=1024
for run in 1; do
    for ppc in 105 210; do
        for dt in 0.5 2.5 4.5; do
            folder=$runFolder/$ppc"_mpi"

            echo "Running MPI" $file $ppc $folder $(date +"%Y-%m-%d %H:%M:%S")

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile

            (( actual_ny=$mesh*$SLURM_JOB_NUM_NODES*4 ))

            sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
            sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
            sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

            sed -i "s|REAL dt = 0.5|REAL dt = ${dt}|" ${currentfilename}
            sed -i "s|STRING opp_dh_mesh_folder  = <dh_folder>|STRING opp_dh_mesh_folder  = /home/dp360/dp360/dc-lant1/phd/box_mesh_gen/dh_meshes/advection/dh_${mesh}/|" ${currentfilename}
            
            srun gpu_launch.sh ${binary} ${currentfilename} > $folder/log_G8_M${mesh}_D${ppc}_ARR${dt}_H0_R${run}.log;

            sed -i "s|BOOL opp_global_move = false|BOOL opp_global_move = true|" ${currentfilename}

            srun gpu_launch.sh ${binary} ${currentfilename} > $folder/log_G8_M${mesh}_D${ppc}_ARR${dt}_H1_R${run}.log;
        done
    done
done

# mesh=512
# for run in 1; do
#     for ppc in 425 850; do
#         for dt in 0.5 2.5 4.5 7.0; do
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
#             sed -i "s|STRING opp_dh_mesh_folder  = <dh_folder>|STRING opp_dh_mesh_folder  = /home/dp360/dp360/dc-lant1/phd/box_mesh_gen/dh_meshes/advection/dh_${mesh}/|" ${currentfilename}
            
#             # srun gpu_launch.sh ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H0_R${run}.log;

#             sed -i "s|BOOL opp_global_move = false|BOOL opp_global_move = true|" ${currentfilename}

#             srun gpu_launch.sh ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H1_R${run}.log;
#         done
#     done
# done

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"

