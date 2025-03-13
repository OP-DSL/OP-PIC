#!/bin/bash --login
#SBATCH --job-name=PIC_8
#SBATCH --nodes=8
#SBATCH --ntasks-per-node=48
#SBATCH --time=3:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=3700
#SBATCH --partition=compute

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

# export OMP_PLACES=cores
# export OMP_PROC_BIND=close
runFolder=$PWD"/NEW_"$SLURM_JOB_NUM_NODES"_ADV_DH_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder -> " $runFolder

binpath='/home/dcs/csrcnj/phd/neso_test2/OP-PIC/app_neso_advection_mdir_cg/bin/'
binary=$binpath"mpi"
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

module purge
# module load GCC/13.3.0
# module load OpenMPI/5.0.3

module load intel/2022a
module load HDF5/1.13.1

CPU_BIND="map_cpu:0,4,8,12,16,20,24,28,32,36,40,44,1,5,9,13,17,21,23,25,29,33,41,45,2,6,10,14,18,22,26,30,34,38,42,46,3,7,11,15,19,27,31,35,37,39,43,47"

configFile="advec.param"
file=$PWD'/'$configFile

# mesh=512
# for run in 1 2; do
#     for ppc in 425 850 1700; do
#             folder=$runFolder/$ppc"_mpi"

#             echo "Running MPI" $file $ppc $folder $(date +"%Y-%m-%d %H:%M:%S")

#             mkdir -p $folder
#             cp $file $folder
#             currentfilename=$folder/$configFile

#             (( actual_ny=$mesh*$SLURM_JOB_NUM_NODES ))

#             sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
#             sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
#             sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

#             srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR1_R${run}.log;
#     done
# done


# mesh=256
# for run in 1 2; do
#     for ppc in 1700 3400 6800; do
#             folder=$runFolder/$ppc"_mpi"

#             echo "Running MPI" $file $ppc $folder $(date +"%Y-%m-%d %H:%M:%S")

#             mkdir -p $folder
#             cp $file $folder
#             currentfilename=$folder/$configFile

#             (( actual_ny=$mesh*$SLURM_JOB_NUM_NODES ))

#             sed -i "s/INT nx = 256/INT nx = ${mesh}/" ${currentfilename}
#             sed -i "s/INT ny = 256/INT ny = ${actual_ny}/" ${currentfilename}
#             sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = ${ppc}/" ${currentfilename}

#             srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR1_R${run}.log;
#     done
# done

mesh=1024
for run in 1 2; do
    for ppc in 105 210; do
        for dt in 0.5 2.5 4.5 7.0; do
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
            sed -i "s|STRING opp_dh_mesh_folder  = <dh_folder>|STRING opp_dh_mesh_folder  = /home/dcs/csrcnj/phd/box_mesh_gen/dh_meshes/advection/dh_${mesh}/|" ${currentfilename}
            
            srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H0_R${run}.log;

            sed -i "s|BOOL opp_global_move = false|BOOL opp_global_move = true|" ${currentfilename}

            srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H1_R${run}.log;
        done
    done
done

mesh=512
for run in 1 2; do
    for ppc in 425 850; do
        for dt in 0.5 2.5 4.5 7.0; do
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
            sed -i "s|STRING opp_dh_mesh_folder  = <dh_folder>|STRING opp_dh_mesh_folder  = /home/dcs/csrcnj/phd/box_mesh_gen/dh_meshes/advection/dh_${mesh}/|" ${currentfilename}
            
            srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H0_R${run}.log;

            sed -i "s|BOOL opp_global_move = false|BOOL opp_global_move = true|" ${currentfilename}

            srun --distribution=block:block --hint=nomultithread --unbuffered --cpu-bind=cores ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H1_R${run}.log;
        done
    done
done

echo "simulation done"

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"