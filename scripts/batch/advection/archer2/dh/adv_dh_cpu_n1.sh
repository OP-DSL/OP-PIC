#!/bin/bash

#SBATCH --job-name=adv_sh1
#SBATCH --time=02:20:00
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=128
#SBATCH --cpus-per-task=1

#SBATCH --account=e723-neptune             
#SBATCH --partition=standard
#SBATCH --qos=standard
#SBATCH --exclusive

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

module load PrgEnv-gnu
module load cray-hdf5-parallel

runFolder=$PWD"/NEW_"$SLURM_JOB_NUM_NODES"_ADV_DH_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=/work/e723/e723/csrcnj/phd/OP-PIC/app_neso_advection_mdir_cg/bin
binary=$binpath"/mpi"

configFile="advec_dh.param"
file=/work/e723/e723/csrcnj/phd/OP-PIC/scripts/batch/advection/archer2/dh/$configFile

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

dh_mesh_folder=/work/e723/e723/csrcnj/phd/box_mesh_gen/dh_meshes/advection/

command=srun --distribution=block:block --hint=nomultithread --unbuffered

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
            
#             srun --distribution=block:block --hint=nomultithread --unbuffered ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H0_R${run}.log;

#             sed -i "s|BOOL opp_global_move = false|BOOL opp_global_move = true|" ${currentfilename}

#             srun --distribution=block:block --hint=nomultithread --unbuffered ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H1_R${run}.log;
#         done
#     done
# done

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
            
            srun --distribution=block:block --hint=nomultithread --unbuffered ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H0_R${run}.log;

            sed -i "s|BOOL opp_global_move = false|BOOL opp_global_move = true|" ${currentfilename}

            srun --distribution=block:block --hint=nomultithread --unbuffered ${binary} ${currentfilename} > $folder/log_G${SLURM_JOB_NUM_NODES}_M${mesh}_D${ppc}_ARR${dt}_H1_R${run}.log;
        done
    done
done