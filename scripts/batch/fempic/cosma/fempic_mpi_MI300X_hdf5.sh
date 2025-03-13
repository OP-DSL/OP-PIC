#!/bin/bash
echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

use_seg_red=0
echo "Using Segmented reductions =" $use_seg_red

runFolder=$PWD"/MPI_SR"${use_seg_red}"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=/cosma/home/do018/dc-lant1/phd/OP-PIC/app_fempic_cg/bin
binary=$binpath"/hip_mpi_hdf5"

configFile="coarse.param"
file=/cosma/home/do018/dc-lant1/phd/OP-PIC/scripts/batch/fempic/cosma/$configFile

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

export HIP_VISIBLE_DEVICES=0,1,2,3,4,5,6,7

plasma_den=1e18
use_hole_fill=1

for gpus in 8 4 2 1; do
    for config in 384000 192000 96000 48000; do
        for gpu_red_arrays in 1024 512 256 128 64 32 16 8 4 2 1; do
            for run in 1 2; do
                (( actual_config=config*gpus ))
                echo "Running MPI BLOCK " $actual_config $run $gpus $gpu_red_arrays

                folder=$runFolder/$config"_mpi"

                mkdir -p $folder
                cp $file $folder
                currentfilename=$folder/$configFile
        
                sed -i "s|STRING hdf_filename = <path_to_hdf5_mesh_files>/box_48000.hdf5|STRING hdf_filename = /cosma/home/do018/dc-lant1/phd/Artifacts/mesh_files/box_${actual_config}.hdf5|" ${currentfilename}

                if [ "$use_seg_red" -eq 1 ]; then
                    sed -i "s/BOOL opp_segmented_red = false/BOOL opp_segmented_red = true/" ${currentfilename}
                fi
                if [ "$use_hole_fill" -eq 1 ]; then
                    sed -i "s/STRING opp_fill = Shuffle_Periodic/STRING opp_fill = HoleFill_All/" ${currentfilename}
                fi
                sed -i "s/REAL plasma_den           = 1.0e18/REAL plasma_den     = ${plasma_den}/" ${currentfilename}
                sed -i "s/INT gpu_reduction_arrays = 16/INT gpu_reduction_arrays = ${gpu_red_arrays}/" ${currentfilename}

                mpirun -np $gpus $binary $currentfilename > $folder/log_G${gpus}_M${actual_config}_D${plasma_den}_ARR${gpu_red_arrays}_R${run}.log;
            done
        done
    done
done

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
