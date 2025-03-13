#!/bin/bash
echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/SEQ_SR"${use_seg_red}"_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")

binpath=/cosma/home/do018/dc-lant1/phd/OP-PIC/app_cabanapic_cg/bin
binary=$binpath"/hip"

configFile="cabana.param"
file=/cosma/home/do018/dc-lant1/phd/OP-PIC/scripts/batch/cabana/cosma/$configFile

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

vel_mult=0.0
expand=1

for config in 750 1500 3000; do
    for run in 1 2 3; do
        for gpus in 1; do

            echo "Running SEQ BLOCK " $config $run $gpus $vel_mult

            folder=$runFolder/$config"_seq"

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile
                   
            (( rnz=$gpus*60*$expand ))
            (( dexp=$gpus*$expand ))
            sed -i "s/INT nz = 60/INT nz = ${rnz}/" ${currentfilename}

            sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
            sed -i "s/INT domain_expansion = 1/INT domain_expansion = ${dexp}/" ${currentfilename}
            # sed -i "s/REAL velY_mult_const = 0.0/REAL velY_mult_const = ${vel_mult}/" ${currentfilename}
            # sed -i "s/REAL velZ_mult_const = 0.0/REAL velZ_mult_const = ${vel_mult}/" ${currentfilename}

            $binary $currentfilename > $folder/log_G${gpus}_D${config}_R${run}.log;
        done
    done
done

expand=2
for config in 750 1500 3000; do
    for run in 1 2 3; do
        for gpus in 1; do

            echo "Running SEQ BLOCK " $config $run $gpus $vel_mult

            folder=$runFolder/$config"_seq"

            mkdir -p $folder
            cp $file $folder
            currentfilename=$folder/$configFile
                   
            (( rnz=$gpus*60*$expand ))
            (( dexp=$gpus*$expand ))
            sed -i "s/INT nz = 60/INT nz = ${rnz}/" ${currentfilename}

            sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${config}/" ${currentfilename}
            sed -i "s/INT domain_expansion = 1/INT domain_expansion = ${dexp}/" ${currentfilename}
            # sed -i "s/REAL velY_mult_const = 0.0/REAL velY_mult_const = ${vel_mult}/" ${currentfilename}
            # sed -i "s/REAL velZ_mult_const = 0.0/REAL velZ_mult_const = ${vel_mult}/" ${currentfilename}

            $binary $currentfilename > $folder/log_G${gpus}_D${config}ex_R${run}.log;
        done
    done
done

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
