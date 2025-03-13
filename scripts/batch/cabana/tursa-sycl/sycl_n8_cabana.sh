#!/bin/bash

#SBATCH --job-name=cabanaS8
#SBATCH --time=00:30:00
#SBATCH --partition=gpu
#SBATCH --account=dp360
#SBATCH --qos=standard

#SBATCH --nodes=8
#SBATCH --cpus-per-task=8
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4

# partitions: gpu, gpu-a100-80, gpu-a100-40 | qos: standard, dev

echo "Start date and time: $(date +"%Y-%m-%d %H:%M:%S")"

runFolder=$PWD"/MPI_"$SLURM_JOB_NUM_NODES"_ADV_S_"$(date +"D_%Y_%m_%d_T_%I_%M_%S")
echo "Creating running folder" $runFolder

binpath=/home/dp360/dp360/dc-lant1/phd/OP-PIC/app_cabanapic_cg/bin
binary=$binpath"/sycl_mpi"
echo "Using Binary -> " $binary

echo "Creating running folder -> " $runFolder
echo "********************************************************"
cd $binpath
gitbranch=$(git branch | sed -n -e 's/^\* \(.*\)/\1/p')
gitcommit=$(git log -n 1 $gitbranch)
echo "Git branch " $gitbranch
echo "Git commit " $gitcommit
cd -
echo "********************************************************"

module load gcc/12.2.0
module load cuda/12.3 

. /home/dp360/dp360/dc-lant1/.spack/opt/spack/linux-rhel8-zen2/gcc-12.2.0/intel-oneapi-compilers-2024.2.1-p25odcwhcvz5ucvps6ort6hfx2kxwioo/setvars.sh --include-intel-llvm
export I_MPI_PMI_LIBRARY=/lib64/libpmi.so

export OMP_NUM_THREADS=1
export OMP_PLACES=cores

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

export OMP_NUM_THREADS=1
export OMP_PLACES=cores
export OMP_PROC_BIND=close

configFile="cabana.param"
file=$PWD/$configFile

base_nz=60
red_arr=16

for ppc in 750 1500; do
    for run in 1 2; do

        echo $file $ppc $(date +"%Y-%m-%d %H:%M:%S")

        folder=$runFolder/$ppc"_mpi"

        mkdir -p $folder
        cp $file $folder
        currentfilename=$folder/$configFile

        (( gpus=$SLURM_JOB_NUM_NODES*4 ))
        (( nz=$gpus*$base_nz ))

        sed -i "s/INT nz = 60/INT nz = ${nz}/" ${currentfilename}
        sed -i "s/INT num_part_per_cell = 1500/INT num_part_per_cell = ${ppc}/" ${currentfilename}
        sed -i "s/INT domain_expansion = 1/INT domain_expansion = ${gpus}/" ${currentfilename}

        sed -i "s/INT gpu_reduction_arrays = 16/INT gpu_reduction_arrays = 1/" ${currentfilename}

        srun --hint=nomultithread --distribution=block:block gpu_launch.sh $binary $currentfilename > $folder/log_G${gpus}_M${ppc}_D1_ARR1_R${run}.log;

        sed -i "s/INT gpu_reduction_arrays = 1/INT gpu_reduction_arrays = ${red_arr}/" ${currentfilename}

        srun --hint=nomultithread --distribution=block:block gpu_launch.sh $binary $currentfilename > $folder/log_G${gpus}_M${ppc}_D1_ARR${red_arr}_R${run}.log;
    
        sed -i "s/INT gpu_reduction_arrays = ${red_arr}/INT gpu_reduction_arrays = 32/" ${currentfilename}

        srun --hint=nomultithread --distribution=block:block gpu_launch.sh $binary $currentfilename > $folder/log_G${gpus}_M${ppc}_D1_ARR32_R${run}.log;
    done
done

rm -rf ./gpu_launch.sh

echo "End date and time: $(date +"%Y-%m-%d %H:%M:%S")"
