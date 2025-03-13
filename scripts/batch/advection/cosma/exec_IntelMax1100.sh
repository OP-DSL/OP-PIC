
cp advec.param max1100/

sed -i "s/INT npart_per_cell = 1000/INT npart_per_cell = 1700/" $PWD/max1100/advec.param

./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M256_D1700_ARR1_r1.log
./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M256_D1700_ARR1_r2.log

sed -i "s/INT npart_per_cell = 1700/INT npart_per_cell = 3400/" $PWD/max1100/advec.param

./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M256_D3400_ARR1_r1.log
./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M256_D3400_ARR1_r2.log

sed -i "s/INT npart_per_cell = 3400/INT npart_per_cell = 6800/" $PWD/max1100/advec.param

./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M256_D6800_ARR1_r1.log
./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M256_D6800_ARR1_r2.log



sed -i "s/INT nx = 256/INT nx = 512/" $PWD/max1100/advec.param
sed -i "s/INT ny = 256/INT ny = 512/" $PWD/max1100/advec.param

sed -i "s/INT npart_per_cell = 6800/INT npart_per_cell = 425/" $PWD/max1100/advec.param

./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M512_D425_ARR1_r1.log
./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M512_D425_ARR1_r2.log

sed -i "s/INT npart_per_cell = 425/INT npart_per_cell = 850/" $PWD/max1100/advec.param

./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M512_D850_ARR1_r1.log
./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M512_D850_ARR1_r2.log

sed -i "s/INT npart_per_cell = 850/INT npart_per_cell = 1700/" $PWD/max1100/advec.param

./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M512_D1700_ARR1_r1.log
./gpu_wrapper.sh $HOME/phd/OP-PIC/app_neso_advection_mdir_cg/bin/sycl $PWD/max1100/advec.param > max1100/run_G1_M512_D1700_ARR1_r2.log

