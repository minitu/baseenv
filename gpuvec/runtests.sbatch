#!/bin/bash

#SBATCH --job-name="gpuvec-tests"
#SBATCH --output="gpuvec-tests.%j.%N.out"
#SBATCH --error="gpuvec-tests.%j.%N.err"
#SBATCH --gres=gpu:v100:1
#SBATCH --partition=solo
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --export=ALL
#SBATCH -t 02:00:00

pwd
# You can set any of these to "none" to turn off that set of tests.
cvector_arch=pgi
chost_arch=power9
#chost_arch=none
cacc_arch=tesla
#srun gpuvec-acc -outfile g-acc.txt

set -x
if [ "$cvector_arch" != "none" ] ; then
    date
    srun ./gpuvec -outfile gpuvec.txt
    date
    srun ./gpuvec-double -outfile gpuvec-double.txt
    date
    srun ./gpuvec-len256k -outfile gpuvec-len256k.txt
    date
    srun ./gpuvec-double-len256k -outfile gpuvec-double-len256k.txt
fi
if [ "$chost_arch" != "none" ] ; then
    date
    srun ./gpuvec-acc -outfile gpuvec-acc.txt
    date
    srun ./gpuvec-acc-len256k -outfile gpuvec-acc-len256k.txt
    date
    srun ./gpuvec-acc-double -outfile gpuvec-acc-double.txt
    date
    srun ./gpuvec-acc-double-len256k -outfile gpuvec-acc-double-len256k.txt
fi
if [ "$cacc_arch" != "none" ] ; then
    date
    srun ./gpuvec-acc-gpu -outfile gpuvec-acc-gpu.txt
    date
    srun ./gpuvec-acc-gpu-len256k -outfile gpuvec-acc-gpu-len256k.txt
    date
    srun ./gpuvec-acc-gpu-double -outfile gpuvec-acc-gpu-double.txt
    date
    srun ./gpuvec-acc-gpu-double-len256k -outfile gpuvec-acc-gpu-double-len256k.txt
    date
    srun ./gpuvec-acc-gpu-kernel -outfile gpuvec-acc-gpu-kernel.txt
    date
    srun ./gpuvec-acc-gpu-kernel-len256k -outfile gpuvec-acc-gpu-kernel-len256k.txt
    date
    srun ./gpuvec-acc-gpu-kernel-double -outfile gpuvec-acc-gpu-kernel-double.txt
    date
    srun ./gpuvec-acc-gpu-kernel-double-len256k -outfile gpuvec-acc-gpu-kernel-double-len256k.txt
fi
