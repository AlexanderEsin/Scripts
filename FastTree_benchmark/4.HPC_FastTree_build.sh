#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=3:mem=16gb
#PBS -J 1-50

module load fasttree

# Make sure permissions are correct
chmod +x /home/ade110/Scripts/FastTree_benchmark/4.HPC_FastTree_build.tcl

# Run script
/home/ade110/Scripts/FastTree_benchmark/4.HPC_FastTree_build.tcl $PBS_ARRAY_INDEX