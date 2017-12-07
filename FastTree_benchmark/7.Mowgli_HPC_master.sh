#!/bin/sh
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=1:mem=50gb
#PBS -J 1-597

module load mowgli/v2.0A
penalty=6

chmod +x /home/ade110/Scripts/FastTree_benchmark/7.Mowgli_HPC_slave.tcl
/home/ade110/Scripts/FastTree_benchmark/7.Mowgli_HPC_slave.tcl $PBS_ARRAY_INDEX $penalty
