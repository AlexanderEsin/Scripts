#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -J 1-7303

module load anaconda3/personal

tclsh8.6 /home/ade110/Scripts/Bacillus/Align_tree/10.Group_fastTree_slave.tcl $PBS_ARRAY_INDEX