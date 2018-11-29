#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=10gb
#PBS -J 1-7304

module load clustalo/1.2.4

/home/ade110/Scripts/Bacillus/Align_tree/9.Group_alignment_slave.tcl $PBS_ARRAY_INDEX