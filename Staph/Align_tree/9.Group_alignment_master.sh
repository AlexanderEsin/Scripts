#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=5gb
#PBS -J 1-5018

module load anaconda3/personal
module load clustalo/1.2.4

tclsh8.6 /home/ade110/Scripts/Staph/Align_tree/9.Group_alignment_slave.tcl $PBS_ARRAY_INDEX