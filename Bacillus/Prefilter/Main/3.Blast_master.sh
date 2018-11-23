#!/bin/sh
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -J 1-2

module load fix_unwritable_tmp
module load blast+/2.2.30

new_number=$(printf %03d $PBS_ARRAY_INDEX)
echo $new_number

/home/ade110/Scripts/Bacillus/Prefilter/Main/3.Blast_slave.tcl $new_number