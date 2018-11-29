#!/bin/sh
#PBS -l walltime=20:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-1670

module load anaconda3/personal
module load blast+/2.2.30

index=$(printf %04d $PBS_ARRAY_INDEX)
echo $index

tclsh8.6 /home/ade110/Scripts/Staph/BLAST/2.Blast_slave.tcl $index