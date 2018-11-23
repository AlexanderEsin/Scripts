#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-1691

module load anaconda3/personal

index=$(printf %04d $PBS_ARRAY_INDEX)
evalue=1e-10

/home/ade110/Scripts/Bacillus/RBBH_MCL/3.RBBH_calculate_slave.tcl $evalue $index