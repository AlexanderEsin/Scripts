#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -J 1-1670

module load anaconda3/personal

index=$(printf %04d $PBS_ARRAY_INDEX)
evalue=1e-150

tclsh8.6 /home/ade110/Scripts/Staph/RBBH_MCL/3.RBBH_calculate_slave.tcl $evalue $index