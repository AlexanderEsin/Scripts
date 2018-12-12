#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=10:mem=50gb
#PBS -J 1-2

module load mowgli/v2.0A

penalty=$PENALTY

tclsh8.6 /home/ade110/Scripts/Bacillus/Mowgli/2.Run_mowgli_slave.tcl $penalty $PBS_ARRAY_INDEX >> trackMow$penalty\.txt