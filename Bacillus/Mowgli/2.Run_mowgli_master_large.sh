#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=120gb
#PBS -J 1-7298

module load anaconda3/personal
module load mowgli/v2.0A

penalty=$PENALTY

tclsh8.6 /home/ade110/Scripts/Bacillus/Mowgli/2.Run_mowgli_slave.tcl $penalty $PBS_ARRAY_INDEX >> trackMow$penalty\.txt