#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=6gb
#PBS -J 1-21

module load mowgli

chmod +x /home/ade110/Scripts/Mowgli/SPR/2.HPC_SPR_Mowgli.tcl

tclsh8.5 /home/ade110/Scripts/Mowgli/SPR/2.HPC_SPR_Mowgli.tcl $PBS_ARRAY_INDEX