#!/bin/sh
#PBS -l walltime=6:00:00
#PBS -l select=1:ncpus=1:mem=6gb
#PBS -J 1-5568

module load mowgli

chmod +x /home/ade110/Scripts/Mowgli/Pipeline/3.HPC_run_mowgli.tcl

tclsh8.5 /home/ade110/Scripts/Mowgli/Pipeline/3.HPC_run_mowgli.tcl $PBS_ARRAY_INDEX