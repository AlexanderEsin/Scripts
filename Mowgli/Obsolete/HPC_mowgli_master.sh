#!/bin/sh
#PBS -l walltime=1:59:59
#PBS -l select=1:ncpus=1:mem=25gb
#PBS -J 1-56

module load mowgli#

chmod +x /home/ade110/Scripts/Mowgli/HPC_run_mowgli.tcl

tclsh8.5 /home/ade110/Scripts/Mowgli/HPC_run_mowgli.tcl $PBS_ARRAY_INDEX