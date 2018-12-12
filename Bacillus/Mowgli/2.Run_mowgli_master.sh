#!/bin/sh
#PBS -l walltime=24:00:00
#PBS -l select=1:ncpus=1:mem=25gb
#PBS -J 1-7298

# Load modules
module load anaconda3/personal
module load mowgli/v2.0A

# Run script $PENALTY is supplied as argument to this shell script
tclsh8.6 /home/ade110/Scripts/Bacillus/Mowgli/2.Run_mowgli_slave.tcl $PENALTY $PBS_ARRAY_INDEX >> /home/ade110/Work/Bacillus/Mowgli/trackMow$PENALTY\.log