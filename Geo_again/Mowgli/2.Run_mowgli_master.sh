#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=125gb
#PBS -J 1-6712

## +1:ncpus=0:mem=1mb
module load mowgli/v2.0A
penalty=$PENALTY

chmod +x /home/ade110/Scripts/Geo_again/Mowgli/2.Run_mowgli_slave.tcl

/home/ade110/Scripts/Geo_again/Mowgli/2.Run_mowgli_slave.tcl $penalty $PBS_ARRAY_INDEX 