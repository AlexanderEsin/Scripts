#!/bin/sh
#PBS -l walltime=70:00:00
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -J 1-32

# -J flag followed by a range of subjobs

#org=Geo_v_all

## Write permission ##
chmod +x $SCRATCH/$org/Unzip_dirs.tcl

## Execute script using the subjob-specific variable (e.g. subjob 1 -> $PBS_ARRAY_INDEX = 1 )
tclsh8.5 $SCRATCH/$org/Unzip_dirs.tcl $PBS_ARRAY_INDEX