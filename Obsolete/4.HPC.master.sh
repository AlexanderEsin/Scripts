#!/bin/sh
#PBS -l walltime=70:00:00
#PBS -l select=1:ncpus=1:mem=2gb
#PBS -J 1-64

module load blast+/2.2.28

org=Geo_v_all

chmod +x $SCRATCH/$org/3.HPC.Blastp.tcl

tclsh8.5 $SCRATCH/$org/3.HPC.Blastp.tcl $PBS_ARRAY_INDEX
