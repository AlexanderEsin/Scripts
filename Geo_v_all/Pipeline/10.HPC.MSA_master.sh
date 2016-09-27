#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=2:mem=40gb
#PBS -J 1-60

module load clustal-omega/1.2 

direct=/home/ade110/Scripts
org=Geo_v_all

chmod +x $direct/$org/10.5.HPC.MSA_slave.tcl

tclsh8.5 $direct/$org/10.5.HPC.MSA_slave.tcl $PBS_ARRAY_INDEX
