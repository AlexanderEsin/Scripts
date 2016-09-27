#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-128

module load blast+/2.2.28

dir=Geo_v_all

chmod +x /home/ade110/Scripts/$dir/3.HPC.Blastp_slave.tcl

tclsh8.5 /home/ade110/Scripts/$dir/3.HPC.Blastp_slave.tcl $PBS_ARRAY_INDEX