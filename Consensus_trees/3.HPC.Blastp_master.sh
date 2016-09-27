#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-32

module load blast+/2.2.28

dir=Consensus_trees

chmod +x /home/ade110/Scripts/$dir/3.5.HPC.Blastp_slave.tcl

tclsh8.5 /home/ade110/Scripts/$dir/3.5.HPC.Blastp_slave.tcl $PBS_ARRAY_INDEX