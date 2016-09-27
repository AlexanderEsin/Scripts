#!/bin/sh
#PBS -l walltime=25:00:00
#PBS -l select=1:ncpus=1:mem=1gb
#PBS -J 1-32

dir=Consensus_trees

chmod +x /home/ade110/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl

tclsh8.5 /home/ade110/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl $PBS_ARRAY_INDEX
