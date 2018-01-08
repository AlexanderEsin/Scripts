#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-64

chmod +x /home/ade110/Scripts/Geo_again/Consensus_trees/Main/5.RBBH_calculate_slave.tcl

group_name=$GROUP
evalue=1e-150

/home/ade110/Scripts/Geo_again/Consensus_trees/Main/5.RBBH_calculate_slave.tcl $group_name $evalue $PBS_ARRAY_INDEX