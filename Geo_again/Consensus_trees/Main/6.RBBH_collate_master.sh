#!/bin/sh
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=1:mem=50gb

chmod +x /home/ade110/Scripts/Geo_again/Consensus_trees/Main/6.RBBH_collate_slave.tcl

group_name=$GROUP
evalue=1e-10

/home/ade110/Scripts/Geo_again/Consensus_trees/Main/6.RBBH_collate_slave.tcl $group_name $evalue