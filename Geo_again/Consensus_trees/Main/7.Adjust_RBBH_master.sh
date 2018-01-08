#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=20gb

chmod +x /home/ade110/Scripts/Geo_again/Consensus_trees/Main/7.Adjust_RBBH_slave.tcl

group_name=$GROUP
evalue=1e-150
echo $evalue

/home/ade110/Scripts/Geo_again/Consensus_trees/Main/7.Adjust_RBBH_slave.tcl $group_name $evalue
