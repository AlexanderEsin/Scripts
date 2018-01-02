#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=20gb

group_name=$GROUP

chmod +x /home/ade110/Scripts/Geo_again/Consensus_trees/4.Blast_out_rearrange_slave.tcl

/home/ade110/Scripts/Geo_again/Consensus_trees/4.Blast_out_rearrange_slave.tcl $group_name