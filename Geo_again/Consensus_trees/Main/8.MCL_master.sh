#!/bin/bash
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=20:mem=240gb
#PBS -l place=pack:excl:group=board

module load mcl/14.137

group_name=$GROUP
ival=2.0
evalue=-150

echo "I-value == $ival"
echo "E-value == $evalue"

direct="/scratch/ade110/Geo_again/Consensus_groups/$group_name/RBBH"

master_RBBH_dir="$direct/Master_RBBH"
mcl_output_dir="$direct/MCL_groups"
mkdir -p $mcl_output_dir

mcl $master_RBBH_dir/Master_RBBH_adj_$evalue\.txt --abc -I $ival -o $mcl_output_dir/MCL_groups_adj_$ival\_$evalue\.txt -te 20 -scheme 7 --show-log=y --abc-neg-log10 -abc-tf "ceil(200)" &> $mcl_output_dir/MCL_adj_$ival\_$evalue\.log
