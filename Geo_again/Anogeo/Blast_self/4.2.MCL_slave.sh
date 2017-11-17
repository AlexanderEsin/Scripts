#!/bin/bash

ival=$1
evalue=$2

echo "I-value == $ival"
echo "E-value == $evalue"

direct="/Users/aesin/Desktop/Geo_again/Anogeo_analysis/RBBH"

master_RBBH_dir="$direct/Master_RBBH"
mcl_output_dir="$direct/MCL_groups"
mkdir -p $mcl_output_dir

mcl $master_RBBH_dir/Master_RBBH_$evalue\.txt --abc -I $ival -o $mcl_output_dir/MCL_groups_$ival\_$evalue\.txt -te 20 -scheme 7 --show-log=y --abc-neg-log10 -abc-tf "ceil(200)"