#!/bin/sh
#PBS -l walltime=1:00:00
#PBS -l select=1:ncpus=1:mem=4gb
#PBS -J 1-1691

chmod +x /home/ade110/Scripts/Geo_again/RBBH_MCL/3.RBBH_calculate_slave.tcl

new_number=$(printf %04d $PBS_ARRAY_INDEX)
evalue=1e-10

echo $evalue
echo $new_number

/home/ade110/Scripts/Geo_again/RBBH_MCL/3.RBBH_calculate_slave.tcl $evalue $new_number