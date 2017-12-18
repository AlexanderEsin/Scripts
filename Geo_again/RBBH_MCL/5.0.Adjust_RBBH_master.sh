#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=4:mem=20gb

chmod +x /home/ade110/Scripts/Geo_again/RBBH_MCL/5.0.Adjust_RBBH_slave.tcl

evalue=1e-10
echo $evalue

/home/ade110/Scripts/Geo_again/RBBH_MCL/5.0.Adjust_RBBH_slave.tcl $evalue
