#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=1:mem=50gb

chmod +x /home/ade110/Scripts/Geo_again/RBBH_MCL/4.RBBH_collate_slave.tcl

evalue=1e-150
echo $evalue

/home/ade110/Scripts/Geo_again/RBBH_MCL/4.RBBH_collate_slave.tcl $evalue
