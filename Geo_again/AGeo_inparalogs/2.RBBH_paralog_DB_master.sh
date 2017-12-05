#!/bin/sh
#PBS -l walltime=3:00:00
#PBS -l select=1:ncpus=1:mem=80gb

chmod +x /home/ade110/Scripts/Geo_again/AGeo_inparalogs/2.RBBH_paralog_DB.tcl

evalue=1e-10
echo $evalue

/home/ade110/Scripts/Geo_again/AGeo_inparalogs/2.RBBH_paralog_DB.tcl $evalue