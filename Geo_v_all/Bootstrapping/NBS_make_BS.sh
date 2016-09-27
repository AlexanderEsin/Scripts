#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=3:mem=16gb
#PBS -J 1-584

module load raxml/2015-02-17 

org=Geo_v_all/Bootstrapping

chmod +x /home/ade110/Scripts/$org/NBS_make_BS.tcl

tclsh8.5 /home/ade110/Scripts/$org/NBS_make_BS.tcl $PBS_ARRAY_INDEX