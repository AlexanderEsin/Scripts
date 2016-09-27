#!/bin/sh
#PBS -l walltime=72:00:00
#PBS -l select=1:ncpus=16:mem=8gb
#PBS -l place=free:shared:group=board
#PBS -J 1-6

module load raxml/2015-02-17 

org=Geo_v_all

chmod +x /home/ade110/Scripts/$org/12.5.G.RAxML_slave.tcl

tclsh8.5 /home/ade110/Scripts/$org/12.5.G.RAxML_slave.tcl $PBS_ARRAY_INDEX