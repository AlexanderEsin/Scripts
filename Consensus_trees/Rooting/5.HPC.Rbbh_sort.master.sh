#!/bin/sh

dir=Consensus_trees/Rooting
org=$1

chmod +x /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl

tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 1 $org &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 2 $org &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 3 $org &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 4 $org &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 5 $org &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 6 $org &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 7 $org &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 8 $org &
