#!/bin/sh

dir=Consensus_trees/Usearch

chmod +x /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl

tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 1 &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 2 &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 3 &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 4 &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 5 &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 6 &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 7 &
tclsh8.5 /users/aesin/desktop/Scripts/$dir/5.5.HPC.Rbbh_sort_slave.tcl 8 &
