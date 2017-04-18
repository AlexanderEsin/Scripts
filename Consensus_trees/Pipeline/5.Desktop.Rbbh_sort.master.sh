#!/bin/sh

dir=Consensus_trees

chmod +x /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl

tclsh8.5 /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl 1 &
tclsh8.5 /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl 2 &
tclsh8.5 /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl 3 &
tclsh8.5 /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl 4 &
tclsh8.5 /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl 5 &
tclsh8.5 /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl 6 &
tclsh8.5 /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl 7 &
tclsh8.5 /users/aesin/Dropbox/Scripts/$dir/5.5.Desktop.Rbbh_sort_slave.tcl 8 &

wait
