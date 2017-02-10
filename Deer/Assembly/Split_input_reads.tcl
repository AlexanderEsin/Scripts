#!/usr/local/bin/tclsh
source ~/Documents/Scripts/General_utils.tcl


set locus_name Beta_local
set direct /users/aesin/desktop/Deer/Assembly/Bos_taurus/$locus_name

set file_name Ovirg2Bos_unique_beta_local_1.fasta

###########################################################################

cd $direct

set output_dir $direct/$locus_name\_input_1
file mkdir $output_dir

set input_trim_name [string range $file_name 0 [string last \. $file_name]-1]

openfile $file_name
split_genes [string trim $data]
set file_split_counter 1

puts "Splitting $file_name into [llength $genes] files with one read per file ..."

foreach gene $genes {
	set out [open $output_dir/$input_trim_name\_seq_$file_split_counter\.fasta w]
	puts $out $gene
	close $out
	incr file_split_counter
}

puts "Splitting done for: $file_name"