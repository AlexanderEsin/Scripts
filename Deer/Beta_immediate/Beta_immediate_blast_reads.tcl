#!/usr/local/bin/tclsh
set direct /users/aesin/desktop/Deer/Assembly/Beta_immediate_blast

###########################################################################

proc openfile {fl} {
    global data
    set in [open $fl r]
    set data [read $in]
    close $in
    return
}

proc split_genes {fasta} {
    global genes
    regsub -all {>} [string trim $fasta] {£>} fasta
    set fasta [split $fasta £]
    regsub -all "{}" $fasta {} genes
    return $genes
}

###########################################################################

cd $direct
set immediate_reads_output $direct/Beta_immediate_reads
file mkdir $immediate_reads_output

openfile HBB_locus_test
set hits_l [split [string trim $data] \n]
set total_number_hits [llength $hits_l]

set beta_mapped_reads_file [glob *unique_beta*]
openfile $beta_mapped_reads_file
set reads_l [split_genes [string trim $data]]
set total_number_reads [llength $reads_l]

set i 1
set done_reads_l {}

foreach hit $hits_l {
	set hit_id [lindex [split [string trim $hit] \t] 0]

	if {[lsearch $done_reads_l $hit_id] > -1} {
		puts "READ already done: $hit_id"
		continue
	} else {
		lappend done_reads_l $hit_id
	}

	set hit_sequences [lsearch -all -inline -glob $reads_l >$hit_id*]

	if {[llength $hit_sequences] > 1} {
		puts "Found more than one read for: $hit_id..."
		set test_seq [lindex $hit_sequences 0]
		set exact_matches [lsearch -all $hit_sequences $test_seq]
		if {[llength $exact_matches] == [llength $hit_sequences]} {
			puts "... but they are identical. Picking first one..."
			set hit_sequence [lindex $hit_sequences 0]

			set out [open $immediate_reads_output/Beta_immediate_read_$i\.fasta w]
			puts $out [string trim $hit_sequence]
			close $out

			incr i
		} else {
			puts "... and they are not identical. Using all ..."
			foreach hit_sequence $hit_sequences {
				set out [open $immediate_reads_output/Beta_immediate_read_$i\.fasta w]
				puts $out [string trim $hit_sequence]
				close $out

				incr i
			}

		}
	} elseif {[llength $hit_sequences] == 0} {
		puts "ERROR: did not find a hit for: $hit_id"
		exit 2
	} else {
		set hit_sequence [lindex $hit_sequences 0]

		set out [open $immediate_reads_output/Beta_immediate_read_$i\.fasta w]
		puts $out [string trim $hit_sequence]
		close $out

		incr i
	}
}

puts "\nTotal number of mapped reads to the 'beta globin' locus: $total_number_reads\nTotal number hits to the immediate HBB proximity: $total_number_hits\nTotal number of reads extracted mapping to the immediate HBB proximity: $i"
