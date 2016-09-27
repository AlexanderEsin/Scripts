set directory /users/aesin/desktop/Deer/Mapping
set species HBB_contigs

##########################

proc openfile {fl} {
    global data
    set in [open $fl r]
    set data [read $in]
    close $in
    return
}

##########################

cd $directory/$species
file mkdir Seed_reads

set sam_1_mis_files [glob -nocomplain HBB*_1*]

foreach read_file $sam_1_mis_files {
	
	set out_fasta_name "[string range $read_file 0 end-4].fasta"
	openfile $read_file
	set reads [split [string trim $data] \n]

	set out_fasta_l {}
	foreach read $reads {
		set read_entries [split [string trim $read] \t]
		set output_line ">[lindex $read_entries 2]:[lindex $read_entries 3]\n[lindex $read_entries 9]"
		lappend out_fasta_l $output_line
	}

	set out [open $directory/$species/Seed_reads/$out_fasta_name w]
	puts $out [join $out_fasta_l \n]
	close $out
}