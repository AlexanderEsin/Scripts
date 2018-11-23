#!/usr/local/bin/tclsh

## In this script we label each protein ID in each proteome with
## a unique suffix corresponding to the genome taxid and start position

## In addition, we remove any protein entries with a "missing protein ID"

source ~/Documents/Scripts/General_utils.tcl

set org			"bacteria"

set direct		/Users/aesin/Desktop/HTa_tonio
set prot_dir	$direct/Proteomes
set raw_in		$prot_dir/Proteome_fastas
set clean_out	$prot_dir/Proteome_clean
file mkdir		$clean_out

# Process the refseq assembly table data so we can find the taxids
set genome_list_dir		$direct/Genomes/$org/Genome_lists
set acc_ass_taxid_tbl	[split [string trim [openfile $genome_list_dir/Acc_ass_taxid_table.tsv]] \n]


progress_init		[llength $acc_ass_taxid_tbl]

set cleaning_table	{"Proteome_file\tAcc_ass\tRaw\tClean\tMissing"}
set process_count	1
set all_prot_ids	{}

foreach assembly $acc_ass_taxid_tbl {

	progress_tick	$process_count

	set acc_ass [lindex $assembly 0]
	set taxid	[lindex $assembly 1]

	set file_name	$acc_ass\_proteome.fasta

	# Read in the fasta file
	set proteome	[openfile $raw_in/$file_name]
	set all_prots	[split_genes $proteome]
	# Count all proteins
	set num_all_p	[llength $all_prots]


	# For each one, check presence of valid ProtID 
	# If present add to "clean" proteome, and a taxid
	set clean_prots		{}
	set all_genome_ids	{}
	set num_miss_prot	0
	foreach protein $all_prots {
		## Header first item is the protID (omit ">")
		set header		[split [string trim [string range $protein 1 [string first \n $protein]]] "|"]
		set protID		[string trim [lindex $header 0]]
		set location	[string trim [lindex $header end]]
		set start		[lindex [split $location \:] 0]

		## If the protein has a valid ID - label it with the taxon ID
		if {$protID ne "missing_protein_id_qualifer"} {

			set new_protID		"$protID\.$taxid\_$start"
			if {[lsearch $all_genome_ids $new_protID] > -1} {
				continue
			}
			set new_protein		[string map [list $protID $new_protID] $protein]
		
			lappend clean_prots $new_protein
			lappend all_genome_ids $new_protID
			lappend all_prot_ids $new_protID
		} else {
			incr num_miss_prot
		}
	}

	# Count "clean" proteins
	set num_clean_p		[llength $clean_prots]
	set clean_proteome	[join $clean_prots \n]

	# Write out the clean proteome - all proteins should have proper IDs
	set out [open $clean_out/$file_name w]
	puts $out $clean_proteome
	close $out

	lappend cleaning_table "$file_name\t$acc_ass\t$num_all_p\t$num_clean_p\t$num_miss_prot"
	incr process_count
}

# Write out the summary counts of all proteins vs clean proteins
set out [open $prot_dir/$org\_cleaning_stats.tsv w]
puts $out [join $cleaning_table \n]
close $out

# Write out a list of all the protIDs
set out [open $prot_dir/$org\_all_protIDs.tsv w]
puts $out [join $all_prot_ids \n]
close $out
