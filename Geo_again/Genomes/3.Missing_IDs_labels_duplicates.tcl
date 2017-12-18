#!/usr/local/bin/tclsh

## In this script we label each protein ID in each proteome with
## a unique suffix corresponding to the genome taxid and start position

## In addition, we remove any protein entries with a "missing protein ID"

source ~/Documents/Scripts/General_utils.tcl

set direct		/Users/aesin/Desktop/Geo_again/Proteomes
set raw_in		$direct/Proteome_raw_fastas
set clean_out	$direct/Proteome_clean_fastas
set dup_out		$direct/Proteome_dup_fastas
file mkdir		$clean_out; file mkdir $dup_out

## Process the refseq assembly table data so we can find the taxids
set refseq_annotation_list	[split [string trim [openfile $direct/../Genomes/assembly_summary_refseq_131117.txt]] \n]

## List of proteome files, initiate ticker
set raw_proteomes	[glob $raw_in/*fasta]
set num_raw			[llength $raw_proteomes]
progress_init		$num_raw

## Output table holder
set cleaning_table	{"Proteome_file\tAcc_ass\tRaw\tClean"}
set process_count	1

foreach raw_proteome $raw_proteomes {

	progress_tick	$process_count

	## Identify the taxid based on the file name 
	## (i.e. accession and assembly)
	set file_name	[file tail $raw_proteome]
	set acc_ass		[string range $file_name 0 end-15]
	set entry		[lsearch -all -inline $refseq_annotation_list *$acc_ass*]
	# Make sure we only get one entry!
	if {[llength $entry] > 1 || [llength $entry] == 0} {
		error "Expecting to have one unique hit"
	}
	set taxid		[lindex [split $entry \t] 5]

	## Read in the fasta file
	set proteome	[openfile $raw_proteome]
	set all_prots	[split_genes $proteome]
	## Count all proteins
	set num_all_p	[llength $all_prots]


	## For each one, check presence of valid ProtID 
	## If present add to "clean" proteome, and a taxid
	set clean_prots		{}
	set duplicate_prots	{}
	set processed_ID	{}
	set num_miss_prot	0
	foreach protein $all_prots {
		## Header first item is the protID (omit ">")
		set header		[split [string trim [string range $protein 1 [string first \n $protein]]] "|"]
		set protID		[string trim [lindex $header 0]]
		set location	[string trim [lindex $header end]]
		set start		[lindex [split $location \:] 0]

		## If the protein has a valid ID - label it with the taxon ID
		if {$protID ne "missing_protein_id_qualifer"} {

			## Check the protein_ID has not yet been processed
			## If it has we add it to the "duplicate proteome"
			## We will use this to populate gene families after MCL
			if {[lsearch $processed_ID $protID] > 0} {
				set new_protID		"$protID\.$taxid\_$start"
				set new_protein		[string map [list $protID $new_protID] $protein]
				lappend duplicate_prots $new_protein
				continue
			}

			set new_protID		"$protID\.$taxid\_$start"
			set new_protein		[string map [list $protID $new_protID] $protein]
		
			lappend clean_prots $new_protein
			lappend processed_ID $protID
		} else {
			incr num_miss_prot
		}
	}

	## Count "clean" proteins
	set num_clean_p		[llength $clean_prots]
	set num_dup_p		[llength $duplicate_prots]
	set clean_proteome	[join $clean_prots \n]
	set dupl_proteome	[join $duplicate_prots \n]

	## Write out the clean proteome - all proteins should have proper IDs
	set out [open $clean_out/$file_name w]
	puts $out $clean_proteome
	close $out

	## Write out the duplicate proteome - all proteins should have proper IDs
	set out [open $dup_out/$file_name w]
	puts $out $dupl_proteome
	close $out


	lappend cleaning_table "$file_name\t$acc_ass\t$num_all_p\t$num_clean_p\t$num_dup_p\t$num_miss_prot"
	incr process_count
}

# Write out the summary counts of all proteins vs clean proteins
set out [open $direct/Raw_v_clean_stats.tsv w]
puts $out [join $cleaning_table \n]
close $out


