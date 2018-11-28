#!/usr/local/bin/tclsh

## We align each "keyed" fasta file using muscle and then
## stitch the alignments together by taxid to create a 
## supermatrix.

source ~/Documents/Scripts/General_utils.tcl

set evalue			"1e-50"
set trunc_eval		[string range $evalue 2 end]

## Define the directories and create new output folder
set direct			/Users/aesin/Desktop/Staph
set consens_dir		$direct/Prefilter/Staphylococcus

set groups_dir		$consens_dir/Family_groups/$trunc_eval
set key_fasta_dir	$groups_dir/Group_fastas_key

set output_dir		$groups_dir/Aligned_groups
file mkdir			$output_dir

## // Align each orthologous group fasta files // ##

set key_fasta_files	[glob $key_fasta_dir/*fasta]
set sorted_files	[lsort -dictionary $key_fasta_files]

set muscle_counter	1

foreach fasta_file $sorted_files {
	set file_name	[file tail $fasta_file]
	set group		[string range $file_name 0 end-6]

	puts -nonewline stdout "Aligning group family $group ... $muscle_counter // [llength $sorted_files]\r"
	flush stdout

	## We want to write out the alignment in fasta format
	## This will be easier to stitch and we can convert
	## to phylip later
	catch {exec muscle -in $fasta_file -out $output_dir/$group\.faa}
	incr muscle_counter
}

## // Use the taxids to combine all the alignments // ##

## Each alignment file contains the same set of
## taxids, so we can just use the set in the first
## file to combine the rest
puts "\n"
set aligned_files	[lsort -dictionary [glob $output_dir/*faa]]
set index_file		[lindex $aligned_files 0]

set index_prots		[split_genes [string trim [openfile $index_file]]]
set taxids			{}
foreach prot $index_prots {
	set taxid	[string trim [string range $prot 1 [string first \n $prot]]]
	lappend taxids $taxid
}

## Using this set of taxids, we explore each alignment
## and create a supermatrix alignment

set supermatrix_align	{}
set super_counter 		1
foreach taxid $taxids {
	set fasta_header	">$taxid"
	set taxid_seqs		{}

	puts -nonewline stdout "Concatenating taxid: $taxid ... $super_counter // [llength $taxids]\r"
	flush stdout

	foreach aligned_file $aligned_files {
		set aligned_prots	[split_genes [string trim [openfile $aligned_file]]]
		set taxid_aligned	[lsearch -all -inline $aligned_prots "*$fasta_header\n*"]

		## Check there are no duplicates - just in case
		if {[llength $taxid_aligned] > 1} {
			error "There should not be more than one unique taxid per alignment"
		} else {
			set taxid_aligned	[lindex $taxid_aligned 0]
		}

		## Extract aligned sequence
		set aligned_seq			[string trim [string range $taxid_aligned [string first \n $taxid_aligned] end]]

		## Unformat the fasta into a single line
		## string map 5x faster than regsub here
		set unformat_seq		[string map {\n ""} $aligned_seq]
		lappend taxid_seqs		$unformat_seq
	}
	## Combine all into a single line supermatrix
	set taxid_super		[join $taxid_seqs ""]
	## Break again into a fasta-format 80 char width
	set taxid_linebreak	[linebreak $taxid_super]
	## Add the header and append to global list
	set taxid_fasta		[join [concat $fasta_header $taxid_linebreak] \n]
	lappend supermatrix_align	$taxid_fasta

	incr super_counter
}

## Format the supermatrix alignment and output
set out		[open $groups_dir/Supermatrix_fasta.faa w]
puts $out	[join $supermatrix_align \n]
close $out


## // Format into phylip using seqret // ##
puts "\nFormatting as phylip ...."
catch {exec seqret $groups_dir/Supermatrix_fasta.faa $groups_dir/Supermatrix_phylip.phy -osformat phylip 2> /dev/null}
puts "Formatting done! Output is in $groups_dir/Supermatrix_phylip.phy"



