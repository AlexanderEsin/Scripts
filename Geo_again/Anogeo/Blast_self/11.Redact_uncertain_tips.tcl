#!/usr/local/bin/tclsh

## The following taxids - see ./Anogeo_trees/ create uncertainty
## in the Anoxy/Geobacillus species tree. We will remove these from
## the supermatrix alignment and reconstruct the trees to see whether
## the uncertainty is removed.

# 1895648
# 1963026
# 1817674

source ~/Documents/Scripts/General_utils.tcl

set evalue			1e-150
# set evalue		[lindex $argv 0]
set trunc_eval		[string range $evalue 2 end]

## Define the directories and create new output folder
set direct			/users/aesin/Desktop/Geo_again/Anogeo_analysis

set groups_dir		$direct/Family_groups/$trunc_eval
set superm_fasta	$groups_dir/Supermatrix_fasta.faa

set output_dir		$groups_dir/Uncertain_redacted
file mkdir			$output_dir

set redact_taxids	[list 1895648 1963026 1817674]


## Open the supermatrix fasta and remove the
## taxids chosen to be redacted.

set aligns			[split_genes [string trim [openfile $superm_fasta]]]
set redacted_align	{}

foreach align $aligns {
	set taxid	[string trim [string range $align 1 [string first \n $align]]]
	if {[lsearch $redact_taxids $taxid] == -1} {
		lappend redacted_align $align
	}
}

## Write out the redacted alignment
set out		[open $output_dir/Redacted_supermatrix_fasta.faa w]
puts $out	[join $redacted_align \n]
close $out

## // Format into phylip using seqret // ##
puts "\nFormatting as phylip ...."
catch {exec seqret $output_dir/Redacted_supermatrix_fasta.faa $output_dir/Redacted_supermatrix_phylip.phy -osformat phylip 2> /dev/null}
puts "Formatting done! Output is in $output_dir/Supermatrix_phylip.phy"


## // Run FastTree on the redacted alignment // ##
set out		[open $output_dir/FastTree_log.txt a]
catch {exec fasttree -out $output_dir/FTtree_redacted.txt $output_dir/Redacted_supermatrix_phylip.phy >&@$out}
close $out