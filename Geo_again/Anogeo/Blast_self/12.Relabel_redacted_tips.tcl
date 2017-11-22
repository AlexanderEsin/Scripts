source ~/Documents/Scripts/General_utils.tcl
package require sqlite3

set evalue			1e-150
#set evalue			[lindex $argv 0]
set trunc_eval		[string range $evalue 2 end]

## Define the directories and create new output folder
set direct			/users/aesin/Desktop/Geo_again/Anogeo_analysis
set prot_db_file	$direct/../Genomes/Anogeo/Anogeo_prot_db

set tree_dir		$direct/Family_groups/$trunc_eval/Uncertain_redacted
set bipart_file		$tree_dir/FTtree_redacted.txt

set out_tree_dir	$direct/Anogeo_trees
file mkdir			$out_tree_dir

## Open the protein database
sqlite3 db1 $prot_db_file

## Read in the tree, and isolate all the branch labels
## with a regex pattern. Tip labels always have the following
## format: [(or,][label][:]

set tree_data		[string trim [openfile $bipart_file]]
set tip_regex		{[\(\,]{1}[0-9]+[:]{1}}
set taxid_tips		[regexp -inline -all $tip_regex $tree_data]

## Don't work on the original data
set relabel_tree	$tree_data

foreach taxid_tip $taxid_tips {
	## The leading character can differ either "(" or ","
	set leading_char	[string range $taxid_tip 0 0]
	set trim_taxid		[string range $taxid_tip 1 end-1]

	## Get both the binomial and the strain
	set binom_strain	[db1 eval {select binomial, strain from t1 where taxid = $trim_taxid limit 1}]
	set binomial		[string trim [lindex $binom_strain 0]]
	set strain			[string trim [lindex $binom_strain 1]]

	## If the strain name is included in the binomial, strip it
	if {[string match "*$strain*" $binomial] == 1} {
		set binom_clean		[string trim [string map [list $strain ""] $binomial]]
	} else {
		set binom_clean		$binomial
	}
	
	## Add strain name to binomial and replace all
	## whitespace with underscores
	set tip_name		[join [list $binom_clean $strain] " "]
	set tip_name_und	[string map {" " _} $tip_name]

	## Readd the leading and lagging (:) characters
	set tip_name_final	"$leading_char$tip_name_und\:"

	## Substitute out the numeric taxid for the binomial
	set relabel_tree	[string map [list $taxid_tip $tip_name_final] $relabel_tree]
	puts "OLD TIP: $taxid_tip\t\tNEW TIP: $tip_name_final"
}

## Write out the relabelled tree to the working RAxML directory
set out		[open $tree_dir/FT_tree_redacted_relab.txt w]
puts $out	[string trim $relabel_tree]
close $out

db1 close