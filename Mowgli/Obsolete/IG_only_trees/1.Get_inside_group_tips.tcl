source ~/Dropbox/Scripts/General_utils.tcl
set direct /users/aesin/desktop
package require sqlite3

## Get list of inside taxa ##
cd $direct/Clean_Proteomes/Inside_group
set in_group_fastas [glob *faa]

## Open sql databse ##
cd $direct/Clean_Proteomes
sqlite3 db1 all_prot_geo_db

set in_binomials {}
set counter 0
## Get binomial for each file ##
foreach file $in_group_fastas {
	db1 eval {SELECT binomial FROM t1 WHERE file_name = $file} {
		set binomial "$binomial"
	}
	lappend in_binomials $binomial
	incr counter
}

puts "Total names found: $counter"
db1 close

## Confirm that all these tips can be found in the species tree used for the Mowgli reconciliation ##
# Go to the mowgli dir and open up the species tree #

cd $direct/Mowgli/Species_tree
set tree_data [openfile Ultrametric_species_tree_tipfix.txt]

foreach binomial $in_binomials {

	regsub -all "__" $binomial {_} binomial

	if {[regexp $binomial $tree_data] == 0} {
		puts $binomial
	}
}

set out [open $direct/Mowgli/Inside_group_tips.txt]
puts $out [join $in_binomials \n]
close $out