set direct ~/desktop/Geobac/Outgroups/First_out_smithii

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################
cd $direct

file mkdir temp

cd $direct/Group_fastas

set pata {\[+?.+?\]+?}
set patb {>+?.+?}
set patc {.+?€+?}

set fams [glob *.txt]
set fams [lsort -dictionary $fams]

set new_glist ""

foreach fam $fams {
	regsub {.txt} $fam {} group_name
	append new_glist "\n$group_name"
	openfile $fam
	if {[regexp "MULTISPECIES" $data] >= 1} {
		regsub -all $group_name $new_glist {} new_glist
	} else {
		set species_names ""
		regsub -all ">" $data "£>" genes
		set genes [split $genes £]
		regsub -all "{}" $genes {} genes
		foreach gene $genes {
			regexp $pata $gene match
			regsub -all {\[} $match {} match
			regsub -all {\]} $match {} match
			append species_names "$match\t"
		}
		set namelist [split $species_names \t]
		regsub -all "{}" $namelist {} namelist
		foreach name $namelist {
			if {[regexp -all $name $namelist hit] >= 2} {
				regsub -all $group_name $new_glist {} new_glist
			}
		}
	}
}
set new_glist [string trim $new_glist]
regsub -all {\n\n[\n]*} $new_glist "\n" new_glist
set out [open $direct/Filtered_glist.tsv w]
puts $out $new_glist
close $out

###########################################################################

cd $direct

set pata {\t+?.+?(\n\n)+?}

openfile Glist.tsv
set groups $data
set groups "\n$groups"

set filtered_glist ""

openfile Filtered_glist.tsv
set fams [split $data \n]
set fams [lsort -dictionary $fams]
regsub -all "{}" $fams {} fams
foreach fam $fams {
	regexp "\n$fam$pata" $groups hit
	set hit [string trim $hit]
	append filtered_glist "\n$hit"
}
set filtered_glist [string trim $filtered_glist]
set out [open $direct/Filtered_glist.tsv w]
puts $out $filtered_glist
close $out

###########################################################################

file rename $direct/Filtered_glist.tsv ~/desktop/Geobac/Filtered_glist_1out.tsv
