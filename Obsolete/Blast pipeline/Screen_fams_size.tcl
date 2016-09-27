###
set direct ~/desktop/Geobac
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct

openfile Master_groups_weight.txt
set glist $data

openfile faa_ids.txt
set prot_master $data

set pata {>}
set patb {\s+?.+?\]+?}
set patd {\(+?.+?\)+?}

set g 1
append new_glist "Group no.\tGenes\n"

set genefams [split $glist \n]
regsub -all "{}" $genefams {} genefams
foreach fam $genefams {
	set genes [split $fam \t]
	set desig_fam ""
	set i 1
	if {[llength $genes] > 21} {
		puts "Too long"
	} elseif {[llength $genes] < 21} {
		puts "Too short"
	} else {
		foreach gene $genes {
			regexp $pata$gene$patb $prot_master match
			regexp $patd $match desig
			append desig_fam "$gene\_$desig\t"
			puts "$i/[llength $genes]"
			incr i
		}
		set out [open $direct/Glist.tsv a]
		puts $out "$g\t$desig_fam"
		close $out
		incr g
	}
}


###########################################################################


