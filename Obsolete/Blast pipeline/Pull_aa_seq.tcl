###
set direct ~/desktop/Test_groups
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

cd $direct

#openfile Master_groups.txt
#set glist $data

set prot_input ""
set prot_master ""
set patx {>+?.+?\]+?}

cd $direct/faa_all
set protein_seqs [glob *.faa]
foreach species $protein_seqs {
	openfile $species
	set data [string trim $data]
	regsub -all {>} $data {£>} data
	set data [split $data £]
	regsub -all "{}" $data {} data
	foreach gene $data {
		regexp $patx $gene hit1
		append prot_input "$hit1\n"
	}
}

set out [open $direct/faa_ids.txt a]
puts $out $prot_input
close $out

###########################################################################

regsub -all {>} $prot_input "£\n>" prot_input
set prot_master "$prot_input\n£"

set pata {>}
set patb {\s+?.+?}
set patc {£+?}
set patd {\(+?.+?\)+?}

set g 1
set new_glist ""
append new_glist "Group no.\tGenes\n"

set genefams [split $glist \n]
regsub -all "{}" $genefams {} genefams
foreach fam $genefams {
	set genes [split $fam \t]
	set desig_fam ""
	set i 1
	if {[llength $genes] > 50} {
		puts "Too long"
	} elseif {[llength $genes] == 1} {
		puts "Too short"
	} else {
		foreach gene $genes {
			regexp $pata$gene$patb$patc $prot_master match
			#regsub -all {£} $match {} match
			#set match [string trim $match]
			regexp $patd $match desig
			#set group "$match\n"
			#set out [open $direct/Ortho/Groups/group_$g.faa a]
			#puts $out $group
			#close $out
			append desig_fam "$gene\_$desig\t"
			puts "$i/[llength $genes]"
			incr i
		}
		append new_glist "$g\t$desig_fam\n"
		incr g
	}
}

set out [open $direct/Ortho/Glist.tsv a]
puts $out $new_glist
close $out
##

