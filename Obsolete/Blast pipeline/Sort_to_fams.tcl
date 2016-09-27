###
set direct ~/desktop/Fungi-Archaea-Bact
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}


cd $direct

set z [glob -nocomplain -type d *]
if {[regexp Ortho $z] == 1} {
} else {
	file mkdir Ortho/Groups
}

cd $direct/Rbbh_all/Bin
set rbbh_list [glob *.txt]
set glist ""
set i 1

foreach file $rbbh_list {
	openfile $file
	set rbbhs $data
	set rhead [string first \n $rbbhs]
	set rbbhs [string range $rbbhs $rhead end]
	set rbbhs [split $rbbhs \n]
	regsub -all "{}" $rbbhs {} rbbhs

	foreach line $rbbhs {
		set geneA [lindex $line 0]
		set geneB [lindex $line 1]

		set patm {.*?\n+?}
		set patn {\n*?.*?}

		if {[regexp -line $patn$geneA$patm $glist hit] == 1 && [regexp -line $geneB $hit] == 0} {
			set trimhit [string trim $hit]
			set rep "$trimhit\t$geneB"
			regsub $trimhit $glist $rep glist
		} elseif {[regexp -line $patn$geneB$patm $glist hit] == 1 && [regexp -line $geneA $hit] == 0} {
			set trimhit [string trim $hit]
			set rep "$trimhit\t$geneA"
			regsub $trimhit $glist $rep glist
		} elseif { [regexp -line $patn$geneA$patm $glist hit] == 0 && [regexp -line $patn$geneB$patm $glist hit2] == 0} {
			append glist "\n$geneA\t$geneB\n"
		}
	}
	puts "$i/[llength $rbbh_list]"
	incr i
}

regsub -all {\n\n} $glist "\n" glist

set out [open $direct/Ortho/orig_glist.txt. a]
puts $out $glist
close $out

set prot_input ""
set prot_master ""

cd $direct/faa_all
set protein_seqs [glob *.faa]
foreach species $protein_seqs {
	openfile $species
	set data [string trim $data]
	append prot_input $data
}

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
	set group ""
	set genes [split $fam \t]
	set desig_fam ""
	foreach gene $genes {
		regexp $pata$gene$patb$patc $prot_master match
		regsub -all {£} $match {} match
		set match [string trim $match]
		regexp $patd $match desig
		append group "$match\n"
		append desig_fam "$gene\_$desig\t"
		}
	append new_glist "$g\t$desig_fam\n"
	set out [open $direct/Ortho/group_$g.faa a]
	puts $out $group
	close $out
	incr g
}


set out [open $direct/Ortho/Groups/Glist.tsv a]
puts $out $new_glist
close $out
##
