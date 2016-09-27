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

###########################################################################

cd $direct

openfile all_prots.txt
set faas $data

openfile hg_can.tsv

set glist $data
set glist [split $glist \n]
regsub -all "{}" $glist {} glist

set pata {>}
set patb {\s+?.+?\n\n+?}
set patc {\(+?.+?\)+?}
set i 1

foreach fam $glist {
	set header [string first \t $fam]
	set groupno [string trim [string range $fam 0 $header]]
	set genes [string trim [string range $fam $header end]]
	set genes [split $genes \t]
	foreach gene $genes {
		regsub " $patc" $gene {} gene
		regexp $pata$gene$patb $faas match
		set match [string trim $match]
		set out [open $direct/Group_fastas/$groupno.txt a]
		puts $out $match
		close $out
	}
	puts "$i/[llength $glist]"
	incr i
}

###########################################################################

##Take the potential candidate (sorted) list and pull out the relevant protein sequences.


