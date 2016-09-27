set direct ~/desktop/Geobac/Geobac_1out_difference

set no_out 1

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

###########################################################################

cd $direct

openfile Filtered_glist_1out.tsv
set database $data

openfile Filtered_glist.tsv
set glist $data

set glist [split $glist \n]
regsub -all "{}" $glist {} glist

set pata {[0-9]+?\t.+?}
set patb {.+?\n+?}
set patc {[0-9]+?\t}

set new_glist ""

foreach group $glist {
	set query [lindex $group 1]
	set x [regexp -line $pata$query$patb $database hit]
	if {$x == 1} {
	} else {
		append new_glist "$group\n"
	}
}

set new_glist [string trim $new_glist]
set out [open $no_out\out_geobac_diff.tsv w]
puts $out $new_glist
close $out

###########################################################################

cd ~/desktop/Geobac/Filtered_geobac_fastas/

file copy all_prots.txt $direct

cd $direct
file mkdir Group_fastas

openfile all_prots.txt
set faas $data

openfile $no_out\out_geobac_diff.tsv

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

file copy Group_fastas Group_fastas_key_ids