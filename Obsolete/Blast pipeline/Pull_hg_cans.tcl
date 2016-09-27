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

openfile glist.tsv

regsub -all "\n\n" $data "\n" data
set glist [split $data \n]
regsub -all "{}" $glist {} glist

set pata {\(+?Bact+?\)+?}
set patb {\(+?Arch+?\)+?}
set patc {\(+?Euk+?\)+?}

set hg_can ""

foreach fam $glist {
	if {[regexp -all $patc $fam] <10 && [regexp -all $patc $fam] > 1} {
		if {[regexp -all $patb $fam] > 10 || [regexp -all $pata $fam] > 1} {
			append hg_can "$fam\n"
		}
	} else {
	}
}

set out [open $direct/hg_can.tsv a]
puts $out $hg_can
close $out

###########################################################################

