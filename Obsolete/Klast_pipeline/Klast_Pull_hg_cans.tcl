proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

cd ~/desktop/test_Klast/ortho/groups
openfile glist.tsv
set header [string first \n $data]
set glist [string range $data $header end]
set glist [split $glist \n]
regsub -all "{}" $glist {} glist

set pata {€+?Bact?€+?}
set patb {€+?Arch?€+?}
set patc {€+?Euk?€+?}
set hg_can ""

foreach fam $glist {
	if {[regexp -all $patc $fam] == 1} {
		if {[regexp -all $patb $fam] == 2 || [regexp -all $pata $fam] == 2} {
			append hg_can "$fam\n"
		} else {
		}
	} else {
	}
}

set out [open ~/desktop/test_Klast/ortho/groups/hg_can.tsv a]
puts $out $hg_can
close $out