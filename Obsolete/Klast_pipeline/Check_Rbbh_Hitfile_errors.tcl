####
set org Archaea
####

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

cd ~/desktop/Klast_$org\_all/Rbbh

openfile RBBH_Hitfile_KLAST.txt

set firstn [string first \n $data]
set header [string range $data 0 $firstn]
set hitlist [string range $data $firstn end]

set hitlist [split $hitlist \n]
regsub -all "{}" $hitlist {} hitlist
set pata "\t+?0\t+?"

set fail ""

foreach hit $hitlist {
	if {[regexp -line $pata $hit] == 1} {
		append fail $hit
	} else {
	}
}