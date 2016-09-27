proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

cd ~/Desktop/Marked_proteomes/Plants
set olist [glob *.faa]
set pata {€+?.+?€+?}

foreach o $olist {
	openfile $o
	regexp $pata $data line
	regsub -all "€" $line {} line
	set line "($line)"
	regsub -all $pata $data $line data2
	set out [open $o w]
	puts $out $data2
	close $out
}