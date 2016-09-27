proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

cd ~/desktop/proteome_refseq_dls/bacteria_curated

openfile bacteria_curated_list.txt
set curated_list [split $data \n]
regsub -all "{}" $curated_list {} curated_list

openfile assembly_summary_refseq.txt
set reflist $data

set pata {.+_+?[0-9]*}
set patb {.+?\n+?}

set named_geobac ""

foreach org $curated_list {
	regexp -line $pata $org hit
	set hit [string range $hit 0 14]
	regexp -line $hit$patb $reflist match
	if {[regexp -nocase "Geobacillus" $match] == 1} {
		set entry [split $match \t]
		set name [lindex $entry 7]
		append named_geobac "$org\t$name\n"
	} else {
	}
}

set out [open ~/desktop/proteome_refseq_dls/bacteria_curated/named_geobac.txt a]
puts $out $named_geobac
close $out 