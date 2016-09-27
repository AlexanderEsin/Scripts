####
set direct ~/Desktop/Klast_Fungi_All
####

cd $direct/marked_faa

set z [glob -nocomplain -type d *]
if {[regexp Out $z] == 1} {
} else {
	file mkdir Out
}

set dblist [glob *.faa]
set i 1

foreach db $dblist {
	puts "$i/[llength $dblist]"
	set glist [glob *.faa]
	foreach g $glist {
		regsub _protein.faa $g {} x
		regsub _protein.faa $db {} y
		if {$g == $db} {
			puts "skip"
		} else {
			puts "$x\_vs_$y"
			catch {exec ~/desktop/software/klastrunner-4.3a-distrib-macos/Klast.sh -p plastp -i $g -d $db -o Out/$x&$y.tsv -e 1e-10 -outfmt 2 -H 1 -Q 1}
		}
	}
	incr i
}
file copy -force -- $direct/marked_faa/Out $direct
file delete -force -- $direct/marked_faa/Out
puts "DONE"
