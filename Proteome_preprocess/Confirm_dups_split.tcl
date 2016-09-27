###
set org Vertebrate_mammalian
set direct ~/desktop/Cleaning_prots/$org
###

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

#Proc to wait for an any-key input
proc anykey {{msg "Checked the Flagged? Press any key: "}} {
    set stty_settings [exec stty -g]
    exec stty raw -echo
    puts -nonewline $msg
    flush stdout
    read stdin 1
    exec stty $stty_settings
    puts ""
}



###########################################################################
set flaglist "Organism\tDup_Number\n"
set donelist {}
set i 1

cd $direct/Orig_faa
set orig_faas [glob *faa]

cd $direct/MSA_check/[lindex $argv 0]
set orgs [glob -type d *]
set counter_file ""

foreach org $orgs {
	set z 0
	set x 1
	cd $direct/Orig_faa
	openfile $org
	set new_faa $data
	cd $direct/MSA_check/[lindex $argv 0]/$org
	set dups [glob -type d *]
	foreach dup $dups {
		cd $direct/MSA_check/[lindex $argv 0]/$org/$dup
		set input [glob *.faa]
		puts "$x / [llength $dups]"
		openfile $input
		regsub -all ">" $data "£>" genes
		set genes [split [string trim $genes] £]
		regsub -all "{}" $genes {} genes
		set genes [lrange $genes 1 end]
		foreach gene $genes {
			set gene [string trim $gene]
			set genelength [string length $gene]
			set index1 [string first $gene $new_faa]
			set index2 [expr $index1 + $genelength]
			set new_faa [string replace $new_faa $index1 $index2]
			incr z
		}
		incr x
	}
	set out [open $direct/Clean_faa/$org w]
	puts $out $new_faa
	close $out
	lappend donelist $org
	puts "$i/[llength $orgs]"
	incr i
	append counter_file "$org\t$z\n"
}

set out [open $direct/Counter_[lindex $argv 0].tsv w]
puts $out $counter_file
close $out

###########################################################################
