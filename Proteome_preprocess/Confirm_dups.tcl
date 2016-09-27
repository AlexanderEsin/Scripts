#!/usr/local/bin/tclsh
###
set org Viral
set direct /users/aesin/desktop/Cleaning_prots/$org
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

proc split_genes {fasta} {
	global genes
	regsub -all {>} [string trim $fasta] {£>} fasta
	set fasta [split $fasta £]
	regsub -all "{}" $fasta {} genes
}

###########################################################################
set flaglist "Organism\tDup_Number\n"
set donelist {}
set i 1

cd $direct/Orig_faa
set orig_faas [glob *faa]

cd $direct/MSA_check
set orgs [lsort -dictionary [glob -type d *]]

foreach org $orgs {

	cd $direct/Orig_faa
	openfile $org
	split_genes $data
	set new_gene_list $genes

	cd $direct/MSA_check/$org
	set dups [glob -type d *]

	foreach dup $dups {
		cd $direct/MSA_check/$org/$dup
		set input [lindex [glob *.faa] 0]
		set output "[string range $input 0 end-4].fas"

		catch {exec muscle -in $input -phyiout $output}

		openfile $output

		if {[regexp -all -- "\-" $data] > 0} {
			cd $direct/MSA_check/$org
			file rename $dup FLAG_$dup
			append flaglist "$org\t$dup\n"
		} else {
			openfile $input
			split_genes $data
			set genes [lrange $genes 1 end]
			foreach gene $genes {
				set id [string range $gene 0 [string first " " $gene]]
				set matches [lsearch -all -glob $new_gene_list $id*]

				if {[llength $matches] == 1} {
					set index [lindex $matches 0]
					set new_gene_list [lreplace $new_gene_list $index $index]
				}
			}
		}
	}
	set out [open $direct/Clean_faa/$org w]
	puts $out [join $new_gene_list {}]
	close $out
	lappend donelist $org
	puts "$i/[llength $orgs]"
	incr i
}
set out2 [open $direct/Flagged.tsv w]
puts $out2 $flaglist
close $out2

###########################################################################

foreach orig $orig_faas {
	if {[string first $orig $donelist] == -1} {
		file copy $direct/Orig_faa/$orig $direct/Clean_faa
	}
}

###########################################################################
