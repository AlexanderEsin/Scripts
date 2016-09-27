proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

#This proc remakes the 80-wide fasta entries from a flat CDS sequence

proc linebreak {s {width 60}} {
   global res
   set res ""
   while {[string length $s]>$width} {
       set res "$res[string range $s 0 [expr $width - 1]]\n"
       set s [string range $s $width end]
   }
   set res "$res$s"
   return $res
}

cd ~/desktop/Geobac_rnaseq

openfile [lindex [glob *gbk] 0]
set x [regexp -all {£} $data]
set prot_number [regexp -all {/translation} $data]

set total_CDS 0
set glist {}

if {$x == 0} {
	regsub -all {LOCUS} [string trim $data] {£LOCUS} new_data
	set split_by_locus [split $new_data £]

	set number_loci [llength $split_by_locus]

	foreach locus $split_by_locus {
		set locus_name [lindex [split [string range $locus [string first "NODE" $locus] [string first "bp" $locus]] " "] 0]
		
		if {[regexp -all { CDS } $locus] > 0} {
			set locus [string range $locus 0 [string first "BASE COUNT" $locus]-1]
			regsub -all { CDS } $locus {£CDS } locus
			set CDS_set [lrange [split $locus £] 1 end]
			set total_CDS [expr $total_CDS + [llength $CDS_set]]

			set CDS_counter 1

			foreach CDS $CDS_set {
				set fasta_name "$locus_name\_CDS$CDS_counter"

				set map_range [string trim [string range $CDS [string first " " $CDS] [string first "/" $CDS]-1]]
				if {[regexp {complement} $map_range] != 1} {
					set start_end_l [split $map_range " .. "]
					set map_range_name "[lindex $start_end_l 0]:[lindex $start_end_l end] Forward"
				} else {
					set map_range [string range $map_range [string first \( $map_range]+1 [string first \) $map_range]-1]
					set start_end_l [split $map_range ".."]
					set map_range_name "[lindex $start_end_l 0]:[lindex $start_end_l end] Reverse"
				}
				
				set features [split $CDS "/"]

				set product_feature [lsearch -inline -glob $features product*]
				regsub -all {\n} $product_feature {} product_feature
				regsub -all {\s+} $product_feature " " product_feature
				set product_name [string range $product_feature [string first \" $product_feature]+1 [string last \" $product_feature]-1]

				set translation_feature [lsearch -inline -glob $features translation*]
				regsub -all {\s+} $translation_feature {} translation_feature
				set aa_seq [string range $translation_feature [string first \" $translation_feature]+1 [string last \" $translation_feature]-1]
				set aa_seq [linebreak $aa_seq]

				set fasta_header ">$fasta_name | $map_range_name | $product_name\n"
				set fasta_entry "$fasta_header$aa_seq"

				lappend glist $fasta_entry

				incr CDS_counter
			}
		}		
	}
}

set out [open DL33genbank_protein_CDS.faa w]
puts $out [join $glist \n]
close $out

puts "TOTAL CDS: $total_CDS"
puts $number_loci
puts $prot_number
puts [llength $glist]