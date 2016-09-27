### This script takes the RBBH_Hitfile and processes some of the data to identify a basic similarity score between Geobacilli and other organisms - dividing them into Bacteria, Archaea, and Eukaryotes. This will only show large scale similarity, e.g. identifying species that are likely close related to the Geobacilli but are not classified as such ###

set direct /users/aesin/Desktop/Consensus_trees
set org Geobac/All_geobac
set rbbh_dir Rbbh_all_-40
set hitfile RBBH_Hitfile_-40.txt

proc openfile {fl} {
	global data
	set in [open $fl r]
	set data [read $in]
	close $in
	return
}

proc commas {var {num 3} {char ,}} {
    set len   [string length $var]
    set first [expr $len - $num]
    set x     {}
    while {$len > 0} {
        # grab left num chars
        set lef [string range $var $first end] 
        if {[string length $x] > 0} {
            set x   "${lef}$char${x}"
        } else {
            set x   $lef
        }
        # grab everything except left num chars
        set var [string range  $var 0 [expr $first -1]]
        set len   [string length $var]
        set first [expr {$len - $num}]
    }
    return $x
}


package require sqlite3

###########################################################################

cd /users/aesin/dekstop/Clean_proteomes
sqlite3 db1 all_prot_geo_db

file mkdir $direct/$org/$rbbh_dir/Rbbh_stats

cd $direct/$org/$rbbh_dir

openfile $hitfile

# Get the hitfile ready to be processed line by line #
regsub -all "\n\n" $data "\n" $data
set data [string range $data [string first \n $data] end]
set hits [split [string trim $data] \n]
regsub -all "{}" $hits {} hits

set bact_hits ""
set arch_hits ""
set euk_hits ""

set line_counter 1
set save_counter 0

foreach hit $hits {
	# Ignore those hits that are between two Geobacilli #
	if {[regexp -all {Geobacillus} $hit] == 1} {
		# Identify the Geobacillus species name #
		set geo_index [string first {Geobacillus_} $hit]
		set geo_name [string trim [string range [string range $hit $geo_index end] 0 [string first \t [string range $hit $geo_index end]]]]

		set line_vals [split $hit \t]

		#Get the subject organism name and its size #
		if {$geo_index == 0} {
			set non_geo_name [lindex $line_vals 1]
			set other_genome_size [lindex $line_vals 8]
		} else {
			set non_geo_name [lindex $line_vals 0]
			set other_genome_size [lindex $line_vals 7]
		}
		# Filter out v. small genomes (as percentage similarity is greatly skewed for these) #
		if {$other_genome_size > 500} {
			set file $non_geo_name\.faa
			# Pull out the binomial name of the organism as appears in the fasta comment line #
			db1 eval {SELECT binomial FROM t1 WHERE file_name = $file} {
				set binomial "$binomial"
			}
			regsub -all {\[} $binomial {} binomial
			regsub -all {\]} $binomial {} binomial
			# Also get the class (Bact, Arch, Euk)
			db1 eval {SELECT class FROM t1 WHERE file_name = $file} {
				set class "$class"
			}
			# The number of RBBH hits for this comparison #
			set hit_number [lindex $line_vals 2]

			# This a vague measure of similarity = number of RBBH hits / size of subject genome #
			set R [expr double($hit_number) / $other_genome_size]

			# Filter results into the appropriate group #
			if {$class == {(Bact)}} {
				append bact_hits "$binomial\t[lindex $line_vals 2]\t$other_genome_size\t$R\n"
			} elseif {$class == {(Arch)}} {
				append arch_hits "$binomial\t[lindex $line_vals 2]\t$other_genome_size\t$R\n"
			} elseif {$class == {(Euk)}} {
				append euk_hits "$binomial\t[lindex $line_vals 2]\t$other_genome_size\t$R\n"
			}
		}

	}
	puts "$line_counter / [llength $hits]"
	incr line_counter
	incr save_counter
}

db1 close

# Output #
cd $direct/$org/$rbbh_dir/Rbbh_stats

set out [open Hits_to_bact.tsv w]
puts $out [string trim $bact_hits]
close $out
set out [open Hits_to_arch.tsv w]
puts $out [string trim $arch_hits]
close $out
set out [open Hits_to_euk.tsv w]
puts $out [string trim $euk_hits]
close $out

###########################################################################

cd $direct/$org/$rbbh_dir

set i 0
set total_val 0
set fname "Master_rbh_weight.txt"
set f [open $fname]
while {[gets $f line] >= 0} {
	set vals [split $line \t]
	set eval [lindex $vals 2]
	if {$eval == 0.0} {
		set eval 1e-200
	}
	set total_val [expr $total_val + $eval]
	incr i
	puts $i
}
close $f

set av_eval [expr $total_val/$i]
set av_eval [format %.3g $av_eval]
set edge_no [expr $i / 2]
set edge_no [commas $edge_no]

set out [open $direct/$org/$rbbh_dir/Rbbh_stats/Edge_number_av_eval.txt w]
puts $out "Number of edges in the Master_rbh_weight file: $edge_no\nThe average e-val is: $av_eval"
close $out

###########################################################################

puts "======== DONE ========"