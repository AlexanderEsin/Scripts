#!/usr/local/bin/tclsh
###### CHEAT AND CHECK #######
proc openfile {fl} {
    global data
    set in [open $fl r]
    set data [read $in]
    close $in
    return
}

proc split_genes {fasta} {
    global genes
    regsub -all {>} [string trim $fasta] {£>} fasta
    set fasta [split $fasta £]
    regsub -all "{}" $fasta {} genes
    return $genes
}
###############################
set contig_file_split_size 25
set min_ident 85

###############################
set direct /users/aesin/desktop/Deer/Assembly
set cheat_op_dir $direct/Cheat_check
file mkdir $cheat_op_dir

set output_dir $direct/Bos_taurus/Alpha/Alpha_output_2
cd $output_dir

set best_sequences_l {}
set done_directories [lsort -dictionary [glob -type d *seq*]]
set number_directories_queried [llength $done_directories]

foreach done_directory $done_directories {
    cd $output_dir/$done_directory

    set best_seq_file [lindex [glob -nocomplain *.best.*] 0]
    if {[llength $best_seq_file] > 0} {
        openfile $best_seq_file

        set number [string range $best_seq_file [string last "\_" $best_seq_file]+1 [string first "\." $best_seq_file]-1]

        set sequences_l [split_genes [string trim $data]]
        if {[llength $sequences_l] > 100} {
            puts "Sequence number $number had too many sequences [llength $sequences_l] and is also likely to be a highly repetitive sequence ..."
            # set folder_to_move [glob -type d $output_dir/*_seq_$number]
            # file rename $folder_to_move $repeat_seq_dir

            continue
        } else {
            set success_counter 0
            foreach best_sequence $sequences_l {
                set header [string range $best_sequence 0 [string first \n $best_sequence]]
                set coverage [string range $header [string last "\_" $header]+1 [string last "\." $header]-1]
                if {$coverage > 5} {
                    set best_sequence [string trim $best_sequence]
                    lappend best_sequences_l $best_sequence
                    incr success_counter
                }
            }
            puts "Directory: $number\t===\t[llength $sequences_l]\t===\t$success_counter"
        }
        # set best_sequence [string trim [lindex $sequences_l 0]]
        # lappend best_sequences_l $best_sequence
    }
}
catch {set number_sequences_included [llength $best_sequences_l]}
set best_sequences_str [join $best_sequences_l \n]
set out [open $cheat_op_dir/Cheat_check.fasta w]
puts $out $best_sequences_str
close $out


cd $cheat_op_dir
catch {exec minimo Cheat_check.fasta -D FASTA_EXP=1 -D MIN_IDENT=$min_ident}

openfile [lindex [glob *fa] 0]
set number_of_contigs [regexp -all {>} [string trim $data]]
set contigs [split_genes $data]; set contig_counter [llength $contigs]

set output_file_number 1
while {$contig_counter > 0} {
    set output {}
    set contigs_written 0
    while {$contigs_written < $contig_file_split_size} {
        lappend output [lindex $contigs 0]
        set contigs [lrange $contigs 1 end]

        set contig_counter [llength $contigs]
        puts "Contig counter: $contig_counter"
        incr contigs_written
        puts "Contigs written out: $contigs_written"
    }
    set out [open Contigs_split_$min_ident\_$output_file_number\.fa w]
    puts $out [join $output ""]
    close $out
    incr output_file_number
}


puts "\nNumber of done directories queried:\t$number_directories_queried\nNumber of sequences included in the contig building:\t$number_sequences_included\nNumber of resulting contigs:\t$number_of_contigs\n"
