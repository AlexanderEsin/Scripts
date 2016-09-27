#!/usr/local/bin/tclsh
set direct /users/aesin/desktop/Deer/Assembly

set locus_name Alpha
set db_dir /users/aesin/desktop/Deer/Assembly/DB

set max_contig_inclusion 3

set max_target_seqs 1000

##########################

set round_no [lindex $argv 0]
if {$round_no == 0 || $round_no == ""} {
    set round_no 1
}
set next_round_no [expr $round_no + 1]
set round_name "_$round_no"

##########################

set input_dir $locus_name\_input$round_name
set output_dir $locus_name\_output$round_name
set done_dir $locus_name\_done_reads$round_name
set repeat_seq_dir $locus_name\_repeats$round_name
set minimus_dir $locus_name\_minimus$round_name
set velvet_log_dir Velvet_logs

###########################################################################

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

proc reverse_dict {lst} {
    global revdict
    set revdict {}
    set lst [lsort -dictionary $lst]
    while {[llength $lst] > 0} {
        set element [lindex $lst end]
        lappend revdict $element
        set lst [lrange $lst 0 end-1]
    }
    return $revdict
}

###########################################################################

cd $direct
set input_dir $direct/$input_dir
set output_dir $direct/$output_dir
set done_dir $direct/$done_dir
set repeat_seq_dir $direct/$repeat_seq_dir
set minimus_dir $direct/$minimus_dir
set velvet_log_dir $output_dir/$velvet_log_dir

file mkdir $output_dir
file mkdir $done_dir
file mkdir $repeat_seq_dir
file mkdir $minimus_dir
file mkdir $velvet_log_dir

##########################

cd $db_dir
set db_file [lindex [glob *db*] 0]
set db_name [string range $db_file 0 [string first "." $db_file]-1]

###########################################################################

cd $input_dir
set input_files [reverse_dict [glob *.fasta]]

foreach file $input_files {
    set output_name [string range $file 0 end-6]
    puts "Now making assembly based on $file ..."
    file mkdir $output_dir/$output_name

    ## Set max_memory (in GBs) below requested max - some users reported spiking above max argument given. 
    ## Default assembler is Velvet (documentation), so no need to specify
    ## Kmer size == 31
    ## Insert length == 270
    ## Exp-coverage is in kmers = i.e. 
    catch {exec /users/aesin/downloads/aTRAM-master/aTRAM_AE.pl -reads $db_dir/$db_name -target $file -output $output_dir/$output_name/Out_atram_$output_name -log_file $output_dir/$output_name/Log_$output_name\_assembly.txt -velvet_log $velvet_log_dir/Velvet_log_$file\.txt -kmer 31 -max_target_seqs $max_target_seqs -ins_length 270 -exp_coverage 8 -cov_cutoff 2 -max_memory 24 -max_processes 9 -iterations 1 >> $output_dir/Log_assemble_short.txt}

    set all_fasta_file [glob -nocomplain $output_dir/$output_name/*$output_name\.all.*]

    if {[llength $all_fasta_file] == 0} {
        set results_file [lindex [glob -nocomplain $output_dir/$output_name/*$output_name\.results.*] 0]
        openfile $results_file
        if {[regexp -all {no results} $data] == 1} {
            puts "No hits were found for this seeder sequence in the readset ==> continuing ..."
            file rename $file $done_dir
            continue
        } else {
            puts "ERROR"
            exit 2
        }
    } else {
        set $all_fasta_file [lindex $all_fasta_file 0]
    }

    openfile $all_fasta_file; set contig_no [llength [split_genes [string trim $data]]]

    if {$contig_no > $max_contig_inclusion} {
        puts "After one Atram iteration $file has been included in more than $max_contig_inclusion contigs ($contig_no) - we will treat this as a highly repetitive sequence. Moving it to $repeat_seq_dir ..."
        file rename $output_dir/$output_name $repeat_seq_dir

        file rename $file $done_dir

        continue

    } else {
        puts "After one Atram iteration $file has been included in fewer than 10 contigs ($contig_no) - running Atram for a full 5 iterations ... "
        catch {exec /users/aesin/downloads/aTRAM-master/aTRAM_AE.pl -reads $db_dir/$db_name -target $file -output $output_dir/$output_name/Out_atram_$output_name -log_file $output_dir/$output_name/Log_$output_name\_assembly.txt -velvet_log $velvet_log_dir/Velvet_log_$file\.txt -kmer 31 -max_target_seqs $max_target_seqs -ins_length 270 -exp_coverage 8 -cov_cutoff 2 -max_memory 24 -max_processes 9 -iterations 5 -start_iteration 1 >> $output_dir/Log_assemble_short.txt}

        file rename $file $done_dir
    }
}

###########################################################################

cd $output_dir

set best_sequences_l {}
set done_directories [glob -type d *seq*]

foreach done_directory $done_directories {
    cd $output_dir/$done_directory

    set best_seq_file [lindex [glob -nocomplain *.best.*] 0]
    if {[llength $best_seq_file] > 0} {
        openfile $best_seq_file

        set number [string range $best_seq_file [string last "\_" $best_seq_file]+1 [string first "\." $best_seq_file]-1]

        set sequences_l [split_genes [string trim $data]]
        if {[llength $sequences_l] > 100} {
            puts "Sequence number $number had too many sequences [llength $sequences_l] and is also likely to be a highly repetitive sequence. Moving it to $repeat_seq_dir ..."
            set folder_to_move [glob -type d $output_dir/*_seq_$number]
            file rename $folder_to_move $repeat_seq_dir

            continue
        }
        set best_sequence [string trim [lindex $sequences_l 0]]

        lappend best_sequences_l $best_sequence
    }
}

set best_sequences_str [join $best_sequences_l \n]
set out [open $minimus_dir/Good_contigs$round_name\.fasta w]
puts $out $best_sequences_str
close $out


cd $minimus_dir
catch {exec minimo Good_contigs$round_name\.fasta -D FASTA_EXP=1 -D MIN_IDENT=95}
set minimus_contigs_file [glob *fa]
openfile $minimus_contigs_file; set contigs [split_genes [string trim $data]]

cd $direct
set next_round_input_dir $direct/$locus_name\_input_$next_round_no
file mkdir $next_round_input_dir
set output_fasta_counter 1

foreach contig $contigs {
    set out [open $next_round_input_dir/$locus_name\_round_$next_round_no\_seq_$output_fasta_counter\.fasta w]
    puts $out [string trim $contig]
    close $out
    incr output_fasta_counter
}

puts "After round $round_no:\t[llength $contigs] contigs have been produced. They can now seed the next round of Atram assembly...."

