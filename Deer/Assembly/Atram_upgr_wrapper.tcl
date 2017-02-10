#!/usr/local/bin/tclsh

#### Packages ####
package require cmdline

#### OPTIONS ####
set species Bos_taurus

set parent /users/aesin/desktop/Deer
set direct $parent/Assembly/$species
set db_dir $parent/Assembly/Atram_DB

# Process the command line
set parameters {
    {locus_name.arg        ""   "Name of locus to use (String)"}
    {iteration.arg         1    "Iteration (Integer)"}
    {max_processes.arg     10    "Maximum number of processes (Integer)"}
    {max_target_seqs.arg   2000 "Maximum number of target sequences for blastn (Integer)"}
    {min_coverage.arg      1.5  "Minimum coverage for inclusion in Velvet assembly (Float)"}
    {max_contig_incl.arg   10   "Maximum number of contigs from iteration 1 of Atram (Integer)"}
    {contig_assembly.arg   1    "Perform the contig assembly using Minimus (Binary operator)"}
    {min_ident.arg         99   "Minimum identity for Minimus overlap (Integer)"}
    {contig_file_split.arg 25   "The number of contigs entries per file after Minimus for UCSC web-blat (Integer)"}
    {chromosome.arg        25   "Chromosome(s) onto which the contigs should map"}
}

array set arg [cmdline::getoptions argv $parameters]

# Verify required parameters
set requiredParameters {locus_name}
foreach parameter $requiredParameters {
    if {$arg($parameter) == ""} {
        puts stderr "Missing required parameter: -$parameter"
        exit 1
    }
}

# Displays the arguments
puts "Atram run with the following arguments:"
parray arg
puts ""

set locus_name              $arg(locus_name)
set round_no                $arg(iteration)
set max_processes           $arg(max_processes)
set max_target_seqs         $arg(max_target_seqs)
set min_coverage            $arg(min_coverage)
set max_contig_inclusion    $arg(max_contig_incl)
set contig_assembly_switch  $arg(contig_assembly)
set min_ident               $arg(min_ident)
set contig_file_split_size  $arg(contig_file_split)
set chromosome_map          $arg(chromosome)

###########################################################################
#####                           PROCS                                ######
###########################################################################
source ~/Documents/Scripts/General_utils.tcl

# Runs blat in current directory - must provide query, database arguments. Output parameter is optional -- all other blat options as default
# Output variable blat_output_file is the path to the output file
proc blat {query database {output blat_output.psl}} {
    global blat_output_file
    catch {exec blat $database $query $output}
    if {[file exists $output] == 1} {
        set blat_output_file [file join [pwd] $output]
    }
}
# Converts the blat output psl file into the UCSC fomat scoring system
proc pslToWebScore {psl_file {psl_to_web_output Psl_to_web_score_output.txt}} {
    global psl_to_web_output_path
    catch {exec pslScore $psl_file > $psl_to_web_output}
    set psl_to_web_output_path [file join [pwd] $psl_to_web_output]
}

# Input a fasta sequences(s) and calculate total sequence length for the sequence(s)
proc totalSeqLen {fasta_seq} {
    global total_seq_len
    set total_seq_len 0
    set seq_l [split_genes $fasta_seq]
    foreach seq $seq_l {
        # Remove the fasta header
        set no_header [string trim [string range $seq [string first \n $seq] end]]
        # Remove all new lines
        regsub -all "\n" $no_header {} no_newlines
        set seq_len [string length $no_newlines]
        set total_seq_len [expr $total_seq_len + $seq_len]
    }
    return $total_seq_len 
}

###########################################################################

##############################
## Process the round number ##
##############################

set next_round_no [expr $round_no + 1]
set round_name "_$round_no"

##########################
##  Set necessary dirs  ##
##########################

cd $direct
set input_dir $direct/$locus_name/$locus_name\_input$round_name
set output_dir $direct/$locus_name/$locus_name\_output$round_name
set done_dir $direct/$locus_name/$locus_name\_done_reads$round_name
set repeat_seq_dir $direct/$locus_name/$locus_name\_repeats$round_name
set blat_dir $direct/$locus_name/$locus_name\_blat_v_genome$round_name
set velvet_log_dir $output_dir/Velvet_logs

file mkdir $output_dir
file mkdir $done_dir
file mkdir $repeat_seq_dir
file mkdir $blat_dir
file mkdir $velvet_log_dir

##########################
##   Find and prep DB   ##
##########################

cd $db_dir
set db_file [lindex [glob *db*] 0]
set db_name [string range $db_file 0 [string first "." $db_file]-1]
set db_name "repeat_trim_ilv.fq"

###########################################################################
##              Main body of the atram assembly pipeline                 ##
###########################################################################

cd $input_dir
set input_files [reverse_dict [glob -nocomplain *.fasta]]

if {[llength $input_files] > 0} {

    foreach file $input_files {

        set output_name [string range $file 0 end-6]
        puts "\n\nNow making assembly based on $file ..."
        file mkdir $output_dir/$output_name

        ## Set max_memory (in GBs) below requested max - some users reported spiking above max argument given. 
        ## Default assembler is Velvet (documentation), so no need to specify
        ## Kmer size == 31
        ## Insert length == 270
        ## Exp-coverage is in kmers = i.e. 
        catch {exec /users/aesin/Documents/Scripts/Deer/aTRAM-master/aTRAM_AE.pl -reads $db_dir/$db_name -target $file -output $output_dir/$output_name/Out_atram_$output_name -log_file $output_dir/$output_name/Log_$output_name\_assembly.txt -velvet_log $velvet_log_dir/Velvet_log_$file\.txt -kmer 31 -max_target_seqs $max_target_seqs -ins_length 270 -exp_coverage 8 -cov_cutoff 2 -max_memory 24 -max_processes $max_processes -iterations 1 >> $output_dir/Log_assemble_short.txt}

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
            if {$round_no == 1} {
                puts "After one Atram iteration $file has been included in more than 10 contigs ($contig_no) - we will treat this as a highly repetitive sequence. Moving it to $repeat_seq_dir ..."
                # Move output directory to "repeat sequences"
                file rename $output_dir/$output_name $repeat_seq_dir
                # Move input file to done folder
                file rename $file $done_dir
            } else {
                puts "After one Atram iteration $file has been included in more than 10 contigs ($contig_no) - the expansion has moved into a highly repetitive sequence. Keep the original contig as the best contig."
                # Delete the output files produced after 1 iteration, but keep the log file
                set output_files [glob -nocomplain $output_dir/$output_name/Out_atram_*]
                foreach output_file $output_files {
                    file delete $output_file
                }
                # Take the contig sequence from the input file and save it as the .best.fasta output file to be included in the minimus assembly for this round
                openfile $file
                set original_contig [string trim $data]
                set out [open $output_dir/$output_name/Unchanged_$output_name\.best.fasta w]
                puts $out $original_contig
                close $out
                # Move input file to done folder
                file rename $file $done_dir
            }
            
            continue

        } else {
            puts "After one Atram iteration $file has been included in fewer than 10 contigs ($contig_no) - running Atram for a full 5 iterations ... "
            catch {exec /users/aesin/Documents/Scripts/Deer/aTRAM-master/aTRAM_AE.pl -reads $db_dir/$db_name -target $file -output $output_dir/$output_name/Out_atram_$output_name -log_file $output_dir/$output_name/Log_$output_name\_assembly.txt -velvet_log $velvet_log_dir/Velvet_log_$file\.txt -kmer 31 -max_target_seqs $max_target_seqs -ins_length 270 -exp_coverage 8 -cov_cutoff 2 -max_memory 24 -max_processes $max_processes -iterations 5 -start_iteration 1 >> $output_dir/Log_assemble_short.txt}

            file rename $file $done_dir
        }
    }
} else {
    puts stdout "Can't find any input files. Exiting..."
    exit 1
}

###########################################################################
###                      Optional contig assembly                        ##
###########################################################################

if {$contig_assembly_switch == 1} {
    ### Prepare the directories for the minimus contig assembly ###
    cd $direct
    set minimus_dir $direct/$locus_name/$locus_name\_minimus$round_name
    file mkdir $minimus_dir

    ### Go the output directory for this round and process all the *best* contig-containing fasta files ###
    cd $output_dir

    set best_sequences_l {}
    set done_directories [lsort -dictionary [glob -nocomplain -type d *$locus_name*]]
    if {[llength $done_directories] == 0} {
        set dir_search [string tolower $locus_name]
        set done_directories [lsort -dictionary [glob -type d *$dir_search*]]
    }
    set number_directories_queried [llength $done_directories]

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
            } else {
                set success_counter 0
                foreach best_sequence $sequences_l {
                    set header [string range $best_sequence 0 [string first \n $best_sequence]]
                    set coverage [string range $header [string last "\_" $header]+1 [string last "\." $header]-1]
                    if {$coverage > $min_coverage} {
                        set best_sequence [string trim $best_sequence]
                        lappend best_sequences_l $best_sequence
                        incr success_counter
                    }
                }
                #puts "Directory: $number\t===\t[llength $sequences_l]\t===\t$success_counter"
            }
        }
    }

    set number_sequences_included [llength $best_sequences_l]
    set best_sequences_str [join $best_sequences_l \n]
    set out [open $minimus_dir/Good_contigs$round_name\.fasta w]
    puts $out $best_sequences_str
    close $out

    puts "Now running Minimus ..."
    cd $minimus_dir
    catch {exec minimo Good_contigs$round_name\.fasta -D FASTA_EXP=1 -D MIN_IDENT=$min_ident}
    set minimus_contigs_file [lindex [glob *fa] 0]
    puts "--> Minimus done!"
    openfile $minimus_contigs_file
    set number_of_contigs [regexp -all {>} [string trim $data]]
    set contigs [split_genes $data]

    # Absolute path to contig file
    set path_to_contigs_file [file join [pwd] $minimus_contigs_file]
    # Path to file containing the whole genome
    set path_to_genome $parent/Genomes/$species/all_chroms.fa
    # Path to file translating chromosome number to GI
    set path_to_chr_no $parent/Genomes/$species/chr_NC_gi.txt

    openfile $path_to_chr_no
    set chr_transl [split [string trim $data] \n]
    set chr_of_interest [split [lsearch -inline -glob $chr_transl $chromosome_map*] \t]
    set chr_transl_gi "gi|[lindex $chr_of_interest 2]|ref|[lindex $chr_of_interest 1]|"
    puts "\n--> Chromosome $chromosome_map in $species corresponds to the following GI label: $chr_transl_gi"

    cd $blat_dir
    puts "\nRunning blat on minimus contigs vs $species genome..."
    blat $path_to_contigs_file $path_to_genome
    puts "Blat complete.\n\tOutput file == $blat_output_file"
    pslToWebScore $blat_output_file
    puts "\nBlat psl output converted to the UCSC web score format.\n\tOutput file == $psl_to_web_output_path"
    openfile $psl_to_web_output_path

    set blat_output_data [split [string trim $data] \n]
    foreach blat_hit $blat_output_data {
        set chr_mapped_to [lindex [split $blat_hit \t] 0]
        if {$chr_mapped_to == $chr_transl_gi} {
            set contig_name [lindex [split [lindex [split $blat_hit \t] 3] \:] 0]
            lappend mapped_correctly_l $contig_name
        }
    }

    set mapped_correctly_l [lsort -unique $mapped_correctly_l]
    set number_contigs_mapped_correctly [llength $mapped_correctly_l]
    puts "\nOut of the original $number_of_contigs contigs given by Minimus:\n\t--> $number_contigs_mapped_correctly contigs correctly mapped to Chromosome $chromosome_map"

    set checked_contigs {}
    foreach mapped_contig $mapped_correctly_l {
        lappend checked_contigs [lsearch -inline -glob $contigs >$mapped_contig*]
    }
    set checked_contigs_fasta [join $checked_contigs ""]
    set total_base_number [commas [totalSeqLen $checked_contigs_fasta]]
    puts "\nContigs cover a total of $total_base_number bases..."


    ### If the min_identity level is set to 99, take the output contigs and make them into the input sequences for the next iteration of this script ###
    if {$min_ident == "99"} {
        puts "\nMin_ident == 99 ---> Writing contigs as input for next round..."
        cd $direct
        set next_round_input_dir $direct/$locus_name/$locus_name\_input_$next_round_no
        file mkdir $next_round_input_dir
        set output_fasta_counter 1
        
        foreach contig $checked_contigs {
            set out [open $next_round_input_dir/$locus_name\_round_$next_round_no\_seq_$output_fasta_counter\.fasta w]
            puts $out [string trim $contig]
            close $out
            incr output_fasta_counter
        }
        puts "Contig writing done!"
    }

    ### Split the contigs into fasta files with 25 sequences per file for input into UCSC Blat ###
    puts "\nWriting contigs to fasta files containing 25 sequences each for UCSC blat web-version..."
    cd $minimus_dir
    set contig_counter [llength $checked_contigs]
    set output_file_number 1
    while {$contig_counter > 0} {
        set output {}
        set contigs_written 0
        while {$contigs_written < $contig_file_split_size} {
            lappend output [lindex $checked_contigs 0]
            set checked_contigs [lrange $checked_contigs 1 end]

            set contig_counter [llength $checked_contigs]
            #puts "Contig counter: $contig_counter"
            incr contigs_written
            #puts "Contigs written out: $contigs_written"
        }
        set out [open Contigs_split_$min_ident\_$output_file_number\.faa w]
        puts $out [join $output ""]
        close $out
        incr output_file_number
    }
    puts "Also done!s\n"


    puts "After round $round_no:\n\tNumber of done directories queried:\t$number_directories_queried\n\tNumber of sequences included in the contig building:\t$number_sequences_included\n\tNumber of resulting contigs:\t$number_of_contigs\n\tNumber of contigs mapping to Chromosome $chromosome_map:\t$number_contigs_mapped_correctly\n\tNumber of bases included in contigs:\t$total_base_number\n"
}

