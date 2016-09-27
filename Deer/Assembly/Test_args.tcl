#!/usr/local/bin/tclsh
package require cmdline
 
# Show argv before processing
puts "Before, argv = '$argv'"

# Process the command line
set parameters {
    {locus_name.arg        ""   "Name of locus to use"}
    {max_processes.arg     1    "Maximum number of processes"}
    {max_target_seqs.arg   2000   "Maximum number of target sequences for blastn"}
    {min_ident.arg         99   "Minimum identity for Minimus overlap"}
    {min_coverage.arg      1.5  "Minimum coverage for inclusion in Velvet assembly"}
    {debug                      "Output extra debug info"}
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
puts ""
parray arg
puts ""

# Show argv after processing
puts "After, argv = '$argv'"
puts "$arg(locus_name)"