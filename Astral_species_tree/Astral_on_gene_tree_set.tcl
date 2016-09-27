source ~/Dropbox/Scripts/Procs/combine_tree_files.tcl

## Set the working directory, the input and output subdirectories ##
set direct /users/aesin/Desktop/Geo_final_trees
set indiv_tree_dir $direct/No_duplicates_geo_mono

set astral_combined_input $direct/Astral/Astral_input
set astral_output $direct/Astral/Astral_output

## Set seed tree directory ##

set seed_dir $direct/Inside_seed_trees

## Make output directory, if it does not exist ##
file mkdir $astral_output
cd $direct

## Set the name for combined & seed tree files ##

set comined_tree_file No_duplicates_geomono.txt
set seed_tree_file Inside_seed_combined.txt

## If the individual gene trees have not yet been combined into a single file, do so now. Specify directory housing gene tree files (newick format), the extension to look for, and the name of the output file ##

if {[file exists $astral_combined_input/$comined_tree_file] == 0} {
	combine_tree_files $indiv_tree_dir .txt $astral_combined_input/$comined_tree_file
}

if {[file exists $astral_combined_input/$seed_tree_file] == 0} {
	combine_tree_files $seed_dir .txt $astral_combined_input/$seed_tree_file
}

## Input and output file names for Astral ##

set input_file $astral_combined_input/$comined_tree_file
set seed_file $astral_combined_input/$seed_tree_file

set output_file No_duplicates_geomono.txt

## Set stdout/stderr channel buffering to none to monitor live output ##

chan configure stdout -buffering none

## Run astral ##
exec >&@stdout java -Xms4g -Xmx30g -jar /users/aesin/desktop/software/Astral/astral.4.7.8.jar -i $input_file -o $astral_output/$output_file

## Optional - seed trees ##
# -e $seed_file


