#!/usr/local/bin/tclsh
set direct /users/aesin/desktop/Geo_v_all/2.0
set group_type Anoxy_geo_only_groups

## Viral proteome location ##
set viral_loc /users/aesin/desktop/Viral_files
## If there is no DB folder at $viral_loc - what is the name of the combined proteome file? ##
set viral_prot_comp Viral_complete.faa

###########################################################################

cd $direct/$group_type
file mkdir Blast_v_viral

cd $direct/$group_type/Blast_v_viral
file mkdir Input
file copy -force $direct/$group_type/Group_fastas $direct/$group_type/Blast_v_viral/Input
file rename $direct/$group_type/Blast_v_viral/Input/Group_fastas $direct/$group_type/Blast_v_viral/Input/Fastas

cd $viral_loc

if {[file exists DB] == 1} {
	file copy -force DB $direct/$group_type/Blast_v_viral/Input
} elseif {[file exists $viral_prot_comp ] == 1} {
	set input_file $viral_prot_comp

	catch {exec makeblastdb -in $input_file -dbtype prot -out $viral_loc/DB/$input_file.db}
	file copy -force $viral_loc/DB $direct/$group_type/Blast_v_viral/Input
} else {
		"ERROR: Cannot find the proteome file!"
		exit 2
}

puts "===DONE==="

###########################################################################

