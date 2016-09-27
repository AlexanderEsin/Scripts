## Run the progressiveMauve command line. The binary is run out of the app directory - added to $PATH ##

cd /users/aesin/desktop/Geo_analysis/Geo_omes/Genome_alignment_2/Geobac_genomes_renamed
genome_files=*

ulimit -n 10000
## Trying to input the extracted geo tree - even with the names matching the files names of the genome files exactly - kept throwing an error. As the mauve guide-tree topology is very similar, I ran it without custom guide tree input ##
#progressiveMauve --scratch-path-1=../Mauve_scratch_1 --scratch-path-2=../Mauve_scratch_2 --output=../Mauve_multi_align/Geobac_19way --input-guide-tree=../Geobac_guide_tree/Species_tree_extracted_19_geo_topo.txt $genome_files

progressiveMauve --output=../Mauve_multi_align/Geobac_19way --disable-cache --output-guide-tree=../Geobac_guide_tree/19_way_geobac_guide_tree_mauve.txt $genome_files > ../Alignment_log.txt

## For reference - the guide tree translation is as such :

# 1		Geobacillus_caldoxylosilyticus.fna
# 2		Geobacillus_kaustophilus_HTA426.fna
# 3		Geobacillus_sp_C56_T3.fna
# 4		Geobacillus_sp_G11MC16.fna
# 5		Geobacillus_sp_G1w1.fna
# 6		Geobacillus_sp_GHH01.fna
# 7		Geobacillus_sp_JF8.fna
# 8		Geobacillus_sp_MAS1.fna
# 9		Geobacillus_sp_WCH70.fna
# 10	Geobacillus_sp_WSUCF1.fna
# 11	Geobacillus_sp_Y412MC52.fna
# 12	Geobacillus_sp_Y412MC61.fna
# 13	Geobacillus_sp_Y41MC1.fna
# 14	Geobacillus_stearothermophilus.fna
# 15	Geobacillus_thermocatenulatus.fna
# 16	Geobacillus_thermodenitrificans_NG80_2.fna
# 17	Geobacillus_thermoglucosidasius_C56_YS93.fna
# 18	Geobacillus_thermoleovorans_CCB_US3_UF5.fna
# 19	Geobacillus_vulcani.fna