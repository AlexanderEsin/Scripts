############	Geobacillus RNA_seq		############

## RNA-seq data file and the DL33 genbank file the reads were mapped to. 
## It's for G. thermoglucosidasius strain DL33 (from the Leak Lab in Bath) growing aerobically in a rich medium (TGP).

##	The transcript file has no headers so format is as follows:
#	Annotation_ID     CDS_sequence_number     CDS_ID    read_count    rpkm

## Ran the custom Make_fasta_gbk.tcl script to create a fasta file containing all the protein seqeunces that are labelled to match the RNA_seq rpkm file for x-reference --> DL33genbank_protein_CDS.faa

## The total number of protein sequences from the RNA_seq mapping		== 4,213
## Total number of proteins in refseq Geobacillus thermoglucosadius 	== 3,711

## Comparing the Geo_only (124) and Anoxy_geo_only (58) gene families frome the Geo_v_all/2.0 dataset

cd ~/desktop/Geo_v_all/2.0/
mkdir -p Geo_only_groups/Blast_v_Gthermo/Input/DB
mkdir -p Anoxy_geo_only_groups/Blast_v_Gthermo/Input/DB

## Manually copy the group fasta files for Geo_only/Anoxy_geo_only e.g. from Blast_v_viral in the same dir into the Input folders. The folder containing the group fastas should be called "Fastas"
## Manually copy the RNA_seq genbank fasta file into the DB directory (DL33genbank_protein_CDS.fasta --> Input/DB)

## Make blast DB for both sets of groups
cd ~/desktop/Geo_v_all/2.0/Geo_only_groups/Blast_v_Gthermo/Input/DB
makeblastdb -in *.faa -dbtype prot -out DL33genbank_protein_CDS

cd ~/desktop/Geo_v_all/2.0/Anoxy_geo_only_groups/Blast_v_Gthermo/Input/DB
makeblastdb -in *.faa -dbtype prot -out DL33genbank_protein_CDS

## Then split the input files (the Group fastas in the "Fastas" folder) so that we can parallelise Blastp
## Run the 2.Desk.Split_query_equal.tcl Script from Scripts/Geo_only/V_viral/ - change the various directory names in the script as appropriate. The default run should split the directory into 6 subfolders in "Split"

## The 3.HPC.Blastp_master.sh runs the Blast of all the groups against the G. thermoglucosidans mRNA-based protein sequence predictions..



