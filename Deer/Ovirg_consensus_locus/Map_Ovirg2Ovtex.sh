## Download the Ovtex scaffold containing the HB locus - NW_018335738.1
cd /users/aesin/Desktop/Deer/Assembly/Ovirg_scaffold_0917
mv NW_018335738.1[1..380642].fa Ovtex_HB_scaff.fa

## Also download the genbank flat file of the NW_018335738.1 scaffold (for annotations)
# Ovtex_HB_scaff.flat

# Make index based on the Ovtex scaffold fasta
mkdir Ovtex_index
bowtie2-build Ovtex_HB_scaff.fasta ./Ovtex_index/Ovtex_HB

# We want to include mixed reads because of cases where one mate is in a non-Ovtex represented region
bowtie2 --local -p 20 -x ./Ovtex_index/Ovtex_HB -1 ./Reads/Ovirg_paired_R1.fastq -2 ./Reads/Ovirg_paired_R2.fastq -S ./Mapped/Ovirg2Ovtex_local.sam

# 163656267 reads; of these:
#   163656267 (100.00%) were paired; of these:
#     153083304 (93.54%) aligned concordantly 0 times
#     2110859 (1.29%) aligned concordantly exactly 1 time
#     8462104 (5.17%) aligned concordantly >1 times
#     ----
#     153083304 pairs aligned concordantly 0 times; of these:
#       140248 (0.09%) aligned discordantly 1 time
#     ----
#     152943056 pairs aligned 0 times concordantly or discordantly; of these:
#       305886112 mates make up the pairs; of these:
#         296295256 (96.86%) aligned 0 times
#         2375321 (0.78%) aligned exactly 1 time
#         7215535 (2.36%) aligned >1 times
# 9.48% overall alignment rate
# For comaprison - without "local"
# 3.21% overall alignment rate


# Extract all mapped reads (regardless of mapped mate)
samtools view -@ 20 -h -b -F 4 -q 1 ./Mapped/Ovirg2Ovtex_local.sam > ./Mapped/Ovirg2Ovtex_onemap_Q1.bam
# Name sort
samtools sort -n -@ 20 -o ./Mapped/Ovirg2Ovtex_onemap_namesort.bam ./Mapped/Ovirg2Ovtex_onemap_Q1.bam
# Add ms and MC tags for markdup to use later
samtools fixmate -m ./Mapped/Ovirg2Ovtex_onemap_namesort.bam ./Mapped/Ovirg2Ovtex_onemap_fixmate.bam
# Coordinate sort
samtools sort -@ 20 -o ./Mapped/Ovirg2Ovtex_onemap_coordsort.bam ./Mapped/Ovirg2Ovtex_onemap_fixmate.bam
# Mark + remove duplicates. Print some stats with -s
samtools markdup -@ 20 -r -s ./Mapped/Ovirg2Ovtex_onemap_coordsort.bam ./Mapped/Ovirg2Ovtex_onemap_remdup.bam

# READ 30051897 WRITTEN 9042985 
# EXCLUDED 0 EXAMINED 30051897
# PAIRED 22736688 SINGLE 7315209
# DULPICATE PAIR 13707882 DUPLICATE SINGLE 7301030
# DUPLICATE TOTAL 21008912

# Read file into Geneious so that we can take the consensus: Local/Ovirg_ovtex_assembly
# read in:					Ovirg2Ovtex_onemap_remdup - Ov_tex_HB_scaffold
# make consensus:			Ovirg_locus_consensus_sequence
# extract consensus:		Ovirg_locus_consensus_sequence_extract
# LastZ map to Ovtex:		NW_018335738_LASTZ_alignment
# Alignment: ref + cons:	Ovirg_ovtex_scaffold_alignment

# Write out the alignment (last file) to folder:	Ovtex_Ovirg_alignment.fasta
# Write out per-nucl identity score (3dp):			Ovtex_Ovirg_alignment_ID.csv

## Proceed to R script for processing
## Zip all mapping files (delete initial SAM file). Delete unzipped copy of reads. Originals in Deer/Ovirg_mapping_initial/Reads/Zipped_read_files
