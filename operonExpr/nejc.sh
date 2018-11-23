# Paper A: https://www.sciencedirect.com/science/article/pii/S0891584916303185
# Paper B: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5026706/

# SRA data downloaded from: https://www.ncbi.nlm.nih.gov/Traces/study/?acc=SRP067583
# Reference genome dowbloaded from: https://www.ncbi.nlm.nih.gov/nuccore/NC_007779.1/

# ------------------------ #
# Download all the individual SRA files based on list of SRAs downloaded from SRA data link above
/Users/aesin/Documents/Scripts/operonExpr/readOneLine.sh SRR_Acc_List.txt

# ------------------------ #
master_dir="/Users/aesin/Desktop/operonExpr"

# ------------------------ #
# Run fastqc on all individual SRA files
cd $master_dir/SRA_files/Individual
mkdir -p $master_dir/FastQC_indiv_raw

fastqc -o $master_dir/FastQC_indiv_raw -t 20 *gz

# Identify, in each set, files with no quality dip in middle read (see allCombined_fastQC)
# SRX1492014 : 275 / 276 / 279
# SRX1492015 : 280 / 281
# SRX1492016 : 284 / 285
# SRX1492017 : 288 / 289 / 292
# SRX1492018 : 293 / 294 / 296
# SRX1492019 : 297 / 298 / 301

# Make them big
gunzip *gz

# Combine read files based on experiment: SRX..14-16 is control, SRX..17-19 + nickel
cat SRR3033275.fastq SRR3033276.fastq SRR3033279.fastq > ../SRX1492014.fastq
cat SRR3033280.fastq SRR3033281.fastq > ../SRX1492015.fastq
cat SRR3033284.fastq SRR3033285.fastq > ../SRX1492016.fastq
cat SRR3033288.fastq SRR3033289.fastq SRR3033292.fastq > ../SRX1492017.fastq
cat SRR3033293.fastq SRR3033294.fastq SRR3033296.fastq > ../SRX1492018.fastq
cat SRR3033297.fastq SRR3033298.fastq SRR3033301.fastq > ../SRX1492019.fastq

# Make them small 
pigz *

# ------------------------ #
# Trim the raw read files with quality 30
cd $master_dir/SRA_files
trim_galore -o ../Trim_reads -q 30 *fastq

# ------------------------ #
# Run fastqc on trimmed combined files
cd ../Trim_reads
mkdir ../FastQC_comb_trim
fastqc --outdir ../FastQC_comb_trim -t 20 SRX*fq

# ------------------------ #
# Use the gffread function to convert ggf3 to gff file
gffread ../Genome/NC_007779.1.gff3 -T -o ../Genome/NC_007779.1.gtf

# ------------------------ #
# ------------------------ #
# Build STAR index
STAR --runThreadN 20 \
--runMode genomeGenerate \
--genomeDir ../Genome/STAR_index \
--genomeFastaFiles  ../Genome/NC_007779.1.fasta \
--sjdbGTFfile ../Genome/NC_007779.1.gtf \
--sjdbOverhang 49 \
--sjdbGTFfeatureExon CDS \
--genomeSAindexNbases 3


# STAR align the reads to the genome and produce count table
cd $master_dir/Trim_reads

trimReads_files=$(find . -name "*.fq" -print0 | xargs -0 ls)
for file in $trimReads_files
do
	dirName=$(basename $(echo $file) | sed 's/_trimmed.fq//g')
	mkdir -p ../Align/$dirName

	STAR --runThreadN 20 \
	--genomeDir ../Genome/STAR_index \
	--readFilesType Fastx \
	--readFilesIn $file \
	--alignIntronMax 1 \
	--quantMode TranscriptomeSAM GeneCounts \
	--outFileNamePrefix ../Align/$dirName/$dirName\_
done



















# # ------------------------ #
# # ------------------------ #
# # Build index
# cd $master_dir/Genome
# bowtie2-build NC_007779.1.fasta NC_007779.1_index

# # ------------------------ #
# # Alignment. % = aligned exactly once
# cd $master_dir/Trim_reads
# bowtie2 -x ../Genome/NC_007779.1_index -q ./SRX1492014_trimmed.fq -S ../Align/SRX1492014_aln.sam -p 20 # 97.00%
# bowtie2 -x ../Genome/NC_007779.1_index -q ./SRX1492015_trimmed.fq -S ../Align/SRX1492015_aln.sam -p 20 # 97.16%
# bowtie2 -x ../Genome/NC_007779.1_index -q ./SRX1492016_trimmed.fq -S ../Align/SRX1492016_aln.sam -p 20 # 97.15%
# bowtie2 -x ../Genome/NC_007779.1_index -q ./SRX1492017_trimmed.fq -S ../Align/SRX1492017_aln.sam -p 20 # 97.80%
# bowtie2 -x ../Genome/NC_007779.1_index -q ./SRX1492018_trimmed.fq -S ../Align/SRX1492018_aln.sam -p 20 # 97.53%
# bowtie2 -x ../Genome/NC_007779.1_index -q ./SRX1492019_trimmed.fq -S ../Align/SRX1492019_aln.sam -p 20 # 97.71


# # ------------------------ #
# # Select just mapped reads, sort and index
# cd $master_dir/Align
# mkdir -p $master_dir/Mapped
# samtools view -b -h -F 4 SRX1492014_aln.sam > ../Mapped/SRX1492014_mapped.bam
# samtools view -b -h -F 4 SRX1492015_aln.sam > ../Mapped/SRX1492015_mapped.bam
# samtools view -b -h -F 4 SRX1492016_aln.sam > ../Mapped/SRX1492016_mapped.bam
# samtools view -b -h -F 4 SRX1492017_aln.sam > ../Mapped/SRX1492017_mapped.bam
# samtools view -b -h -F 4 SRX1492018_aln.sam > ../Mapped/SRX1492018_mapped.bam
# samtools view -b -h -F 4 SRX1492019_aln.sam > ../Mapped/SRX1492019_mapped.bam


# samtools sort ../Mapped/SRX1492014_mapped.bam -o ../Mapped/SRX1492014_sorted.bam -@ 20
# samtools sort ../Mapped/SRX1492015_mapped.bam -o ../Mapped/SRX1492015_sorted.bam -@ 20
# samtools sort ../Mapped/SRX1492016_mapped.bam -o ../Mapped/SRX1492016_sorted.bam -@ 20
# samtools sort ../Mapped/SRX1492017_mapped.bam -o ../Mapped/SRX1492017_sorted.bam -@ 20
# samtools sort ../Mapped/SRX1492018_mapped.bam -o ../Mapped/SRX1492018_sorted.bam -@ 20
# samtools sort ../Mapped/SRX1492019_mapped.bam -o ../Mapped/SRX1492019_sorted.bam -@ 20

# samtools index ../Mapped/SRX1492014_sorted.bam
# samtools index ../Mapped/SRX1492015_sorted.bam
# samtools index ../Mapped/SRX1492016_sorted.bam
# samtools index ../Mapped/SRX1492017_sorted.bam
# samtools index ../Mapped/SRX1492018_sorted.bam
# samtools index ../Mapped/SRX1492019_sorted.bam

# # ------------------------ #
# # Use the gffread function to convert ggf3 to gff file
# gffread ../Genome/NC_007779.1.gff3 -T -o ../Genome/NC_007779.1.gff


# # ------------------------ #
# # Use the gffread function to convert ggf3 to gff file
# htseq-count -f "bam" -s "yes" SRX1492014_sorted.bam ../Genome/NC_007779.1.gff > SRX1492014_sorted.tab

# paste -d "\t" SRX1492014_sorted.tab SRX1492015_sorted.tab ...  > 


















