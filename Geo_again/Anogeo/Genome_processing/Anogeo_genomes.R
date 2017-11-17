library(dplyr)

## 12/11/2017
## Filter out all the "complete genomes" from the Bacterial Refseq db: ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/
## Download assembly_summary_refseq.txt from above url: wget ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/bacteria/assembly_summary.txt; mv assembly_summary.txt assembly_summary_refseq_13117.txt
## Manually remove the hashtag infront of columns names (line 2)

setwd("/users/aesin/Desktop/Geo_again/Genomes/")
dir.create("Anogeo/Anogeo_selection")

# We get 100,270 organism entries
bact_genome_tbl <- read.table("assembly_summary_refseq_13117.txt", sep = "\t", as.is = TRUE, header = TRUE, fill = TRUE, quote = "")

# Semantically find all Geobacillus and Anoxybacillus taxids
all_geobac_entries	<- bact_genome_tbl[which(grepl(pattern = "geobacillus", ignore.case = T, bact_genome_tbl$organism_name) == TRUE),]
all_anobac_entries	<- bact_genome_tbl[which(grepl(pattern = "anoxybacillus", ignore.case = T, bact_genome_tbl$organism_name) == TRUE),]
all_anogeo_entries	<- rbind(all_geobac_entries, all_anobac_entries)

# A total of 116 genomes representing Geobacillus and Anoxybacillus
total_entries	<- nrow(all_anogeo_entries)
# A total of 87 unique taxids
unique_anogeo	<- unique(all_anogeo_entries$taxid)

assembly_levels	<- c("Complete Genome", "Chromosome", "Scaffold", "Contig")

best_assembled_anogeo <- lapply(unique_anogeo, function(taxid) {
	subset_data	<- all_anogeo_entries[which(all_anogeo_entries$taxid == taxid),]

	for (ass_level in assembly_levels) {
		assembled_genomes	<- subset_data[which(subset_data$assembly_level == ass_level),]

		if (nrow(assembled_genomes) == 1) {
			picked_assembly <- assembled_genomes
			break
		} else if (nrow(assembled_genomes) > 1) {
			latest_date		<- max(as.Date(subset_data$seq_rel_date))
			picked_assembly	<- subset_data[which(subset_data$seq_rel_date == latest_date),]
			if (nrow(picked_assembly) > 1) {
				picked_assembly	<- picked_assembly[sample(1:nrow(picked_assembly), 1),]
			}
			break
		}
	}

	return(picked_assembly)
})

best_assembled_anogeo <- bind_rows(best_assembled_anogeo)
names(best_assembled_anogeo)	<- names(bact_genome_tbl)
write.table(x = best_assembled_anogeo, file = "./Anogeo/Anogeo_selection/best_assembled_anogeo.tsv", sep = "\t", quote = FALSE, row.names = FALSE)

## Then picked a maximum of two genomes per annotated species. Complete genomes were preferred, otherwise by best assembly (chromosome, scaffold, contig) and by date.
## See Anogeo_strain_select.xlsx, yellow highlighted entries are not included --> output ./Anogeo_selection/Assembled_anogeo_trimmed.tsv

# Read in the manually trimmed Anogeo genome info
anogeo_curated <- read.table("./Anogeo/Anogeo_selection/Assembled_anogeo_trimmed.tsv", sep = "\t", as.is = TRUE, header = TRUE, fill = TRUE, quote = "")
dir.create("./Anogeo/Anogeo_genomes_raw/Genomic_gbff")
dir.create("./Anogeo/Anogeo_proteomes/Proteome_fasta_raw")

temp	<- lapply(1:nrow(anogeo_curated), function(row_index) {
	genome_info	<- anogeo_curated[row_index,]
	genome_name	<- genome_info$organism_name
	ftp_path	<- genome_info$ftp_path

	acc_ass		<- basename(ftp_path)
	file_name	<- paste0(acc_ass, "_genomic.gbff.gz")
	dl_file		<- paste0(ftp_path, "/", file_name)
	out_file	<- paste0(getwd(), "./Anogeo/Anogeo_genomes_raw/Genomic_gbff/", file_name)

	download.file(dl_file, out_file, method = "wget", quiet = TRUE)
	return("DL_done")
})
rm(temp)

## Convert the gbff files into proteome fasta files in ./Genomes/Anogeo_proteomes/Proteome_fasta_raw
system2("/Users/aesin/Documents/Scripts/Geo_again/Anogeo_proteome_convert.sh")



