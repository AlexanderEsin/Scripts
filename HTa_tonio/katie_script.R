library(stringr)
# Take all annotated sequences
#histone_aligned_annot <- getAnnot(histone_aligned_fasta_length)
# Remove seed from annotations
#histone_aligned_seed <- gsub('seed', '', histone_aligned_annot)
# Regex to remove all seed underscores by taing everything after the first letter/number
#histone_aligned_annot <- str_match(histone_aligned_seed, "[A-Za-z0-9]+.+")
# Regex to remove the alignment positions following name WP_xxxx/69-67 basically after the /
#raw_fasta_headers <- str_replace(names(histone_aligned_fasta_length), "/+.+", "")
# Make into df
#all_fasta_names <- data.frame(header = raw_fasta_headers, stringsAsFactors = FALSE)
### NOT UNIQUEEEEEE GO BACK AND SORT THIS OUT! PROBABLY FROOM DOUBLETS

# Take all annotated sequences
histone_aligned_annot <- getAnnot(histone_aligned_fasta_length)
# Regex to remove the alignment positions following name WP_xxxx/69-67 basically after the /
raw_fasta_headers <- str_replace(names(histone_aligned_fasta_length), "/+.+", "")
# Remove seed from annotations
histone_aligned_seed <- gsub('seed', '', raw_fasta_headers)
# Regex to remove all seed underscores by taing everything after the first letter/number
histone_aligned_annot_comp <- str_match(histone_aligned_seed, "[A-Za-z0-9]+.+")

head(histone_aligned_annot_comp)
# Make into df
all_fasta_names <- data.frame(header = histone_aligned_annot_comp, stringsAsFactors = FALSE)


## Combine all the labelling into one dataframe
arch_labels_raw <- getAnnot(archaeal_histone_monomers)
# split the archaea by the first bracket ">No_annot [Methanosphaera stadtmanae DSM 3091]"
arch_split_by_brack <- str_split(arch_labels_raw, "\\[")

# Process species names: Take the second part (i.e. after bracket above) and remove the last bracket  ">WP_014965732.1  MULTISPECIES  histone 1027373 " "Nitrosopumilus sp. AR]"
species_names <- unlist(lapply(arch_split_by_brack, function(header) {header[2]}))
species_names <- str_replace(species_names, "\\]", "")

# The rest
# Take everything before the bracket (the rest)
other_crap <- unlist(lapply(arch_split_by_brack, function(header) {header[1]}))
# Make df of the stuff before the bracket
name_func_df <- as.data.frame(str_split_fixed(other_crap, " ", 2), stringsAsFactors = FALSE)
# Make Entry and Protein.names columns
names(name_func_df) <- c("Entry", "Protein.names")
# In the Entry column remove > (replace with space)
name_func_df$Entry <- str_replace(name_func_df$Entry , ">+", "")
# Put species names into organism column
name_func_df$Organism <- species_names

# Combine the euk and arch label dfs
# Get the uniprot things of interest for the eukaryotic histones i.e. entry, protien name and organism name
uniprot_table <- select(uniprot, "Entry", "Protein.names", "Organism")
all_label_df <- rbind(uniprot_table, name_func_df)


## Replace fasta names with proper labels
new_fasta_names <- left_join(all_fasta_names, all_label_df, by = c("header" = "Entry"))

# Format the header into one line separated by tab
formattedHeader <- unlist(lapply(1:nrow(new_fasta_names), function(index) {
  # Get the index of the labels and call header
  header <- as.character(new_fasta_names[index,])
  # paste the header and collapse by tab
  header_format <- paste(header, collapse = "\t")
  # paste > infront of the header
  header_format <- paste0(">", header_format)
  # header_format is output
  return(header_format)
}))

# Put back into original fasta
histone_aligned_fasta_maybe <- histone_aligned_fasta_length
names(histone_aligned_fasta_maybe) <- formattedHeader

# Write out complete fasta file with names
write.fasta(sequences = histone_aligned_fasta_maybe, names = names(histone_aligned_fasta_maybe), file.out = file.path(mafft_linsi_realigned_dir, "histone_aligned_names_complete.fa"))

