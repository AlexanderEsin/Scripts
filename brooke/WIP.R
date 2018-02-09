## Load packages
## If you don't have these packages, install them from the console using: install.packages("stringr") etc...
## For Biostrings & ShortRead see: http://bioconductor.org/packages/release/bioc/html/Biostrings.html and http://bioconductor.org/packages/release/bioc/html/ShortRead.html
## For ggseqlogo: devtools::install_github("omarwagih/ggseqlogo")

if (!require("pacman")) install.packages("pacman")
pacman::p_load("parallel", "devtools", "stringr", "dplyr", "ggplot2", "gplots", "reshape2", "seqinr", "ggpubr", "Biostrings", "ShortRead", "motifRG", "ggseqlogo", "wesanderson", "dagLogo", "BioSeqClass")

source("/Users/aesin/Documents/Scripts/brooke/seqAnalysis_functions.R")

### READ IN DATA ###
## Change the working directory to where you have your fastq files
## Get a list of fastq files from a given directory
fastq_file_path <- "/Users/aesin/Desktop/brooke"
fastq_file_list <- dir(path = fastq_file_path, pattern = "*merged2.fastq")

## Read the fastq files in. Lapply applies a function - over the elements in a list.
fastq_data_list <- lapply(fastq_file_list, function(fastq_file) {
	# the readFastq function does exactly what you think
	fastq_data <- readFastq(dirPath = fastq_file_path, pattern = fastq_file)
	return(fastq_data)
})

## Rename the elements in the fastq data list to match our P-set names
# Our file names have the format: "P[number]_merged2.fastq. Here we use a pattern search approach to just extract the P[number] part
# The pattern: [\\D]{1}[\\d]+ can be broken down as: [\\D] = any character (not digit), {1} = match previous pattern once, [\\d] = any digit, + = match previous pattern as many times as possible
# Many functions in R automatically apply over a vector: our vector here is three file names, so it automatically finds the pattern in all three elements of the vector
PSet_name_list <- str_extract(string = fastq_file_list, pattern = "([\\D]{1}[\\d]+)")
numPSets		<- length(PSet_name_list)
names(fastq_data_list) <- PSet_name_list


## From the fastq data list above, get just the DNA sequence and put this in a new list
## At the same time, calculate the number of raw reads - we'll want to compare our filtered data to this downstream
seq_data_list <- lapply(fastq_data_list, function(p_set) {
	seq_only <- sread(p_set)
	num_reads <- length(seq_only)
	# For each fastq dataset we return a list containing two objects - the reads and the number of the reads
	# We can name the list elements and the final output will be a nested list object (R loves this)
	return(list(seq_data = seq_only, num_raw_reads = num_reads))
})

## This data list inherits the names from the input (i.e. the P1, P4 etc
## designations) - so no need to rename again. To parse nested lists, we use the
## $ operator. E.g. to see how many reads we have in the P4 set:
print(seq_data_list$P4$num_raw_reads)


### FILTER READ DATA ###

## Here we combine several necessary filtering steps:
# 1. Reads must 324 nucleotides in length
# 2. Reads must contain the correct middle consensus sequence
# 3. Reads with Ns are removed

# Define the consensus DNA sequence we expect each read to contain
midConsensus	<- "GGAGATAGTACCTCGTCTACACCT"
consensus.seq	<- DNAString(midConsensus)

## Again, apply the operations over the data list
filtSeq_data_list <- lapply(PSet_name_list, function(p_set_name) {

	message(paste0("Filtering ", p_set_name, " ... "))
	# Use the $ operator to select the sequence data
	PSet_reads		<- seq_data_list[[p_set_name]]$seq_data
	PSet_num_raw	<- seq_data_list[[p_set_name]]$num_raw_reads
		
	# Filter by length. The square brackets are a way of indexing the entire set
	# of reads, e.g. PSet_reads[1] will give us the first read, PSet_reads[1:10]
	# will give us the first 10 reads. The which function asks the question which
	# indices of this object agree with the following statement. In this case:
	# which reads have 324 bases (width is the nucleotide length)
	len_filter <- PSet_reads[which(width(PSet_reads) == 324)]
	
	# Same concept as above, but the index is a little more complex.
	# vcountPattern(pattern = consensus.seq, subject = len_filter) asks whether a
	# read contains a particular sequence. If it contains it (which is what we
	# want) the function gives us a 1 - so we select for all "1"s
	cons_seq_filter <- len_filter[which(vcountPattern(pattern = consensus.seq, subject = len_filter) == 1)]

	# Extra step to make sure the consensus sequence appears exactly at positions 151-174
	cons_pos_seq_filter	<-  cons_seq_filter[which(unlist(startIndex(vmatchPattern(consensus.seq, cons_seq_filter))) == 151)]
	
	# Almost identical to above, but here if a read contains an "N", we get a "1",
	# so we select for all "0"s
	N_filter <- cons_pos_seq_filter[which(vcountPattern(pattern = DNAString(x = "N"), subject = cons_pos_seq_filter) == 0)]
	
	# Finally, we want to see how our filtering affects the number of reads at
	# each step. length(reads) gives us the number of reads. The data.frame()
	# function is a way of creating a table - a data frame in R. I define each
	# column name, and give each columns a value - the number of reads. The
	# stringsAsFactors statement at the end tells R to just treat the raw values
	# as raw values. The details are not too important.
	filter_df <- data.frame(RawReads = PSet_num_raw, LenFiltReads = length(len_filter), ConsFiltReads = length(cons_pos_seq_filter), NFiltReads = length(N_filter), stringsAsFactors = FALSE)
	
	# As above - for each P-set return the final filtered read set and a table of
	# the read numbers at each filter step
	return(list(filt_seq_data = N_filter, filter_nums = filter_df))
})
names(filtSeq_data_list)	<- PSet_name_list

# We can confirm that no Ns remain in these sequence sets by running e.g.:
# alphabetFrequency(filtSeq_data_list$P1$filt_seq_data, collapse = TRUE)

## We want to combine the data frames from individual P-sets into an overall table which we can then plot
# Extract the dataframes from the list, one by one - and add a column to the 
# table to name each row according to the P-set. Here we apply not over the
# object, but the number of elements in the list: this is because of stupid R
# behaviour where you can't extract the name of an object inside the lapply
# function. And we need to know which row correspond to P1 and which to P4 - etc...
filterNums_list <- lapply(1:length(filtSeq_data_list), function(p_set_index) {
	# use the index to get the P-set name
	p_set_name <- names(filtSeq_data_list)[p_set_index]
	
	# use the the index to get the data frame with the read numbers
	data_row <- filtSeq_data_list[[p_set_index]]$filter_nums
	
	# add a column to the data frane (P-set) so we know what numbers correspond to which set
	data_with_name <- cbind(PSet = p_set_name, data_row, stringsAsFactors = FALSE)
	return(data_with_name)
})

## This list contains elements - each is a row with the same columns. We bind all these rows together into a table
filterNums_df <- bind_rows(filterNums_list)
## ggplot is the best way to plot in R - but can be confusing. First, it
## requires a different data "shape". The melt function takes care of that you
## can see the difference by just calling the objects and seeing how the tables look
filterNums_melt <- melt(data = filterNums_df, id.vars = "PSet", variable.name = "FilterType", value.name = "NumReads")

## Call the plot function: the internal aes() function assigns what data we want where, the geom_bar() tells ggplot we want a bar chart. Assign it to a variable
filterNums_plot <- ggplot(data = filterNums_melt, aes(x = PSet, y = NumReads, fill = FilterType)) + geom_bar(position = "dodge", stat = "identity")

## We can call the variable to see the plot in the rstudio window
filterNums_plot


### UNIQUE READS AND CROSS-SET OVERLAP ###

## For each P-set - get only unique read sequences (appearing at least once)
uniqSeq_data_list <- lapply(filtSeq_data_list, function(p_set) {
	unique_reads <- unique(p_set$filt_seq_data)
	return(list(uniq_seq_data = unique_reads, UniqReads = length(unique_reads)))
})

## As above - create a data frame row with the P-set name and number of unique reads
uniqNums_list	<- lapply(1:length(uniqSeq_data_list), function(p_set_index) {
	p_set_name <- names(uniqSeq_data_list)[p_set_index]
	num_uniq_reads <- uniqSeq_data_list[[p_set_index]]$UniqRead
	data_with_name <- data.frame(PSet = p_set_name, UniqReads = num_uniq_reads, stringsAsFactors = FALSE)
	return(data_with_name)
})
uniqNums_df		<- bind_rows(uniqNums_list)

## Create a dataframe containing the number of filtered reads and unique reads for each P-set
# We merge the two data frames together by the PSet column - but we only want
# the names and the number of final filtered reads from the "filterNums_df" data
# frame. We take just those columns by using the [ ] index system. For data
# frames it is: df[rows,columns] - so we take all the rows (empty) and columns 1
# and 5 (c(1,5)).
filtVsUniq_df <- merge(filterNums_df[,c(1,5)], uniqNums_df, by = "PSet")

## Plot with ggplot as before
filtVsUniq_melt <- melt(data = filtVsUniq_df, id.vars = "PSet", variable.name = "ReadType", value.name = "NumReads")
filtVsUniq_plot <- ggplot(data = filtVsUniq_melt, aes(x = PSet, y = NumReads, fill = ReadType)) + geom_bar(position = "dodge", stat = "identity")


## Let's find how many unique sequences in any set appear in any other set
# Here we have two sapply function (basically like lapply above), iterating over
# each set of unique P-set reads. For the unique reads in each P-set we compare
# how many are identical to reads in other P-sets. The result is a matrix of all
# the comparisons e.g. P1 vs P4, P1 vs P8 etc...
crossSetOverlap_length_matrix <- sapply(uniqSeq_data_list, function(x) sapply(uniqSeq_data_list, function(y) {
	# Don't bother comparing the set to itself - waste of time. If the same,
	# produce an NA value in the matrix
	if (!identical(x$uniq_seq_data, y$uniq_seq_data)) {
	# intersect() gives us the reads that are identical between the sets
	# length() tells us how many are identical
	length(BiocGenerics::intersect(x$uniq_seq_data, y$uniq_seq_data))
	} else {
	NA
	}
}))

## Plot a heatmap of this data
# We use a recordPlot function here so that we can replay the plot at any point,
# not only when we call the function (i.e. here)
heatmap.2(x = crossSetOverlap_length_matrix, Rowv = FALSE, symm = TRUE, dendrogram = "none", trace = "none")
crossSetOverlap_length_hmap <- recordPlot()
plot.new()
crossSetOverlap_length_hmap


## Now we get a list of all sequences that appear in other sequences
# The same as above, but we don't calculate the length - here we want to keep
# the actual sequences
crossSetOverlap_matrix <- sapply(uniqSeq_data_list, function(x) sapply(uniqSeq_data_list, function(y) {
	if (!identical(x$uniq_seq_data, y$uniq_seq_data)) {
	BiocGenerics::intersect(x$uniq_seq_data, y$uniq_seq_data)
	} else {
	NA
	}
}))

## We use the reduce function here to collapse the matrix into a list of sequences - we lose the information on which two sets shared this sequence, but this is not important here
# The is.na() function selects all those matrix elements that are NA
# (self-comparisons) - the ! beforehand selects the inverse I.e. we do not
# include the NA value in our output list
allOverlap_seqs <- Reduce(c, c(crossSetOverlap_matrix[!is.na(crossSetOverlap_length_matrix)]))

# Each shared sequence is represented twice (e.g. from P1 vs P4 & from P4 vs P1)
# Using unique() we remove these duplicate sequences and also those
# sequences that are multiply present as they are shared between P1, P4 and P8
allOverUniq_seqs <- unique(allOverlap_seqs)

## In this lapply function we essentially check each set of original (not 
## unique) sequences for how many times these "shared" sequences appear within 
## that data. Output is a list for each P-set with the ID and number of times
## shared sequence appear in the P-set
seqCopyNum_list <- lapply(1:length(filtSeq_data_list), function(p_set_index) {
	# Name of the P-set we are comparing to
	p_set_name <- names(filtSeq_data_list)[p_set_index]
	# The non-unique set of sequences in this P-set
	pSetFilt_seqs <- filtSeq_data_list[[p_set_index]]$filt_seq_data
	# Select all indices of sequences that overlap (%in%) with our shared sequences set
	filtRepeat_seqs <- pSetFilt_seqs[which(pSetFilt_seqs %in% allOverUniq_seqs)]
	
	# Tabulate it - e.g. a shared sequence could be present one or 10,0000 times in the filtered P-set
	filtRepeat_df <- as.data.frame(table(filtRepeat_seqs), stringsAsFactors = FALSE)
	# Name the columns
	names(filtRepeat_df) <- c("Sequence", p_set_name)
	return(filtRepeat_df)
})

## Use reduce to combine all the list elements from above into a single data frame - this should have exactly the same number of rows as the number of shared sequences we identified
# If a sequence is not present in a P-set, we use an NA to fill the blank spot -
# then convert it to 0 in the second line
presenceAcrossPsets_df <- Reduce(function(x ,y) full_join(x, y, by = "Sequence"), seqCopyNum_list)
presenceAcrossPsets_df[is.na(presenceAcrossPsets_df)] <- 0

## Use the sweep() function to convert these raw numbers of occurence into proportions of all the times the sequence occurs:
# compare the output of: head(presenceAcrossPsets_df) and head(fractionAcrossPsets_df)
fractionAcrossPsets_df <- cbind(Sequence = presenceAcrossPsets_df[,1], round(sweep(presenceAcrossPsets_df[,-1], 1, rowSums(presenceAcrossPsets_df[,-1]), "/"), digits = 3), stringsAsFactors = FALSE)

## (Almost) finally we want to figure out the P-set in which the sequence is most represented - a "winner"
# For each sequence...
shared_seq_winner <- apply(fractionAcrossPsets_df[,-1], 1, function(proportions) {
	# ... figure out which P-set has the largest proportion of that sequence
	P_with_max <- which(proportions == max(proportions))
	# If two P-sets contain the same proportion of the sequence, return NA
	if (length(P_with_max) > 1) {
	return(NA)
	} else {
	# If there is a winner, is it by a big enough margin?
	# Here I arbitrarily used 90% but we can change this
	max_proportion <- proportions[P_with_max]
	# If the proportion  is large enough, return the P-set name of the winner
	if (max_proportion >= 0.9) {
		return(names(proportions[P_with_max]))
	# ... otherwise return NA
	} else {
		return(NA)
	}
	}
})

## Attach the winner identities as a column to the sequence ID and proportion table
# We can define a new column by simply coming up with a name after the $ sign
# and providing an object of the correct size
fractionAcrossPsets_df$winner <- shared_seq_winner

# Try head(fractionAcrossPsets_df) to see what the data looks like

## Super finally :) - we need to remove the shared sequences that are not in the winner P-set from all other P-sets
# (or form all P-sets if there's no winner - but obviously we can change this behaviour if you want)
finalFilter_data_list <- lapply(1:length(filtSeq_data_list), function(p_set_index) {
	# P_set name
	p_set_name <- names(filtSeq_data_list)[p_set_index]
	p_set_seq <- filtSeq_data_list[[p_set_index]]$filt_seq_data
	# All the sequences in the proportion shared data frame that we want to remove from this P-set
	# i.e. any rows where this P-set is NOT the winner
	to_remove_seqs <- fractionAcrossPsets_df$Sequence[-which(fractionAcrossPsets_df$winner == p_set_name)]
	
	# Remove every instance of those sequences
	trim_filtSeq_seqs <- p_set_seq[-which(p_set_seq %in% to_remove_seqs)]
 
	## Count the number of sequences after we remove the shared sequences (non-winners) and add do filtering table
	trim_filtSeq_num <- length(trim_filtSeq_seqs)
	update_data_frame <- cbind(filtSeq_data_list[[p_set_index]]$filter_nums, ShareFiltReads = trim_filtSeq_num)
	
	return(list(filt_seq_data = trim_filtSeq_seqs, filter_nums = update_data_frame))
})
# Rename final object with the P-set names
names(finalFilter_data_list) <- names(filtSeq_data_list)


## As at lines 97-126: plot the effect of filtering on the number of reads
finalFilterNums_list <- lapply(1:length(finalFilter_data_list), function(p_set_index) {
	p_set_name <- names(finalFilter_data_list)[p_set_index]
	data_row <- finalFilter_data_list[[p_set_index]]$filter_nums
	data_with_name <- cbind(PSet = p_set_name, data_row, stringsAsFactors = FALSE)
	return(data_with_name)
})

finalFilterNums_df <- bind_rows(finalFilterNums_list)
finalFilterNums_melt <- melt(data = finalFilterNums_df, id.vars = "PSet", variable.name = "FilterType", value.name = "NumReads")
finalFilterNums_plot <- ggplot(data = finalFilterNums_melt, aes(x = PSet, y = NumReads, fill = FilterType)) + geom_bar(position = "dodge", stat = "identity")

## We can call the variable to see the plot in the rstudio window
finalFilterNums_plot

finalUniqSeq_data_list <- lapply(finalFilter_data_list, function(p_set) {
	unique_reads <- unique(p_set$filt_seq_data)
	return(list(UniqSeqData = unique_reads, NumUniqReads = length(unique_reads)))
})

# sanity_check <- sapply(finalUniqSeq_data_list, function(x) sapply(finalUniqSeq_data_list, function(y) {
#   # Don't bother comparing the set to itself - waste of time. If the same,
#   # produce an NA value in the matrix
#   if (!identical(x$uniq_seq_data, y$uniq_seq_data)) {
#     # intersect() gives us the reads that are identical between the sets
#     # length() tells us how many are identical
#     length(BiocGenerics::intersect(x$uniq_seq_data, y$uniq_seq_data))
#   } else {
#     NA
#   }
# }))



# When translating and assinging codons, assume the bacterial code
bacterial_code	<- 11
BacGeneCode		<- getGeneticCode(as.character(bacterial_code))

## Prepare a "random" background proteome using the library assembly rules, unless one already exists
if (file.exists(file.path(fastq_file_path, "RandomDNA_sequences.rds"))) {
	message("Random DNA sequence file found - reading in...", appendLF = FALSE)
	randomDNA_bg	<- readRDS(file.path(fastq_file_path, "RandomDNA_sequences.rds"))
	message("\rRandom DNA sequence file found - reading in... done")
} else {
	message(paste0("Random DNA sequence file not found - creating and saving to", fastq_file_path, "..."), appendLF = FALSE)
	randomDNA_bg	<- makeRandSeqs(seqNumber = 500000)
	saveRDS(randomDNA_bg, file = file.path(fastq_file_path, "RandomDNA_sequences.rds"), compress = "gzip")
	message(paste0("Random DNA sequence file not found - creating and saving to", fastq_file_path, "... done"))
}

message("Translating random DNA background sequences...", appendLF = FALSE)
randomAA_bg				<- Biostrings::translate(randomDNA_bg, genetic.code = BacGeneCode)
message("\rTranslating random DNA background sequences... done")

PSet_name_list_wRand	<- c("RAND", PSet_name_list)
numPSets_wRand			<- length(PSet_name_list_wRand)



#### //// Nucleotide-level analysis //// ####


### /// GC-content differences between P-sets /// ###

## The expected GC-content per base in 324-long library prep
randGC_perPos		<- rowSums(t(consensusMatrix(randomDNA_bg, as.prob = T)[2:3,]))
## The actual GC-content per base between the P-sets
PSet_GCcont_list <- lapply(PSet_name_list, function(p_set_name) {
	# Calculate the frequencies of G and C nucleotides seperately for each sequence
	nucleot_seqs	<- finalUniqSeq_data_list[[p_set_name]]$UniqSeqData
	GC_perPos		<- rowSums(t(consensusMatrix(nucleot_seqs, as.prob = T)[2:3,]))
	# Normalise for the GC content based on the library prep
	GC_perPosNorm	<- GC_perPos - randGC_perPos
	# Make DF
	GC_perPos_df	<- data.frame(Position = seq(1:324), fractGC = GC_perPosNorm, stringsAsFactors = FALSE)
	names(GC_perPos_df)[2]	<- p_set_name
	return(GC_perPos_df)
})
names(PSet_GCcont_list)	<- PSet_name_list
PSet_GCcont_df		<- Reduce(function(df1, df2) merge(df1, df2, by = "Position"), PSet_GCcont_list)
PSet_GCcont_melt	<- melt(PSet_GCcont_df, id.vars = "Position", variable.name = "PSet", value.name = "deltaGC_Content")

PSet_GCcont_line	<- ggplot(data = PSet_GCcont_melt, aes(x = Position, y = deltaGC_Content, color = PSet)) +
	geom_smooth(span = 0.08, se = FALSE, method = "loess")

# http://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1001115

statComparisons		<- lapply(combn(unique(PSet_GCcont_melt$PSet), 2, simplify = FALSE), paste0)
PSet_GCcont_viol	<- ggplot(data = PSet_GCcont_melt, aes(x = PSet, y = deltaGC_Content, fill = PSet)) +
	geom_violin() +
	scale_fill_manual(values = wes_palette("Darjeeling", n = numPSets, type = "continuous")) +
	stat_compare_means(comparisons = statComparisons, paired = TRUE, method = "wilcox.test", p.adjust = "bonferroni")


### Codon usage bias

## We calculate codon usage using the Codon Adaptability Index for E. coli. The
## codons used in each (unique) sequence are compared to an index (the CAI) of
## optimal codons derived from highly expressed genes in E. coli (I think). We
## compare how the distributions of codons in each sequence in each P-set
## differs to the others.

# Read in Genscript CAI data
GS_cai_tbl	<- read.table(file.path(fastq_file_path, "GenScript_ecoli_caiTbl.txt"), sep = "\t", header = TRUE, colClasses = "character")
GS_cai_tbl	<- GS_cai_tbl[order(GS_cai_tbl$Triplet),]
gs_cai_indx <- as.numeric(GS_cai_tbl$Fraction)

# Pre-run some cai calculations that were part of the old cai() function
stops		<- which("Stp" == aaa(seqinr::translate(s2c(c2s(words())), numcode = bacterial_code)))
singulets   <- which(sapply(syncodons(words(), numcode = bacterial_code), length) == 1)
exclude		<- c(stops, singulets)
gs_cai_adj	<- gs_cai_indx[-exclude]
gs_cai_adj[gs_cai_adj < 1e-04] <- 0.01

## Calculate the CAI index for the random library
message("Calculating CAI for the Random DNA set (sampling 200,000 sequences) ... ")
randomDNA_bg_samp	<- sample(randomDNA_bg, 200000)
randomDNA_char		<- str_split(string = tolower(as.character(randomDNA_bg)), pattern = "")
CAIRandom_result	<- cai.faster_vf(randomDNA_char, w = gs_cai_adj, exclude = exclude)
CAIRandom_df		<- data.frame(PSet = rep("RAND", length(CAIRandom_result)), CAI.Value = CAIRandom_result, stringsAsFactors = FALSE)

# For each sequence in each P-set, calculate the CAI
message("WARNING: CAI calculation takes a while - about 1-2 minutes per P-set")
PSetCAI_genscr_list <- lapply(PSet_name_list, function(p_set_name) {
	message(paste0("Calculating CAI for P-set: ", p_set_name))
	
	# The cai() function takes character vectors split by each individual base (lower case)
	charSplit_seqs	<- str_split(string = tolower(as.character(finalUniqSeq_data_list[[p_set_name]]$UniqSeqData)), pattern = "")
	cai_result		<- cai.faster_vf(charSplit_seqs, w = gs_cai_adj, exclude = exclude)
	
	# Convert result to a dataframe structure
	cai_df	<- data.frame(PSet = rep(p_set_name, length(cai_result)), CAI.Value = cai_result, stringsAsFactors = FALSE)
	return(cai_df)
})

# Bind all the individual P-set results into a single dataframe
PSetCAI_genscr_df		<- bind_rows(c(list(CAIRandom_df), PSetCAI_genscr_list))
PSet_factors			<- unique(PSetCAI_genscr_df$PSet)
PSetCAI_genscr_df$PSet	<- factor(PSetCAI_genscr_df$PSet, levels =  c("RAND", sort(PSet_factors[PSet_factors != "RAND"])))

# Find the maximum CAI value here so we can plot the p-values at the right level
y_max		<- max(PSetCAI_genscr_df$CAI.Value)
pval_level	<- ceiling(y_max / 0.02) * 0.02

# Plot a boxplot of the CAI scores distributions and test difference
statComparisons		<- lapply(combn(unique(PSetCAI_genscr_df$PSet), 2, simplify = FALSE), paste0)
PSetCAI_genscr_plot	<- ggplot(data = PSetCAI_genscr_df, aes(x = PSet, y = CAI.Value, fill = PSet)) +
	geom_violin() +
	scale_fill_manual(values = rev(wes_palette("FantasticFox", n = numPSets_wRand, type = "continuous"))) +
	theme(panel.background = element_blank(), axis.line = element_line(color = "black")) +
	stat_compare_means(ref.group = "RAND", method = "t.test", label.y = pval_level)




### /// Positions of STOP codons within the sequences /// ###

# Define the stop codon
stop_AA		<- "*"

# For all the unique P-set sequences: translate, check for stop codons, get their locations
stop_AllPosits_list <- lapply(PSet_name_list, function(p_set_name) {
	# Get the (unique) nucloetide sequences for each P-set and translate
	nucleot_seqs	<- finalUniqSeq_data_list[[p_set_name]]$UniqSeqData
	transl_seqs		<- Biostrings::translate(nucleot_seqs, genetic.code = BacGeneCode)
	
	# Get all sequences that contain at least 1 stop codon
	stop_AllSeqs	<- transl_seqs[which(vcountPattern(stop_AA, transl_seqs) > 0),]

	# Count how many stop codons are in the dataset, and the number of sequences containing the stop codons
	stop_RawNum		<- sum(vcountPattern(stop_AA, transl_seqs)) / (length(transl_seqs) * mean(width(transl_seqs)))
	stop_InSeqNum	<- length(stop_AllSeqs) / length(transl_seqs)
	
	# For sequences with stop codons, calculate their position in the sequence
	stop_PositsFull	<- vmatchPattern(stop_AA, stop_AllSeqs)
	stop_Posits		<- unlist(stop_PositsFull)@start
	
	# Return a dataframe with two columns: the P-set name and the STOP positions
	stop_data_list	<- list(
		stop_posits_df = data.frame(PSet = rep(p_set_name, length(stop_Posits)), StopPosition = stop_Posits, stringsAsFactors = FALSE),
		stop_nums_df = data.frame(PSet = p_set_name, Raw.Num = stop_RawNum, InSeq.Num = stop_InSeqNum, stringsAsFactors = FALSE))
	return(stop_data_list)
})

# Bind all the individual P-set STOP count results into a single dataframe
stop_AllNums_df		<- bind_rows(lapply(stop_AllPosits_list, function(pset) return(pset$stop_nums_df)))
stop_AllNums_melt	<- melt(stop_AllNums_df, id.vars = "PSet", variable.name = "Count.Type", value.name = "Counts")

# Barplot to compare how many stop codons are present in the Pset as a
# proportion of all codons and how many sequences contains stop codons as a
# proportion of all sequences.
facet_labels		<- c(`Raw.Num` = "Proportion of codons that are STOP", `InSeq.Num` = "Proportion of sequences with at least 1 STOP")
stop_AllNums_plot	<- ggplot(data = stop_AllNums_melt, aes(x = PSet, y = Counts, fill = PSet)) +
	geom_bar(stat = "identity", position = "dodge") +
	scale_fill_manual(values = wes_palette("Cavalcanti", n = numPSets, type = "continuous")) +
	facet_wrap(~Count.Type, scales = "free", labeller = as_labeller(facet_labels)) + 
	theme(panel.background = element_blank(), axis.line = element_line(color = "black"), strip.background = element_rect(fill = "white"), strip.text = element_text(colour = "black", size = 10))


# Bind all the individual P-set STOP position results into a single dataframe
stop_AllPosits_df	<- bind_rows(lapply(stop_AllPosits_list, function(pset) return(pset$stop_posits_df)))
statComparisons		<- lapply(combn(unique(stop_AllPosits_df$PSet), 2, simplify = FALSE), paste0)
# Boxplot of the STOP position distribution across P-sets with pairwise significance tests
stop_AllPosits_plot	<- ggplot(data = stop_AllPosits_df, aes(x = PSet, y = StopPosition, fill = PSet)) +
	geom_violin(scale = "area") +
	ylim(0, 140) +
	scale_fill_manual(values = wes_palette("Darjeeling", n = numPSets, type = "continuous")) +
	geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 51, ymax = 58), color = NA, size = 0.5, fill = "white") +
	geom_hline(yintercept = 108, color = "red", linetype = "dashed") +
	geom_hline(yintercept = 51, color = "blue", linetype = "dashed") +
	geom_hline(yintercept = 58, color = "blue", linetype = "dashed") +
	ggtitle("Positions of stop codons in sequences containing stop codons") +
	theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "grey40")) +
	stat_compare_means(comparisons = statComparisons, method = "wilcox.test", p.adjust = "bonferroni")




### /// Positions of Shine-Dalgarno sequence in P-sets /// ###

# Define the strictest SD consensus sequence
SD_seq		<- DNAString("AGGAGGT")

# For each P-set (including Random) calculate how many SD sequences are present, how many sequences contain SD sites, and get a list of these SD sequences
SD_AllPosits_list <- lapply(PSet_name_list_wRand, function(p_set_name) {

	if (p_set_name != "RAND") {
		# Get the (unique) nucloetide sequences for each P-set
		nucleot_seqs	<- finalUniqSeq_data_list[[p_set_name]]$UniqSeqData
	} else {
		nucleot_seqs	<- randomDNA_bg	
	}
	
	# Get all sequences that contain at an SD sequence with max 0 nucleotide mismatch
	SD_AllSeqs		<- nucleot_seqs[which(vcountPattern(SD_seq, nucleot_seqs, max.mismatch = 0) > 0),]
	# Count how many SD sequences are in the dataset, and the number of sequences containing SD sequences
	SD_RawNum		<- sum(vcountPattern(SD_seq, nucleot_seqs, max.mismatch = 0))
	SD_InSeqNum		<- length(SD_AllSeqs) / length(nucleot_seqs)
	
	# For sequences with stop codons, calculate their position in the sequence
	SD_PositsFull	<- vmatchPattern(SD_seq, SD_AllSeqs, max.mismatch = 0)
	SD_Posits		<- unlist(SD_PositsFull)@start
	
	# Return a dataframe with two columns: the P-set name and the SD positions
	SD_data_list	<- list(
		SD_contain_seqs = SD_AllSeqs,
		SD_posits_df = data.frame(PSet = rep(p_set_name, length(SD_Posits)), SDPosition = SD_Posits, stringsAsFactors = FALSE),
		SD_nums_df = data.frame(PSet = p_set_name, Raw.Num = SD_RawNum, InSeq.Num = SD_InSeqNum, stringsAsFactors = FALSE))
	return(SD_data_list)
})
# Rename list and get the numberical SD stats (number and number of sequences containing SD)
names(SD_AllPosits_list)	<- PSet_name_list_wRand
SD_AllNums_df	<- bind_rows(lapply(SD_AllPosits_list, function(pset) return(pset$SD_nums_df)))


## Combine SD position list into single dataframe 
SD_AllPosits_df		<- bind_rows(lapply(SD_AllPosits_list, function(pset) return(pset$SD_posits_df)))
statComparisons		<- lapply(combn(unique(SD_AllPosits_df$PSet), 2, simplify = FALSE), paste0)

## Reorder so RAND comes first
PSet_factors			<- unique(SD_AllPosits_df$PSet)
SD_AllPosits_df$PSet	<- factor(SD_AllPosits_df$PSet, levels =  c("RAND", sort(PSet_factors[PSet_factors != "RAND"])))

# Boxplot of the SD position distribution across P-sets with pairwise significance tests
SD_AllPos_violin	<- ggplot(data = SD_AllPosits_df, aes(x = PSet, y = SDPosition, fill = PSet)) +
	geom_violin(scale = "area") +
	scale_fill_manual(values = rev(wes_palette("FantasticFox", n = numPSets_wRand, type = "continuous"))) +
	scale_y_continuous(breaks = round(seq(1, 324, by = 9), 1)) +
	geom_hline(yintercept = 324, color = "red", linetype = "dashed") +
	geom_hline(yintercept = 151, color = "blue", linetype = "dashed") +
	geom_hline(yintercept = 174, color = "blue", linetype = "dashed") +
	ggtitle("Positions of SD seqs in sequences containing SD seqs") +
	theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "grey40")) +
	stat_compare_means(ref.group = "RAND", method = "wilcox.test", label.y = 350)


# If we find a region of interest, e.g. the massive density spike in P1 at position 115,
# we can plot a seqLogo to see if there are any other interesting features nearby (visually)
# See the functions script for details...
P1_115_plot	<- plotSDSeqLogo("P1", 115)
P1_115_plot



### /// Peptide-level analysis /// ###

## Translate our per-set unique sequences
finalAA_data 	<- lapply(PSet_name_list, function(p_set_name) {
	nucleot_seqs	<- finalUniqSeq_data_list[[p_set_name]]$UniqSeqData

	transl_seqs		<- Biostrings::translate(nucleot_seqs, genetic.code = BacGeneCode)
	noStp_seqs		<- transl_seqs[which(vcountPattern("*", transl_seqs) == 0)]

	num_all			<- length(transl_seqs)
	num_noStp		<- length(noStp_seqs)
	print(c(num_all, num_noStp))
	return(list(AllSeqs = transl_seqs, NoStpSeqs = noStp_seqs, allAASeqNum = num_all, NoStpSeqNum = num_noStp))
})
names(finalAA_data) <- PSet_name_list


## Format the background proteome for the dagLogo function
dagProteome_bg	<- prepareProteome(fasta = randomAA_bg)

## Format the AA sequences for the dagLogo function
dagAAdata_list	<- lapply(PSet_name_list, function(p_set_name, bgProteome = dagProteome_bg) {
	message(paste0("Formatting input sequences for the enrichment analysis: ", p_set_name, "..."), appendLF = FALSE)

	noStopSeqs_char	<- as.character(finalAA_data[[p_set_name]]$NoStpSeqs)
	dagFormatSeqs	<- formatSequence_AE(seq = noStopSeqs_char, proteome = bgProteome)
	message(paste0("\rFormatting input sequences for the enrichment analysis: ", p_set_name, "... done"))

	return(dagFormatSeqs)
})
names(dagAAdata_list)	<- PSet_name_list

## Perform the dagLogo test
dagAAtest_list <- lapply(PSet_name_list, function(p_set_name, bgProteome = dagProteome_bg) {
	dagPeptideData	<- dagAAdata_list[[p_set_name]]
	
	message(paste0("Peforming enrichment analysis for: ", p_set_name))

	message("\tBuilding backgroud model...", appendLF = FALSE)
	built_bg_model	<- buildBackgroundModel_AE(dagPeptideData, proteome = bgProteome)
	message("\r\tBuilding backgroud model... done")
	
	message("\tAnalysing amino acid enrichment...", appendLF = FALSE)
	byAAEnrich_test			<- testDAU_AE(dagPeptideData, built_bg_model)
	message("\r\tAnalysing amino acid enrichment... done")

	message("\tAnalysing AA type enrichment...", appendLF = FALSE)
	byTypeEnrich_test		<- testDAU_AE(dagPeptideData, built_bg_model, group = "classic")
	message("\r\tAnalysing AA type enrichment... done")
	return(list(AAEnrich = byAAEnrich_test, typeEnrich = byTypeEnrich_test))
})
names(dagAAtest_list)	<- PSet_name_list


## // Hydropathy // ##

message("Calculating per-site hydrophobicity score for random set...", appendLF = FALSE)
randomAA_bg_samp		<- sample(randomAA_bg, 200000)
RAND_hydro_score		<- featureHydro_AE(as.character(randomAA_bg_samp))
RAND_hydro_perSeq_df	<- data.frame(PSet = rep("RAND", length(rowMeans(RAND_hydro_score))), SeqHydro = rowMeans(RAND_hydro_score), row.names = NULL, stringsAsFactors = FALSE)
message("\rCalculating per-site hydrophobicity score for random set... done")

hydroAAtest_list <- lapply(PSet_name_list, function(p_set_name) {
	message(paste0("Calculating per-site hydrophobicity score: ", p_set_name, "..."), appendLF = FALSE)
	PSet_hydro_score	<- featureHydro_AE(as.character(finalAA_data[[p_set_name]]$NoStpSeqs))
	rownames(PSet_hydro_score)	<- NULL
	
	PSet_hydro_perSeq_vec	<- rowMeans(PSet_hydro_score)
	PSet_hydro_perSeq_df	<- data.frame(PSet = rep(p_set_name, length(PSet_hydro_perSeq_vec)), SeqHydro = PSet_hydro_perSeq_vec, stringsAsFactors = FALSE)
	
	PSet_hydro_norm		<- colMeans(PSet_hydro_score) - colMeans(RAND_hydro_score)
	PSet_hydro_norm_mat	<- t(as.matrix(PSet_hydro_norm))

	rownames(PSet_hydro_norm_mat)	<- p_set_name
	colnames(PSet_hydro_norm_mat)	<- as.character(1:ncol(PSet_hydro_score))

	message(paste0("\rCalculating per-site hydrophobicity score: ", p_set_name, "... done"))

	return(list(PSet_hydro_perSite_norm = PSet_hydro_norm_mat, PSet_hydro_perSeq = PSet_hydro_perSeq_df))
})
names(hydroAAtest_list)	<- PSet_name_list
hydroAAtest_list$RAND$PSet_hydro_perSeq	<- RAND_hydro_perSeq_df


## Calculate and plot a heatmap comparing the hydropathy per site between the Psets
all_perSite_hydroNorm	<- do.call(rbind, lapply(hydroAAtest_list, function(p_set) {
	return(p_set$PSet_hydro_perSite_norm)
}))

perSite_hydroNorm_hmap	<- pheatmap::pheatmap(
    mat = all_perSite_hydroNorm * 10,
    color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdBu")))(100),
    cluster_rows = FALSE,
    cluster_cols = FALSE
)

## 
all_perSeq_hydro		<- bind_rows(lapply(hydroAAtest_list, function(x) return(x$PSet_hydro_perSeq)))
all_perSeq_hydro$PSet	<- factor(all_perSeq_hydro$PSet, levels = PSet_name_list_wRand)

pairComparisons			<- lapply(combn(PSet_name_list, 2, simplify = FALSE), paste0)

perSeq_hydro_violin		<- ggplot(data = all_perSeq_hydro, aes(x = PSet, y = SeqHydro, fill = PSet)) +
	geom_violin(scale = "area") +
	scale_fill_manual(values = rev(wes_palette("FantasticFox", n = numPSets_wRand, type = "continuous"))) +
	scale_y_continuous(breaks = seq(-0.4, 0.2, by = 0.05)) +
	ggtitle("Per-sequence hydropathy in P-Sets") +
	theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "grey40")) +
	stat_compare_means(comparisons = pairComparisons, method = "wilcox.test", p.adjust = "bonferroni") + 
	stat_compare_means(ref.group = "RAND", method = "wilcox.test", label.y.npc = "bottom")



## // RNA folding // ##


numCores	<- detectCores() - 2
randSample	<- 50000

RNAfold_path	<- file.path(fastq_file_path, "RNAfold")
dir.create(RNAfold_path, showWarnings = FALSE)

RNAfold_full_list	<- lapply(PSet_name_list_wRand, function(p_set_name) {

	message(paste0("Performing RNA fold analysis for PSet: ", p_set_name))

	## Split data into chunks based on number of cores ##
	if (p_set_name == "RAND") {
		message(paste0("\tFor the random background, sequences sampled = ", randSample))
		PSet_nucleot	<- sample(randomDNA_bg, randSample)
	} else {
		PSet_nucleot	<- finalUniqSeq_data_list[[p_set_name]]$UniqSeqData
	}

	message(paste0("\tSplitting DNA data into ", numCores, " chunks..."), appendLF = FALSE)
	PSet_chunks			<- chunk(PSet_nucleot, numCores)
	names(PSet_chunks)	<- paste("chunk", 1:numCores, sep = "")
	message(paste0("\r\tSplitting DNA data into ", numCores, " chunks... done"))

	## Create the output directory to store RNAfold predictions
	p_set_path	<- file.path(RNAfold_path, p_set_name)
	dir.create(p_set_path, showWarnings = FALSE)

	message(paste0("\tCalculating RNA folding on ", numCores, " cores (this can take a while)..."), appendLF = FALSE)

	clust	<- makeCluster(numCores, type = "FORK")
	outputFiles	<- parLapply(clust, 1:length(PSet_chunks), function(chunk_index) {
		chunk_name		<- names(PSet_chunks)[chunk_index]

		chunk_fileName	<- paste0(chunk_name, "_RNAfold_in.fa")	
		chunk_inPath	<- file.path(p_set_path, chunk_fileName)
		chunk_outPath	<- file.path(p_set_path, paste0(chunk_name, "_RNAfold_out.txt"))

		if (!file.exists(chunk_outPath)) {
			chunk_DNA		<- PSet_chunks[[chunk_index]]
			chunk_RNA		<- RNAStringSet(chunk_DNA)
			names(chunk_RNA)	<- 1:length(chunk_RNA)

			writeXStringSet(chunk_RNA, file = chunk_inPath)
			system2('RNAfold', args = c("--filename-delim=/", '--noPS', paste0("--outfile=", chunk_outPath), paste0("--infile=", chunk_inPath)))
		}
		return(chunk_outPath)
	})
	stopCluster(clust)
	message(paste0("\r\tCalculating RNA folding on ", numCores, " cores (this can take a while)... done"))

	PSet_foldFull	<- bind_rows(lapply(outputFiles, function(file) struct <- readRNAfold(file)))

	PSet_deltaGs	<- unlist(lapply(outputFiles, function(file) deltaG <- readRNAfold(file)$deltaG))
	PSet_deltaGs_df	<- data.frame(PSet = rep(p_set_name, length(PSet_deltaGs)), foldDeltaG = PSet_deltaGs, stringsAsFactors = FALSE)

	return(list(FullFold = PSet_foldFull, deltaGs = PSet_deltaGs_df))
})
names(RNAfold_deltaG_list)	<- PSet_name_list_wRand

RNAfold_deltaG_df		<- bind_rows(lapply(RNAfold_full_list, function(x) return(x$deltaGs))
RNAfold_deltaG_df$PSet	<- factor(RNAfold_deltaG_df$PSet, levels = PSet_name_list_wRand)

pairComparisons			<- lapply(combn(PSet_name_list, 2, simplify = FALSE), paste0)
RNAfold_deltaG_violin	<- ggplot(data = RNAfold_deltaG_df, aes(x = PSet, y = foldDeltaG, fill = PSet)) +
	geom_violin(scale = "area") +
	scale_fill_manual(values = rev(wes_palette("FantasticFox", n = numPSets_wRand, type = "continuous"))) +
	# scale_y_continuous(breaks = seq(-0.4, 0.2, by = 0.05)) +
	ggtitle("Per-sequence RNA-folding free energy across P-sets") +
	theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "grey40")) +
	stat_compare_means(comparisons = pairComparisons, method = "wilcox.test", p.adjust = "bonferroni") + 
	stat_compare_means(ref.group = "RAND", method = "wilcox.test", label.y.npc = "bottom")






























