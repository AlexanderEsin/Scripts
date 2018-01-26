## Load packages
## If you don't have these packages, install them from the console using: install.packages("stringr") etc...
## For Biostrings & ShortRead see: http://bioconductor.org/packages/release/bioc/html/Biostrings.html and http://bioconductor.org/packages/release/bioc/html/ShortRead.html
## For ggseqlogo: devtools::install_github("omarwagih/ggseqlogo")
library(Biostrings)
library(ShortRead)
library(stringr)
library(dplyr)
library(ggplot2)
library(gplots)
library(reshape2)
library(seqinr)
library(ggpubr)
library(devtools)
if("ggseqlogo" %in% rownames(installed.packages()) == FALSE) {
	devtools::install_github("omarwagih/ggseqlogo")
}
library(ggseqlogo)

## AE: look into DNAShapeR

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
consensus.seq <- DNAString(x = "GGAGATAGTACCTCGTCTACACCT")

## Again, apply the operations over the data list
filtSeq_data_list <- lapply(seq_data_list, function(p_set) {
	# Use the $ operator to select the sequence data
	reads_seq <- p_set$seq_data
	
	# Filter by length. The square brackets are a way of indexing the entire set
	# of reads, e.g. reads_seq[1] will give us the first read, reads_seq[1:10]
	# will give us the first 10 reads. The which function asks the question which
	# indices of this object agree with the following statement. In this case:
	# which reads have 324 bases (width is the nucleotide length)
	len_filter <- reads_seq[which(width(reads_seq) == 324)]
	
	# Same concept as above, but the index is a little more complex.
	# vcountPattern(pattern = consensus.seq, subject = len_filter) asks whether a
	# read contains a particular sequence. If it contains it (which is what we
	# want) the function gives us a 1 - so we select for all "1"s
	cons_seq_filter <- len_filter[which(vcountPattern(pattern = consensus.seq, subject = len_filter) == 1)]
	
	# Almost identical to above, but here if a read contains an "N", we get a "1",
	# so we select for all "0"s
	N_filter <- cons_seq_filter[which(vcountPattern(pattern = DNAString(x = "N"), subject = cons_seq_filter) == 0)]
	
	# Finally, we want to see how our filtering affects the number of reads at
	# each step. length(reads) gives us the number of reads. The data.frame()
	# function is a way of creating a table - a data frame in R. I define each
	# column name, and give each columns a value - the number of reads. The
	# stringsAsFactors statement at the end tells R to just treat the raw values
	# as raw values. The details are not too important.
	filter_df <- data.frame(RawReads = p_set$num_raw_reads, LenFiltReads = length(len_filter), ConsFiltReads = length(cons_seq_filter), NFiltReads = length(N_filter), stringsAsFactors = FALSE)
	
	# As above - for each P-set return the final filtered read set and a table of
	# the read numbers at each filter step
	return(list(filt_seq_data = N_filter, filter_nums = filter_df))
})

# We can confirm that no Ns remain in these sequence sets by running e.g.:
alphabetFrequency(filtSeq_data_list$P1$filt_seq_data, collapse = TRUE)

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
uniqNums_list <- lapply(1:length(uniqSeq_data_list), function(p_set_index) {
	p_set_name <- names(uniqSeq_data_list)[p_set_index]
	num_uniq_reads <- uniqSeq_data_list[[p_set_index]]$UniqRead
	data_with_name <- data.frame(PSet = p_set_name, UniqReads = num_uniq_reads, stringsAsFactors = FALSE)
	return(data_with_name)
})
uniqNums_df <- bind_rows(uniqNums_list)

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

#### Nucleotide-level analysis
### GC-content differences between P-sets
PSet_GC_content_list <- lapply(PSet_name_list, function(p_set_name) {
	# Calculate the frequencies of G and C nucleotides seperately for each sequence
	PSet_GC_sep_freq	<- alphabetFrequency(finalUniqSeq_data_list[[p_set_name]]$UniqSeqData, as.prob = TRUE)[,c("C", "G")]
	# Combine these to get the GC combined frequency
	PSet_GC_tot_freq	<- rowSums(PSet_GC_sep_freq)
	# Return as dataframe with two columns: P-set name and GC frequency
	PSet_GC_tot_df		<- data.frame(PSet = rep(p_set_name, length(PSet_GC_tot_freq)), GC.Content = PSet_GC_tot_freq)
	return(PSet_GC_tot_df)
})
# Bind all the individual P-set results into a single dataframe
PSet_GC_content_df	<- bind_rows(PSet_GC_content_list)

# A list of pairwise comparisons we'll want to do for stats 
PairComp_list		<- list(c("P1", "P4"), c("P4", "P8"), c("P1", "P8"))

# Boxplot distributions of GC content for each P-set
PSet_GC_content_plot	<- ggplot(data = PSet_GC_content_df, aes(x = PSet, y = GC.Content, color = PSet)) +
	geom_boxplot() +
	theme(panel.background = element_blank(), axis.line = element_line(color = "black")) +
	stat_compare_means(comparisons = PairComp_list, method = "wilcox.test", p.adjust = "bonferroni")


### Codon usage bias

# When translating and assinging codons, assume the bacterial code
bacterial_code	<- 11
BacGeneCode		<- getGeneticCode(as.character(bacterial_code))

## We calculate codon usage using the Codon Adaptability Index for E. coli. The
## codons used in each (unique) sequence are compared to an index (the CAI) of
## optimal codons derived from highly expressed genes in E. coli (I think). We
## compare how the distributions of codons in each sequence in each P-set
## differs to the others.

# Read in Genscript CAI data
GS_cai_tbl	<- read.table(file.path(fastq_file_path, "GenScript_ecoli_caiTbl.txt"), sep = "\t", header = TRUE, colClasses = "character")
GS_cai_tbl	<- GS_cai_tbl[order(GS_cai_tbl$Triplet),]
gs_cai_indx <- as.numeric(GS_cai_tbl$Fraction)

# The default function for calculating CAI is slow. We can make it faster by
# calculating some constant variables outside the function (i.e. we don't need
# to recalculate the same value for each sequence processed)
cai.faster  <- function (seq, w, exclude) {
	nncod <- uco(seq)
	nncod <- nncod[-exclude]
	sigma <- nncod %*% log(w)
	exp(sigma/sum(nncod))
}
# The vectorize function allows us to provide a vector or list as an argument,
# rather than looping over each sequence
cai.faster_vf <- Vectorize(cai.faster, vectorize.args = "seq", SIMPLIFY = TRUE)

# Pre-run some cai calculations that were part of the old cai() function
stops		<- which("Stp" == aaa(translate(s2c(c2s(words())), numcode = bacterial_code)))
singulets   <- which(sapply(syncodons(words(), numcode = bacterial_code), length) == 1)
exclude		<- c(stops, singulets)
gs_cai_adj	<- gs_cai_indx[-exclude]
gs_cai_adj[gs_cai_adj < 1e-04] <- 0.01

# For each sequence in each P-set, calculate the CAI
message("WARNING: CAI calculation takes a while - about 1-2 minutes per P-set")
PSetCAI_genscr_list <- lapply(PSet_name_list, function(p_set_name) {
	message(paste0("Calculating CAI for P-set: ", p_set_name))
	
	# The cai() function takes character vectors split by each individual base (lower case)
	charSplit_seqs <- str_split(string = tolower(as.character(finalUniqSeq_data_list[[p_set_name]]$UniqSeqData)), pattern = "")
	cai_result <- cai.faster_vf(charSplit_seqs, w = gs_cai_adj, exclude = exclude)
	
	# Convert result to a dataframe structure
	cai_df <- data.frame(PSet = rep(p_set_name, length(cai_result)), CAI.Value = cai_result)
	return(cai_df)
})

# Bind all the individual P-set results into a single dataframe
PSetCAI_genscr_df	<- bind_rows(PSetCAI_genscr_list)
# Plot a boxplot of the CAI scores distributions and test difference
PSetCAI_genscr_plot	<- ggplot(data = PSetCAI_genscr_df, aes(x = PSet, y = CAI.Value, color = PSet)) +
	geom_boxplot() +
	theme(panel.background = element_blank(), axis.line = element_line(color = "black")) +
	stat_compare_means(comparisons = PairComp_list, method = "t.test", p.adjust = "bonferroni")



### Look at positions of STOP codons within the sequences
# Define the stop codon
stop_AA		<- "*"

# For all the unique P-set sequences: translate, check for stop codons, get their locations
stop_AllPosits_list <- lapply(PSet_name_list, function(p_set_name) {
	# Get the (unique) nucloetide sequences for each P-set and r
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
# Normalise the numbers to 1 - where 1 is the highest scoring P-set in each category
normal_to_1			<- function(x){((x)/(max(x)))}
stop_AllNums_norm	<- data.frame(PSet = stop_AllNums_df[,1], apply(stop_AllNums_df[,-1], 2, normal_to_1))
# Melt dataframe for ggplot
stop_AllNums_melt	<- melt(stop_AllNums_norm, id.vars = "PSet", variable.name = "Count.Type", value.name = "Counts")

# Barplot to compare how many stop codons are present in the Pset as a
# proportion of all codons and how many sequences contains stop codons as a
# proportion of all sequences. In both cases these fractions are normalised to
# the largest Pset for each type of count
stop_AllNums_plot	<- ggplot(data = stop_AllNums_melt, aes(x = PSet, y = Counts, fill = Count.Type)) +
	geom_bar(stat = "identity", position = "dodge") +
	labs(y = "Normalised fraction ofall seqs in PSet") +
	scale_fill_discrete(
		name="Type of count",
		breaks=c("Raw.Num", "InSeq.Num"),
		labels=c("Total Stop codons", "Num. seqs with Stops")) +
	theme(legend.title.align = 0.5)


# Bind all the individual P-set STOP position results into a single dataframe
stop_AllPosits_df	<- bind_rows(lapply(stop_AllPosits_list, function(pset) return(pset$stop_posits_df)))
# Boxplot of the STOP position distribution across P-sets with pairwise significance tests
stop_AllPositis_plot	<- ggplot(data = stop_AllPosits_df, aes(x = PSet, y = StopPosition, color = PSet)) +
	geom_boxplot() +
	ylim(0, 140) +
	geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 51, ymax = 58), color = "black", linetype = "solid", size = 0.5, fill = NA) +
	geom_hline(yintercept = 108, color = "red", linetype = "dashed") +
	ggtitle("Positions of stop codons in sequences containing stop codons") +
	theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "grey40")) +
	stat_compare_means(comparisons = PairComp_list, method = "wilcox.test", p.adjust = "bonferroni")




### Look at positions of Shine-Dalgarno sequence in P-sets

# Define the strictest SD consensus sequence
SD_seq		<- DNAString("AGGAGGT")

# For each P-set calculate how many SD sequences are present, how many sequences contain SD sites, and get a list of these SD sequences
SD_AllPosits_list <- lapply(PSet_name_list, function(p_set_name) {
	# Get the (unique) nucloetide sequences for each P-set and r
	nucleot_seqs	<- finalUniqSeq_data_list[[p_set_name]]$UniqSeqData
	
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
names(SD_AllPosits_list)	<- PSet_name_list
SD_AllNums_df	<- bind_rows(lapply(SD_AllPosits_list, function(pset) return(pset$SD_nums_df)))

## Combine SD position list into single dataframe 
SD_AllPosits_df	<- bind_rows(lapply(SD_AllPosits_list, function(pset) return(pset$SD_posits_df)))

# Boxplot of the SD position distribution across P-sets with pairwise significance tests
SD_AllPos_boxplot	<- ggplot(data = SD_AllPosits_df, aes(x = PSet, y = SDPosition, color = PSet)) +
	geom_boxplot() +
	ylim(0, 450) +
	geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 151, ymax = 174), color = "black", linetype = "solid", size = 0.5, fill = NA) +
	geom_hline(yintercept = 324, color = "red", linetype = "dashed") +
	ggtitle("Positions of SD seqs in sequences containing SD seqs") +
	theme(panel.background = element_blank(), plot.title = element_text(hjust = 0.5), axis.line = element_line(colour = "grey40")) +
	stat_compare_means(comparisons = PairComp_list, method = "wilcox.test", p.adjust = "bonferroni")

# Density plot of SD positions in all P-sets
SD_AllPos_densplot	<- ggplot(data = SD_AllPosits_df, aes(SDPosition, color = PSet)) +
	geom_density()

# Break the 324-bp sequences containing SD into 9 x 36 bp chunks so we can plot seqlogos
SD_SeqBreaks_list <- lapply(PSet_name_list, function(p_set_name) {
	nucleot_seqs	<- SD_AllPosits_list[[p_set_name]]$SD_contain_seqs
	nucleo_chars	<- as.character(nucleot_seqs)
	
	splitVector		<- split(1:324, ceiling(seq_along(1:324)/36))
	splitNucleoCassette <- lapply(splitVector, function(cassette) {
		cassette_range	<- c(min(cassette), max(cassette))
		sequence_split 	<- lapply(nucleo_chars, function(char_seq) {
			str_sub(char_seq, cassette_range[1], cassette_range[2])
		})
		return(unlist(sequence_split))
	})
})
names(SD_SeqBreaks_list) <- PSet_name_list

# Only worth plotting seqLogos for P-values and chunks we think might be of interest - e.g. at the density spike in P1
# In P1 - the spike is between 100-150 bps from start, so we plot bins 3,4,5 to cover that area
P1_SD_spike_logo <- ggseqlogo(data = SD_SeqBreaks_list$P1[3:5], method = "probability", ncol = 1, nrow = 3)




