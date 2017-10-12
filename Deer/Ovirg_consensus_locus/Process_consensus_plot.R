# / Libraries / #
library(scales)
library(ggplot2)
library(seqinr)
library(stringr)

# / Function that produces a color palette to mimic ggplot / #
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}

# / Directory containing input data / #
setwd("/Users/aesin/Desktop/Deer/Assembly/Ovirg_scaffold_0917/")

# / Read in Data / #
HB_locus_align	<- read.alignment(file = "Ovtex_Ovirg_alignment.fasta", format = "fasta")
Annot_data_raw	<- read.csv(file = "Ovtex_Ovirg_alignment_annotation.csv", header = TRUE)
Ident_data		<- read.csv(file = "Ovtex_Ovirg_alignment_ID.csv", header = TRUE)

# / Ovirg consensus is the first of the two sequences in the alignment / #
# / Extract all the bases / #
Ovirg_seq_l		<- HB_locus_align$seq[[1]]
Ovirg_seq_raw	<- unlist(str_split(Ovirg_seq_l, ""))

# / Make a dataframe of sequences, positions, and per-base identity / #
raw_ident_df <- data.frame(
	Base		= Ovirg_seq_raw,
	Position	= Ident_data$Position,
	Ident_score = Ident_data$Identity,
	stringsAsFactors = FALSE
)


# / Different sens cutoffs give different number of "bad" windows.
# The lower the cutoff, the more bad positions, and the more
# stringent the inclusion in the good windows / #

# Cutoff:	25 = 528 bad windows
#			50 = 368 bad windows
#			75 = 273 bad windows

# / Define a window and an a cutoff. If a window contains more
# gaps than the sensitivity cutoff, then the window is given a 
# negative value and will be excluded. If it's included, then
# calculate mean identity over all the standard bases. / #
mean_ident_100_list	 <- by(raw_ident_df, ceiling(1:length(raw_ident_df$Position) / 100), function(window) {

	# / Define cutoff / #
	sens_cutoff	<- 50
	# Count gaps and poor-quality bases
	count_gaps	<- nrow(window[(window$Base %in% c("?","-")),])
	count_N		<- nrow(window[window$Base == "n",])

	# / Assign different -ve values for potential sorting / #
	if (count_gaps > sens_cutoff) {
		return(-0.1)
	} else if (count_N > sens_cutoff) {
		return(-0.2)
	} else {
		# / Get identity scores for all standard bases / #
		idents	<- window$Ident_score[!(window$Base %in% c("?","n", "-"))]
		# / If no. good bases > size minus sens cutoff, include window / #
		if (length(idents) >= (100 - sens_cutoff)) {
			mean_ident	<- (1 - mean(idents))
			return(mean_ident)
		}
		# / Return a -ve if num. gaps + bad bases > sens cutoff / #
		return(-0.15)		
	}
}, simplify = FALSE)

# / Data frame to hold the 100bp windows and their mean identity scores / #
av_ident_df	<- data.frame(
	Position		= seq(1:length(mean_ident_100_list)) * 100,
	Average.Ident	= as.numeric(mean_ident_100_list)
)

# / Process the annotation data - for each CDS region, we round down
# for the minimum (upstream) and up for maximum (upstream) to get
# slightly generous cassettes. / #
Annot_list		<- lapply(1:nrow(Annot_data_raw), function(index) {
	row		<- Annot_data[index,]
	minimum	<- 100 * floor(row[2] / 100)
	maximim	<- 100 * ceiling(row[3] / 100)
	return(data.frame(row[1], minimum, maximim))
})
Annot_data_df <- bind_rows(Annot_list)

# / A vector of all the 100bp windows that contains genes / #
gene_positions <- unlist(apply(Annot_data_df, 1 ,function(gene) {
	seq(as.numeric(gene[2]), as.numeric(gene[3]), by = 100)
}))


# / First create another col to hold data for missing (gaps, poor qual)
# data. Any val < 0 in ident becomes -0.025 in second col. All good
# windows are NAs in Missing.Value col. Then convert all bad values
# in original Average.Ident col to NA - this way we can easily plot
# the two datasets seperately. /#
av_ident_df$Missing.Value	<- apply(av_ident_df, 1, function(row) ifelse(row[2] < 0, -0.025, NA))
av_ident_df$Average.Ident	<- apply(av_ident_df, 1, function(row) ifelse(row[2] < 0, NA, row[2]))


# / Same principle as above for the gene positions. If a window is part
# of a gene, add a value to CDS.Region col, otherwise NA. Inverse for
# intergene col, but value is the same because we want to plot them on the
# same line. / #
av_ident_df$CDS.Region	<- apply(av_ident_df, 1, function(row) ifelse(row[1] %in% gene_positions, -0.05, NA))
av_ident_df$Intergene	<- apply(av_ident_df, 1, function(row) ifelse(row[1] %in% gene_positions, NA, -0.05))


# / Subset data by position, to plot short ranges / #
av_ident_sub	<- av_ident_df[which(av_ident_df$Position > 100000),]

# / We plot in different steps, so set up color palette here / #
palette = gg_color_hue(2)

# / Run plot / #
ident_plot <- ggplot(data = av_ident_sub, aes(x = Position, y = Average.Ident)) +
	# / Points and a smooth loess line of nucleotide discordance / #
	geom_point(size = 0.2, color = palette[1]) +
    stat_smooth(geom = "smooth", method = "loess", span = 0.1, size = 1.25, se = F, color = palette[1]) +

    # / Lines connect adjacent windows with poor resolution (gaps/poor mapping). 
    # Solitary poor windows can be visualised with the geo_point command below / #
    geom_line(aes(x = Position, y = Missing.Value), size = 10, color = "black") +
    # geom_point(aes(x = Position, y = Missing.Value), size = 0.2, color = "black") +

    # / A standard cassette representation of gene positions / #
    geom_line(aes(x = Position, y = CDS.Region), size = 10, color = palette[2]) +
    geom_line(aes(x = Position, y = Intergene), size = 0.5, color = palette[2]) +

    ylab("Mean discordance over 100bp window") +
    # / Put commas in the x-axis position labels / #
    scale_x_continuous(labels = comma)



# / T-test whether HBB discordance different to overall / #
all_valid_pos	<- av_ident_df[av_ident_df$Missing == FALSE,]
hbb_posits		<- subset(all_valid_pos, Position > 325499 & Position < 328501)
t.test(hbb_posits$Average.Ident, all_valid_pos$Average.Ident)
# p = 0.0006576

# / Expand HBB cassette (2kb +/-) / #
hbb_posits		<- subset(all_valid_pos, Position > 323499 & Position < 330501)
t.test(hbb_posits$Average.Ident, all_valid_pos$Average.Ident)
# p = 2.2e-08














