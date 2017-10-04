library(seqinr)
library(stringr)

HB_locus_align	<- read.alignment(file = "Ovtex_Ovirg_alignment.fasta", format = "fasta")
Ident_data		<- read.csv(file = "Ovtex_Ovirg_alignment_ID.csv", header = T)

# In the alignment file the consensus sequence for Ovirg is first
Ovirg_seq_l		<- HB_locus_align$seq[[1]]
Ovirg_seq_raw	<- unlist(str_split(Ovirg_seq_l, ""))

raw_ident_df <- data.frame(Base = y, Position = Ident_data$Position, Ident_score = Ident_data$Identity, stringsAsFactors = F)
# raw_ident_sub <- raw_ident_df[which(raw_ident_df$Position > 190201 & raw_ident_df$Position < 190401),]

mean_ident_100_list	 <- by(raw_ident_df, ceiling(1:length(raw_ident_df$Position) / 100), function(x) {
	sens_cutoff	<- 50
	# print(x)
	count_gaps	<- nrow(x[(x$Base %in% c("?","-")),])
	count_N		<- nrow(x[x$Base == "n",])
	if (count_gaps > sens_cutoff) {
		return(-0.1)
	} else if (count_N > sens_cutoff) {
		return(-0.2)
	} else {
		idents	<- x$Ident_score[!(x$Base %in% c("?","n", "-"))]
		if (length(idents) >= (100 - sens_cutoff)) {
			mean_ident	<- (1 - mean(idents))
			return(mean_ident)
		}
		return(-0.15)		
	}
}, simplify = F)

av_ident_df	<- data.frame(Position = seq(1:length(mean_ident_100_list))*100, Average.Ident = as.numeric(mean_ident_100_list))

## // Optional // ##
av_ident_df$Average.Ident	<- apply(av_ident_df, 1, function(row) {
	# message(row)
	if (row[2] < 0) {
		return(-0.1)
	} else {
		return(row[2])
	}
})

av_ident_df$Missing	<- apply(av_ident_df, 1, function(row) {
	# message(row)
	if (row[2] < 0) {
		TRUE
	} else {
		FALSE
	}
})

# with sens 25 missing = 528 positions
# with sens 50 missing = 368 positions
# with sens 75 missing = 273 positions

av_ident_sub	<- av_ident_df[which(av_ident_df$Position > 100000),]
# av_ident_sub	<- av_ident_sub[which(av_ident_sub$Position < 371000),]

ident_plot <- ggplot(data = av_ident_sub, aes(x = Position, y = Average.Ident, color = Missing)) +
    stat_smooth(data = subset(av_ident_sub, Average.Ident > 0), aes(x = Position, y = Average.Ident, color = Missing), geom = "smooth", method = "loess", span = 0.1, size = 0.5, se = F) +
    geom_point(size = 0.25) +
    geom_vline(xintercept = 298000, color = "green") +
    geom_vline(xintercept = 325500, color = "black") +
    geom_vline(xintercept = 328500, color = "black") +
    # geom_vline(xintercept = 281200, color = "black") +
    geom_vline(xintercept = 371000, color = "green") +
    ylab("Mean discordance over 100bp window")


all_valid_pos	<- av_ident_df[av_ident_df$Missing == FALSE,]
# Actual HBB
hbb_posits		<- subset(all_valid_pos, Position > 325499 & Position < 328501)
t.test(hbb_posits$Average.Ident, all_valid_pos$Average.Ident)

# / p = 0.0006576

# Slighty larger cassette HBB (2kb +/-)
hbb_posits		<- subset(all_valid_pos, Position > 323499 & Position < 330501)
t.test(hbb_posits$Average.Ident, all_valid_pos$Average.Ident)

# / p = 2.2e-08


