library(ggplot2)
library(reshape2)
library(plyr)
library(RSvgDevice)

MyMerge <- function(x, y){
	df <- merge(x, y, by = "row.names", all = TRUE, sort = FALSE)
	rownames(df)  <- df$Row.names
	df$Row.names  <- NULL
	return(df)
}

plot_bar_fraction <- function (data_table, x, y, group, ylimit = 100) {
	ggplot(data = data_table, aes_string(x = x, y = y, group = group, fill = x)) +
		geom_bar(position = "dodge", stat = "identity") +
		# geom_text(aes(label = round(Fraction.of.Groups, digits = 2), vjust = -1)) +
		facet_wrap(group, nrow = 4) +
		theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
		scale_y_continuous("Fraction of Groups", limits = c(0, ylimit)) +
		theme(axis.text = element_text(size = 12, colour = "black"), axis.title.x = element_text(size = 16, colour = "black", vjust = -0.8), axis.title.y = element_text(size = 16, colour = "black", vjust = 1), plot.title = element_text(size = 18, colour = "black", vjust = 0.5))
}

###########################################################################
direct <- "/users/aesin/desktop/Mowgli/External_HGT/New_parser/Donor_edges"

min_taxa = 0
penalty_list <- 4:6
rank <- "Class"

###########################################################################
## Make directories ##
setwd(direct)
dir.create("Plots", showWarnings = FALSE)

###########################################################################
## For each penalty, read in the donor ranks ##

setwd(paste0(direct, "/", rank))

ordered_ranks_list <- vector("list", length(penalty_list))
list_counter = 1

for (penalty in penalty_list) {
	## Read in and order table ##
	ranks_tab <- table(read.table(file = paste0("T", penalty, "_donor_", tolower(rank), "_min_", min_taxa, ".tsv"), header = FALSE, sep = "\t")[,2])
	ranks_tab_ordered <- as.data.frame(ranks_tab[order(-ranks_tab)])

	## Write out the tables individually ##
	name_donor_table <- paste0(penalty, "_donors")
	assign(name_donor_table, ranks_tab_ordered)

	## Place into the list ##
	ordered_ranks_list[[list_counter]] <- ranks_tab_ordered
	colnames(ordered_ranks_list[[list_counter]])[1] <- paste0("Donors_T", penalty)
	list_counter = list_counter + 1
}

## Merge the individual donor tables into one and convert any NAs into 0s ##
merged_donors <- Reduce(MyMerge, ordered_ranks_list)
merged_donors[is.na(merged_donors)] <- 0

## Make the first column the Phylum names ##
new_names <- gsub(" group", "", rownames(merged_donors))
merged_donors <- cbind(new_names, merged_donors)
rownames(merged_donors) <- NULL
colnames(merged_donors)[1] <- rank


## Make the data a fraction of the total donors ##
merged_donors_fraction <- cbind(merged_donors[,1], (sweep(merged_donors[,-1], 2, colSums(merged_donors[,-1]), "/") * 100))
colnames(merged_donors_fraction)[1] <- rank

merged_fraction_no_top <- merged_donors_fraction[-1,]
top_name <- as.character(merged_donors_fraction[1,1])

## Melt the tables and also order the levels by the greatest fraction first ##
merged_donors_fraction.molten <- melt(merged_donors_fraction, value.name="Fraction.of.Groups", variable.name="Penalty", na.rm = TRUE)
merged_fraction_no_top.molten <- melt(merged_fraction_no_top, value.name="Fraction.of.Groups", variable.name="Penalty", na.rm = TRUE)

merged_donors_fraction.molten.ordered <- merged_donors_fraction.molten
merged_fraction_no_top.molten.ordered <- merged_fraction_no_top.molten

merged_donors_fraction.molten.ordered[[1]] <- factor(merged_donors_fraction.molten[[1]], levels = unique(merged_donors_fraction.molten[[1]][order(-merged_donors_fraction.molten$Fraction.of.Groups)]))
merged_fraction_no_top.molten.ordered[[1]] <- factor(merged_fraction_no_top.molten[[1]], levels = unique(merged_fraction_no_top.molten[[1]][order(-merged_fraction_no_top.molten$Fraction.of.Groups)]))


## Apply the plot functions and label the plots ##
plot_ranks_fraction <- plot_bar_fraction(merged_donors_fraction.molten.ordered, "Penalty", "Fraction.of.Groups", rank)
#devSVG(file = "/users/aesin/Dropbox/LSR/Presentation/Donor_groups/Class_donors.svg", width = 20, height = 20)
plot_ranks_fraction + ggtitle(paste0(rank, " donors into Anoxy/Geobacillus with minimum sequence cutoff of ", min_taxa))
#dev.off()

#ggsave(path = paste0(direct, "/Plots"), filename = paste0("Donor_", rank, "_across_penalties_min_context_", min_taxa, ".png"), plot = plot_ranks_fraction, width = 24, height = 16)
#ggsave(path = paste0(direct, "/Plots"), filename = paste0("Donor_", rank, "_across_penalties_min_context_", min_taxa, ".pdf"), plot = plot_ranks_fraction, width = 24, height = 16)

# no_top_y_scale <- round_any(max(merged_fraction_no_top.molten.ordered$Fraction.of.Groups), f = ceiling, 10)
# plot_ranks_fraction_no_top <- plot_bar_fraction(merged_fraction_no_top.molten.ordered, "Penalty", "Fraction.of.Groups", rank, ylimit = no_top_y_scale)
# plot_ranks_fraction_no_top <- plot_ranks_fraction_no_top + ggtitle(paste0(rank, " donors into Anoxy/Geobacillus with minimum sequence cutoff of ", min_taxa, " - ", top_name, " removed"))


# ggsave(path = paste0(direct, "/Plots"), filename = paste0("Donor_", rank, "_across_penalties_no_", top_name, "_min_context_", min_taxa, ".png"), plot = plot_ranks_fraction_no_top, width = 24, height = 16)
# ggsave(path = paste0(direct, "/Plots"), filename = paste0("Donor_", rank, "_across_penalties_no_", top_name, "_min_context_", min_taxa, ".pdf"), plot = plot_ranks_fraction_no_top, width = 24, height = 16)

# ## Get the total donors for each penalty ##
# totals <- colSums(merged_donors[,-1])


