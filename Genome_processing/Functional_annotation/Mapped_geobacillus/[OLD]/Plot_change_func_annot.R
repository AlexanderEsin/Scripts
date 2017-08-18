library(ggplot2)
library(reshape2)
library(gridExtra)

MyMerge <- function(x, y){
	df <- merge(x, y, by = "row.names", all = TRUE, sort = FALSE)
	rownames(df)  <- df$Row.names
	df$Row.names  <- NULL
	return(df)
}

###########################################################################
IG_included = TRUE
All_scenarios = FALSE

## Define subset of HGT predictions to be used ##
penalty_folder <- "t3_t4_t5"

###########################################################################

if (All_scenarios == TRUE) {
	Scenario_ID <- "All_scenarios"
	merged_relative_HGT_title <- "All_scenarios"
} else {
	Scenario_ID <- "Scenarios_1_2"
	merged_relative_HGT_title <- "Scenarios 1 & 2 ONLY"
}

if (IG_included == TRUE) {
	IG_ID <- "IG"
} else {
	IG_ID <- "NO_IG"
}

IG_folder <- paste0("/", IG_ID, "/", Scenario_ID)

###########################################################################
## Prepare the list of penalties to be queried downstream ##
penalty_list <- strsplit(penalty_folder, "_")[[1]]
penalty_num <- length(penalty_list)

## Prepare list for the ordered HGT tables, but also write them out to individual variables ##
ordered_hgt_table_list <- vector("list", penalty_num)
ordered_vert_table_list <- vector("list", penalty_num)
list_entry_counter <- 1

###########################################################################
## Read in total and merged tables ##

## The functional annotation for the complete set of gene families/groups/clusters ##
setwd("/users/aesin/desktop/Geo_analysis/Geo_ortholog_nucl/Functional_annotation")
total_cog_tab <- table(read.table(file = "Narrow_COG_list.txt", sep = "\n"))
total_cog_tab_ordered <- as.data.frame(total_cog_tab[order(-total_cog_tab)])
colnames(total_cog_tab_ordered)[1] <- "Total"
# Add to our table lists #
ordered_hgt_table_list[[list_entry_counter]] <- total_cog_tab_ordered
ordered_vert_table_list[[list_entry_counter]] <- total_cog_tab_ordered
list_entry_counter = list_entry_counter + 1

## The functional annotation for the sets predicted as HGT or Vertical by ALL the penalty reconciliation, e.g.: 1 1 1 1 or 0 0 0 0 ##
setwd(paste0("/users/aesin/desktop/Geo_analysis/HGT_sets/Functional_annotation/Merged", IG_folder, "/", penalty_folder))

hgt_cog_tab <- table(read.table(file = "HGT_narrow_COG.tsv", sep = "\n"))
hgt_cog_tab_ordered <- as.data.frame(hgt_cog_tab[order(-hgt_cog_tab)])
colnames(hgt_cog_tab_ordered)[1] <- "All_HGT"
# Add to our table list #
ordered_hgt_table_list[[list_entry_counter]] <- hgt_cog_tab_ordered
ordered_vert_table_list[[list_entry_counter]] <- hgt_cog_tab_ordered
list_entry_counter = list_entry_counter + 1

vertical_cog_tab <- table(read.table(file = "Vertical_narrow_COG.tsv", sep = "\n"))
vertical_cog_tab_ordered <- as.data.frame(vertical_cog_tab[order(-vertical_cog_tab)])
colnames(vertical_cog_tab_ordered)[1] <- "All_Vertical"
# Add to our table list #
ordered_hgt_table_list[[list_entry_counter]] <- vertical_cog_tab_ordered
ordered_vert_table_list[[list_entry_counter]] <- vertical_cog_tab_ordered
list_entry_counter = list_entry_counter + 1

###########################################################################
## Read in the table corresponding to the refined dataset containing transfers from only OUTSIDE the IG ##

setwd("/users/aesin/desktop/Mowgli/External_HGT/New_parser/Functional_annotation")
maximum_penalty <- max(as.numeric(gsub(pattern = "[a-z]+", replacement = "", x = penalty_list)))
refined_cog_tab <- table(read.table(file = paste0("T", maximum_penalty, "_narrow_COG.tsv"), sep = "\n"))
refined_cog_tab_ordered <- as.data.frame(refined_cog_tab[order(-refined_cog_tab)])
colnames(refined_cog_tab_ordered)[1] <- "Refined_external_only"
## Cloned the HGT table list ##
ordered_hgt_refined_tbl_list <- ordered_hgt_table_list
ordered_hgt_refined_tbl_list[[list_entry_counter]] <- refined_cog_tab_ordered
ordered_hgt_refined_tbl_list <- Filter(length, ordered_hgt_refined_tbl_list)

###########################################################################
## Read in the individual HGT and vertical tables ##
for (penalty in penalty_list) {
	penalty <- toupper(penalty)
	setwd(paste0("/users/aesin/desktop/Geo_analysis/HGT_sets/Functional_annotation/", penalty, IG_folder))

	## Read in and order tables ##
	read_in_hgt_table <- table(read.table(file = "HGT_narrow_COG.tsv", sep = "\n"))
	ordered_hgt_table <- as.data.frame(read_in_hgt_table[order(-read_in_hgt_table)])
	read_in_vert_table <- table(read.table(file = "Vertical_narrow_COG.tsv", sep = "\n"))
	ordered_vert_table <- as.data.frame(read_in_vert_table[order(-read_in_vert_table)])

	## Write out the tables individually ##
	name_hgt <- paste0(penalty, "_hgt_tab_ordered")
	name_vert <- paste0(penalty, "_vert_tab_ordered")
	assign(name_hgt, ordered_hgt_table)
	assign(name_vert, ordered_vert_table)

	## Put the table into the list, and rename the columns accordingly for merging ##
	ordered_hgt_table_list[[list_entry_counter]] <- ordered_hgt_table
	colnames(ordered_hgt_table_list[[list_entry_counter]])[1] <- paste0("HGT_", penalty)
	ordered_vert_table_list[[list_entry_counter]] <- ordered_vert_table
	colnames(ordered_vert_table_list[[list_entry_counter]])[1] <- paste0("Vert_", penalty)

	list_entry_counter = list_entry_counter + 1
}

###########################################################################
## Merge the values into a single table ##
merged_hgt <- Reduce(MyMerge, ordered_hgt_table_list)
merged_vert <- Reduce(MyMerge, ordered_vert_table_list)

# Convert any NAs to 0 #
merged_hgt[is.na(merged_hgt)] <- 0
merged_vert[is.na(merged_vert)] <- 0

## Output has the COG labels as rownames -> remove rownames and add a column corresponding to the COG functional categories ##
COG_classes <- rownames(merged_hgt)
COG_classes <- rownames(merged_vert)

rownames(merged_hgt) <- NULL
rownames(merged_vert) <- NULL

merged_hgt <- cbind(COG_classes, merged_hgt)
merged_vert <- cbind(COG_classes, merged_vert)

## Rearrange the columns such that All_HGT & All_Vertical appear at the end (independent of column / penalty number) ##
merged_hgt <- cbind(merged_hgt[,-3:-4], merged_hgt[,3:4])
merged_vert <- cbind(merged_vert[,-3:-4], merged_vert[,3:4])

###########################################################################

## Pull out the penalty-specific column names, coerce into data frame ##
penalty_cols_hgt <- colnames(merged_hgt[ ,3:(ncol(merged_hgt)-2)])
penalty_cols_vert <- colnames(merged_vert[ ,3:(ncol(merged_vert)-2)])

HGT_sums <- colSums(merged_hgt[penalty_cols_hgt])
Vert_sums <- colSums(merged_vert[penalty_cols_vert])

hgt_penalty_nums <- gsub(pattern = "HGT_T","", penalty_cols_hgt)
vert_penalty_nums <- gsub(pattern = "Vert_T","", penalty_cols_vert)

HGT_num_df <- data.frame("Penalty" = hgt_penalty_nums, "Groups_w_HGT" = HGT_sums, row.names = NULL)
Vert_num_df <- data.frame("Penalty" = vert_penalty_nums, "Groups_w_Vert" = Vert_sums, row.names = NULL)

###########################################################################
## Make the plots ##
as_fraction_total_hgt <- sweep(merged_hgt[,-1],2,colSums(merged_hgt[,-1]), "/")*100
as_fraction_total_vert <- sweep(merged_vert[,-1],2,colSums(merged_vert[,-1]), "/")*100

final_table_hgt <- cbind(COG = merged_hgt[,1], as_fraction_total_hgt)
final_table_vert <- cbind(COG = merged_vert[,1], as_fraction_total_vert)

final_hgt.molten <- melt(final_table_hgt, value.name="Percentage.of.group", variable.name="COG.category", na.rm = TRUE)
final_vert.molten <- melt(final_table_vert, value.name="Percentage.of.group", variable.name="COG.category", na.rm = TRUE)

final_hgt.molten.ordered <- final_hgt.molten
final_vert.molten.ordered <- final_vert.molten

final_hgt.molten.ordered$COG <- factor(final_hgt.molten$COG, levels = unique(final_hgt.molten$COG[order(-final_hgt.molten$Percentage.of.group)]))
final_vert.molten.ordered$COG <- factor(final_vert.molten$COG, levels = unique(final_vert.molten$COG[order(-final_vert.molten$Percentage.of.group)]))

#plot_1 <- ggplot(data = final.molten.ordered, aes(x = COG.category, y = Percentage.of.group, group = COG)) + geom_bar(stat = "identity") + facet_wrap("COG", nrow = 1) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

plot_hgt <- ggplot(data = final_hgt.molten.ordered, aes(x = COG.category, y = Percentage.of.group, group = COG)) + 
			geom_bar(position = "dodge", stat = "identity", colour = "black", fill = "white") + 
			facet_wrap("COG", nrow = 1) + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			ggtitle(paste0("Proportion of putative HGT events per COG over a penalty range - ", merged_relative_HGT_title)) + 
			geom_line(data = subset(final_hgt.molten.ordered, COG.category %in% penalty_cols_hgt), colour = "red", size = 1)

plot_vert<- ggplot(data = final_vert.molten.ordered, aes(x = COG.category, y = Percentage.of.group, group = COG)) + 
			geom_bar(position = "dodge", stat = "identity", colour = "black", fill = "white") + 
			facet_wrap("COG", nrow = 1) + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			ggtitle(paste0("Proportion of putative vertical events per COG over a penalty range - ", merged_relative_HGT_title)) + 
			geom_line(data = subset(final_vert.molten.ordered, COG.category %in% penalty_cols_vert), colour = "red", size = 1)

###########################################################################
## Get the table which shows how many groups were original reconciled to produce the data output ## 
table_name <- paste0(IG_ID, "_intersect/", Scenario_ID, "/", "Groups_reconciled_per_penalty_", penalty_folder, ".tsv")
penalty_group_num_df <- read.table(paste0("/users/aesin/desktop/Mowgli/Mowgli_outputs/", table_name), header = TRUE)
names(penalty_group_num_df)[1] <- "Penalty"


merged_penalty_df <- merge(penalty_group_num_df, HGT_num_df, by = "Penalty")
merged_penalty_df <- merge(merged_penalty_df, Vert_num_df, by = "Penalty")
merged_penalty_df.molten <- melt(merged_penalty_df, id.vars = "Penalty")

plot_vert_hgt_total_num <- ggplot(merged_penalty_df.molten, aes(x = Penalty, y = value, colour = variable)) + geom_line() + geom_point()

###########################################################################
## Make a plot to show the relative (fold) enrichment of HGT events per functional category ##
final_table_hgt_2 <- final_table_hgt
final_table_hgt_2$relative_HGT <- (final_table_hgt_2$All_HGT / final_table_hgt_2$All_Vertical)
relative_hgt_ext <- final_table_hgt_2[1:19, c(1, ncol(final_table_hgt_2))]
relative_hgt_ext$COG <- factor(relative_hgt_ext$COG, levels = unique(relative_hgt_ext$COG[order(-relative_hgt_ext$relative_HGT)]))
plot_all_HGT_all_Vert <- ggplot(relative_hgt_ext, aes(x = COG, y = relative_HGT, group = 1)) + geom_line()

## It would be good to see how this changes across penalty ranges ##
extract_hgt <- final_table_hgt[1:19, c(-1:-2,-ncol(final_table_hgt))]
vert_rearrange <- cbind(final_table_vert[,-(ncol(final_table_hgt)-1):-ncol(final_table_hgt)], final_table_vert[,ncol(final_table_hgt):(ncol(final_table_hgt)-1)])
extract_vert <- vert_rearrange[1:19, c(-1:-2,-(ncol(vert_rearrange)))]
x <- extract_hgt / extract_vert
x_named <- cbind(COG_classes[1:19], x)
names(x_named)[1] <- "COG"

x_named$COG <- factor(x_named$COG, levels = unique(x_named$COG[order(-x_named$All_HGT)]))
x_named.molten <- melt(data.frame(as.matrix(x_named), row.names = NULL), id.vars = "COG")
x_named.molten$value = as.numeric(x_named.molten$value)

if (IG_included == TRUE) {
	title_IG <- "(IG only groups included)"
} else {
	title_IG <- "(IG only groups excluded)"
}


plot_relative_HGT <- ggplot(x_named.molten, aes(x = COG, y = value, group = variable, colour = variable)) + geom_line(size = 1.3) + scale_y_continuous("Relative HGT: Fraction of all groups assigned to HGT / fraction assigned to Vertical") + ggtitle(paste0("HGT signal by COG ", title_IG)) + scale_color_discrete(name = "HGT Penalties") + theme(plot.title = element_text(size = 20), legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 14)) + guides(colour = guide_legend(title.hjust = 0.5))

###########################################################################
## For each conditions (each penalty - T +/- IG) take the 
top_HGT_level <- paste0(penalty, "_", IG_ID)
COG_relative_col <- x_named[, c(1, ncol(x_named))]
names(COG_relative_col)[2] <- top_HGT_level
assign(as.character(top_HGT_level), COG_relative_col)

###########################################################################
## This section is only run once - it combines the outputs of several runs of the above ##
# library(RColorBrewer)

# MergeByCOG <- function(x, y){
# 	df <- merge(x, y, by = "COG", all = TRUE, sort = FALSE)
# 	return(df)
# }
# ## Merge the output of all the T (IG + NO_IG) values into a single able ##
# merged_relative_HGT <- Reduce(MergeByCOG, list(T4_IG, T4_NO_IG, T5_IG, T5_NO_IG, T6_IG, T6_NO_IG, T10_IG, T10_NO_IG, T20_IG, T20_NO_IG))

# ## Melt and make the y value a continuous variable ##
# merged_relative_HGT.molten <- melt(data.frame(as.matrix(merged_relative_HGT), row.names = NULL), id.vars = "COG")
# merged_relative_HGT.molten$value = as.numeric(merged_relative_HGT.molten$value)

# ## Set up a colour ramp palette ##
# col_ramp <- brewer.pal(9,"YlOrRd")
# col_palette <- colorRampPalette(col_ramp[1:9])(10)

# ## Make the plot ##
# plot_merged_relative_HGT <- ggplot(merged_relative_HGT.molten, aes(x = COG, y = value, group = variable, color = variable)) + 
# 		geom_point(size = 4) +
# 		scale_y_continuous("Relative HGT: Fraction of all groups assigned to HGT / fraction assigned to Vertical") +
# 		ggtitle(paste0("HGT signal by COG where HGT is true across all penalties - ", merged_relative_HGT_title)) + 
# 		theme(panel.background = element_rect(fill = "gray55"), panel.grid.major = element_line(color = "grey"), plot.title = element_text(size = 20), legend.key = element_rect(fill = "gray55"), legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 14, colour = "black")) +
# 		guides(colour = guide_legend(title.hjust = 0.5)) +
# 		scale_color_manual(values = col_palette, name = "HGT Penalty Ranges")

# plot_merged_relative_HGT

###########################################################################
## Plot the refined / hgt / vertical / total graph ##

###########################################################################
## Merge the values into a single table ##
merged_refined <- Reduce(MyMerge, ordered_hgt_refined_tbl_list)

# Convert any NAs to 0 #
merged_refined[is.na(merged_refined)] <- 0

## Output has the COG labels as rownames -> remove rownames and add a column corresponding to the COG functional categories ##
COG_classes <- rownames(merged_refined)

rownames(merged_refined) <- NULL

merged_refined <- cbind(COG_classes, merged_refined)

## Rearrange the columns such that All_HGT & All_Vertical appear at the end (independent of column / penalty number) ##
merged_refined <- cbind(merged_refined[,-3:-4], merged_refined[,3:4])


###########################################################################
## Make the plots ##
as_fraction_total_refined <- sweep(merged_refined[,-1],2,colSums(merged_refined[,-1]), "/")*100

final_table_refined <- cbind(COG = merged_refined[,1], as_fraction_total_refined)
final_refined.molten <- melt(final_table_refined, value.name="Percentage.of.group", variable.name="COG.category", na.rm = TRUE)

final_refined.molten.ordered <- final_refined.molten
final_refined.molten.ordered$COG <- factor(final_refined.molten$COG, levels = unique(final_refined.molten$COG[order(-final_refined.molten$Percentage.of.group)]))


plot_refined <- ggplot(data = final_refined.molten.ordered, aes(x = COG.category, y = Percentage.of.group, group = COG, fill = COG.category)) + 
			geom_bar(position = "dodge", stat = "identity") + 
			facet_wrap("COG", nrow = 1) + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			ggtitle(paste0("Proportion of putative HGT events per COG - ", merged_relative_HGT_title, " with penatly = ", maximum_penalty))



## Only total, external HGT and vertical ##
final_table_refined_2 <- final_table_refined[-4]
final_table_refined_2$COG <- droplevels(final_table_refined_2$COG)

# Drop unused levels #

final_table_refined_2$COG <- factor(final_table_refined_2$COG, levels = unique(final_table_refined_2$COG[order(-final_table_refined_2$Total)]))
names(final_table_refined_2)[2:4]<- c("All groups", "Predicted HGT", "Predicted Vertical")



final_refined_2.molten <- melt(final_table_refined_2, value.name="Percentage.of.group", variable.name="COG.category", na.rm = TRUE, factorsAsStrings = FALSE)
final_refined_2.molten$Percentage.of.group <- as.numeric(final_refined_2.molten$Percentage.of.group)

plot_refined_2 <- ggplot(data = final_refined_2.molten, aes(x = COG.category, y = Percentage.of.group, group = COG, fill = COG.category)) + 
			geom_bar(position = "dodge", stat = "identity") + 
			facet_wrap("COG", nrow = 1) + 
			theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
			ylab("Fraction of set assigned to each COG Category") +
			xlab("COG Category") +
			ggtitle(paste0("Proportion of putative HGT events per COG - ", merged_relative_HGT_title, " with penatly = ", maximum_penalty)) +
			scale_fill_discrete(name="COG Category") +
			scale_y_continuous(expand = c(0, 0), limits = c(0, 41))

devSVG(file = "/users/aesin/Dropbox/LSR/Presentation/COG_plot.svg")

ggplot(data = final_refined_2.molten, aes(x = COG.category, y = Percentage.of.group, group = COG, fill = COG.category)) + 
			geom_bar(position = "dodge", stat = "identity") + 
			facet_wrap("COG", nrow = 1) + 
			theme(axis.line.x = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
			ylab("Fraction of set assigned to each COG Category") +
			xlab("COG Category") +
			ggtitle(paste0("Proportion of putative HGT events per COG - ", merged_relative_HGT_title, " with penatly = ", maximum_penalty)) +
			scale_fill_discrete(name="COG Category") +
			scale_y_continuous(expand = c(0, 0), limits = c(0, 41))
dev.off()



# ###########################################################################

# ## Chi-squared tests ##
# ## We remove the B and Z categories as uninformative ##

# # HGT vs Vertical #
# hgt_vert <- merged_hgt[1:19, (ncol(merged_hgt)-1):(ncol(merged_hgt))]
# chisq.test(hgt_vert)
# # p-value = 1.5498e-35

# # Intra-HGT comparison #
# hgt_3_4 <- merged[1:19, 3:4]
# hgt_4_5 <- merged[1:19, 4:5]
# hgt_3_5 <- merged[1:19, c(3, 5)]
# hgt_3_6 <- merged[1:19, c(3, 6)]

# chisq.test(hgt_3_4); # p-value 0.9655
# chisq.test(hgt_4_5); # p-value 0.9996
# chisq.test(hgt_3_5); # p-value 0.6622
# chisq.test(hgt_3_6); # p-value 0.7583

# # Compare the distribution of S-categorized groups #
# s_not_s <- merged[1,]
# not_s <- merged[-1,]
# s_not_s <- rbind(s_not_s, data.frame(COG="Not_S",t(colSums(not_s[,-1]))))
# s_not_s_hgt_vert <- s_not_s[, (ncol(s_not_s)-1):(ncol(s_not_s))]

# ## There is no significant difference in the distribution between the S and not-S categories
# chisq.test(s_not_s_hgt_vert); #p-value 0.6462













