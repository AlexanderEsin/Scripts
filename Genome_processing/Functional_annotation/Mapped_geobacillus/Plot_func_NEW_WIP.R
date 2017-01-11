library(ggplot2)
library(reshape2)
library(gridExtra)
library(RSvgDevice)
library(stringr)
library(RColorBrewer)

###########################################################################
## Functions ##

ReadCogTable <- function(file, name) {
	cog_tab <- table(read.table(file = file, sep = "\n"))
	cog_tab_ordered <- as.data.frame(cog_tab[order(-cog_tab)])
	colnames(cog_tab_ordered)[1] <- name
	return(cog_tab_ordered)
}

MyMerge <- function(x, y){
    df <- merge(x, y, by = "row.names", all = TRUE, sort = FALSE)
    rownames(df)  <- df$Row.names
    df$Row.names  <- NULL
    return(df)
}

MergeByCOG <- function(x, y){
	df <- merge(x, y, by = "COG", all = TRUE, sort = FALSE)
	return(df)
}

FormatPenaltyTables <- function(df){
	## Convert any NAs to 0 ##
	df[is.na(df)] <- 0
	## Output has the COG labels as rownames -> remove rownames and add a column corresponding to the COG functional categories ##
	COG_classes <- rownames(df)
	rownames(df) <- NULL
	df <- cbind(COG_classes, df)
	## Rearrange the columns such that All_HGT & All_Vertical appear at the end (independent of column / penalty number) ##
	df <- cbind(df[,-3:-4], df[,3:4])
	return(df)
}

FormatLongDistTables <- function(df) {
	## Convert any NAs to 0 ##
	df[is.na(df)] <- 0
	## Output has the COG labels as rownames -> remove rownames and add a column corresponding to the COG functional categories ##
	COG <- rownames(df)
	rownames(df) <- NULL
	df <- cbind(COG, df)
	return(df)
}

PerPenaltyRawCogNums <- function(formatted_penalty_tbl, name) {
	## Extract just those columns corresponding to the raw predictions ##
	raw_penalty_cols <- colnames(formatted_penalty_tbl[ ,3:(ncol(formatted_penalty_tbl)-2)])
	assign(paste0(name, "_raw_pen_col_names"), raw_penalty_cols, envir = .GlobalEnv)
	## Take a sum of the predictions across all COG classes ##
	column_sums <- colSums(formatted_penalty_tbl[raw_penalty_cols])
	## Rename the columns ##
	penalty_col_names <- gsub(pattern = paste0(name, "_"),"", raw_penalty_cols)
	## Put into a df ##
	HGT_num_df <- data.frame("Penalty" = penalty_col_names, "Col_sums" = column_sums, row.names = NULL)
	colnames(HGT_num_df)[2] <- paste0("Groups_w_", name)
	return(HGT_num_df)
}

FormattedTblForPlot <- function(formatted_penalty_tbl, name) {
	as_fraction_total <- sweep(formatted_penalty_tbl[,-1],2,colSums(formatted_penalty_tbl[,-1]), "/")*100

	as_fraction_tbl <- cbind(COG = formatted_penalty_tbl[,1], as_fraction_total)
	global_name <- paste0(name, "_as_fraction_tbl")
	assign(global_name, as_fraction_tbl, envir = .GlobalEnv)

	as_fraction_molten <- melt(as_fraction_tbl, value.name="Percentage.of.group", variable.name="COG.category", na.rm = TRUE)

	as_fraction_molten_ordered <- as_fraction_molten

	as_fraction_molten_ordered$COG <- factor(as_fraction_molten$COG, levels = unique(as_fraction_molten$COG[order(-as_fraction_molten$Percentage.of.group)]))
	return(as_fraction_molten_ordered)
}

PlotBarsRaw <- function(molten_table, name) {

	plot <-	ggplot(data = molten_table, aes(x = COG.category, y = Percentage.of.group, group = COG)) + 
			geom_bar(position = "dodge", stat = "identity", colour = "black", fill = "white") + 
			facet_wrap("COG", nrow = 1) + 
			theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
			ggtitle(paste0("Proportion of putative HGT events per COG over a penalty range - ")) + 
			geom_line(data = subset(molten_table, COG.category %in% eval(parse(text = paste0(name, "_raw_pen_col_names")))), colour = "red", size = 1)

	return(plot)
}

PlotRelativeEnrichmentBoxes <- function(prop_list, name) {

	merged_prop_df <- Reduce(MergeByCOG, prop_list)
	merged_prop_molten <- melt(data.frame(as.matrix(merged_prop_df), row.names = NULL), id.vars = "COG")
	merged_prop_molten$value = as.numeric(merged_prop_molten$value)

	## Set up a colour ramp palette ##
	col_ramp <- brewer.pal(9,"YlOrRd")
	col_palette <- colorRampPalette(col_ramp[1:9])(length_orig_pen_list)

	plot_relative_HGT <- ggplot(merged_prop_molten, aes(x = COG, y = value, group = variable, color = variable)) + 
			geom_point(size = 4) +
			scale_y_continuous("Relative HGT: Fraction of all groups assigned to HGT / fraction assigned to Vertical") +
			ggtitle(paste0("COG enrichment for HGT for the ", name, " set")) + 
			theme(panel.grid.major = element_line(color = "grey"), plot.title = element_text(size = 20), legend.key = element_rect(fill = "gray55"), legend.text = element_text(size = 16), legend.title = element_text(size = 18), axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14), axis.text = element_text(size = 14, colour = "black")) +
			guides(colour = guide_legend(title.hjust = 0.5)) +
			scale_color_manual(values = col_palette, name = "HGT Penalty Ranges")

	## Draw boxes ##
	coordinate_df <- data.frame(x1 = numeric(), x2 = numeric(), y1 = numeric(), y2 = numeric())

	for (i in 1:length(merged_prop_df$COG)) {
		COG_needed <- levels(merged_prop_molten$COG)[i]
		row_number <- which(merged_prop_df$COG == COG_needed)
		#print(COG_needed)

		max_y <- max(merged_prop_df[row_number,-1])
		#print(max_y)
		min_y <- min(merged_prop_df[row_number,-1])
		#print(min_y)

		## Add vertical lines ##
		coordinate_df[nrow(coordinate_df)+1, ] <- c((i - 0.2), (i - 0.2), min_y, max_y)
		coordinate_df <- rbind(coordinate_df, c((i + 0.2), (i + 0.2), min_y, max_y))

		## Add horizontal lines ##
		coordinate_df <- rbind(coordinate_df, c((i - 0.2), (i + 0.2), min_y, min_y))
		coordinate_df <- rbind(coordinate_df, c((i - 0.2), (i + 0.2), max_y, max_y))

	}

	# devSVG(file = "file_name", width = 15, height = 10)

	plot_relative_HGT_boxes <- plot_relative_HGT + geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = coordinate_df, inherit.aes = FALSE)
	return(plot_relative_HGT_boxes)
}


###########################################################################
IG_included = TRUE
All_scenarios = FALSE

#penalty_list <- c(3, 4, 5, 6, 8, 10, 20)
penalty_list <- c(3, 4, 5, 6, 8, 10)
length_orig_pen_list <- length(penalty_list)
const_prop_list <- vector("list", length_orig_pen_list - 1)
long_dist_prop_list <- vector("list", length_orig_pen_list - 1)
no_more_hor_list <- vector("list", length_orig_pen_list - 1)

if (All_scenarios == TRUE) {
	Scenario_ID <- "All_scenarios/"
	merged_relative_HGT_title <- "All_scenarios"
} else {
	Scenario_ID <- "Scenarios_1_2/"
	merged_relative_HGT_title <- "Scenarios 1 & 2 ONLY"
}

if (IG_included == TRUE) {
	IG_ID <- "Full/"
} else {
	IG_ID <- "NO_IG/"
}

direct <- "/users/aesin/Desktop/"
consistent_dir	<- paste0(direct, "Geo_analysis/HGT_Functional_annotation/Consistent_HGT/", IG_ID, Scenario_ID)
mowg_const_dir	<- paste0(direct, "Mowgli/Consistent_HGT_Vertical/", paste0("Consistent_", IG_ID), Scenario_ID)
per_penalty_dir <- paste0(direct, "Geo_analysis/HGT_Functional_annotation/Per_penalty/", IG_ID, Scenario_ID)
long_dist_dir	<- paste0(direct, "Geo_analysis/HGT_Functional_annotation/Long_distance_HGT/", IG_ID, Scenario_ID)

no_more_hor_dir	<- paste0(consistent_dir, "No_more_HGT/")

i = 1

while (length(penalty_list) > 1) {
	###########################################################################
	## Prepare the list of penalties to be queried downstream ##
	penalty_num <- length(penalty_list)

	## Prepare list for the ordered HGT tables, but also write them out to individual variables ##
	ordered_hor_table_list <- vector("list", penalty_num)
	ordered_ver_table_list <- vector("list", penalty_num)
	hor_col_counter <- 1
	ver_col_counter <- 1
	col_counter <- 1

	###############################
	## Read in total (complete set of gene families/groups/clusters) COG annotation for ALL groups ##
	setwd("/users/aesin/desktop/Geo_analysis/Geo_ortholog_nucl/Functional_annotation")
	total_cog_tab <- ReadCogTable("Narrow_COG_list.txt", "Total")

	# Add to our table lists #
	ordered_hor_table_list[[col_counter]] <- total_cog_tab
	ordered_ver_table_list[[col_counter]] <- total_cog_tab
	col_counter = col_counter + 1

	###############################
	### Read in the "consistent" sets based on our penalty list. We take the HGT at the highest penalty, and the Vertical at the LOWEST penalty ###
	consistent_folder_name <- str_c(penalty_list, collapse = "_")
	hor_const_folder <- paste0(consistent_dir, "HGT/", consistent_folder_name)
	ver_const_folder <- paste0(consistent_dir, "Vertical/", consistent_folder_name)

	if (file.exists(hor_const_folder) == FALSE || file.exists(ver_const_folder) == FALSE) {
		stop(paste0("Either the HGT or Vertical consistent folder is missing for the penalty range: ", paste(penalty_list, collapse = " ")))
	}

	## Read in the consistent tables ##
	hor_const_cog_tab <- ReadCogTable(paste0(hor_const_folder, "/Narrow_COG.tsv"), "Consistent_HGT")
	ver_const_cog_tab <- ReadCogTable(paste0(ver_const_folder, "/Narrow_COG.tsv"), "Consistent_Ver")

	## Enter them into the table lists ##
	ordered_hor_table_list[[col_counter]] <- hor_const_cog_tab
	ordered_ver_table_list[[col_counter]] <- hor_const_cog_tab

	ordered_hor_table_list[[col_counter+1]] <- ver_const_cog_tab
	ordered_ver_table_list[[col_counter+1]] <- ver_const_cog_tab
	col_counter = col_counter + 2

	###############################
	### Read in the individual HGT and vertical tables ###
	for (penalty in penalty_list) {
		## Get the path name for each penalty and check it exists ##
		penalty_folder <- paste0("T", penalty)
		penalty_folder_path <- paste0(per_penalty_dir, penalty_folder)

		if (file.exists(penalty_folder_path) == FALSE) {
			stop(paste0("The per-penalty folder for penalty: ", penalty), " is missing")
		}
		
		## Prepare the column names and read in the per-penalty "raw" COGs ##
		colname_hor <- paste0("HGT_", penalty)
		colname_ver <- paste0("Ver_", penalty)
		hor_penalty_tab <- ReadCogTable(paste0(penalty_folder_path, "/HGT_narrow_COG.tsv"), colname_hor)
		ver_penalty_tab <- ReadCogTable(paste0(penalty_folder_path, "/Vertical_narrow_COG.tsv"), colname_ver)

		## Enter them into the table lists ##
		ordered_hor_table_list[[col_counter]] <- hor_penalty_tab
		ordered_ver_table_list[[col_counter]] <- ver_penalty_tab

		col_counter = col_counter + 1

		## Assign to individual tables ##
		name_hor <- paste0("hor_cog_tab_", penalty)
		name_ver <- paste0("ver_cog_tab_", penalty)
		assign(name_hor, FormatLongDistTables(hor_penalty_tab))
		assign(name_ver, FormatLongDistTables(ver_penalty_tab))

	}

	###############################
	## Merge the values into a single table ##
	merged_hor_df <- Reduce(MyMerge, ordered_hor_table_list)
	merged_ver_df <- Reduce(MyMerge, ordered_ver_table_list)

	## Format the merged tables ##
	formatted_hor_df <- FormatPenaltyTables(merged_hor_df)
	formatted_ver_df <- FormatPenaltyTables(merged_ver_df)

	###############################
	## Make a df containing the total number of COGs assigned to the per-penalty RAW data ##
	hor_raw_num_df <- PerPenaltyRawCogNums(formatted_hor_df, "HGT")
	ver_raw_num_df <- PerPenaltyRawCogNums(formatted_ver_df, "Ver")

	###############################
	## Prepare the data in the formatted tables for plotting ##

	hor_raw_const_for_plot <- FormattedTblForPlot(formatted_hor_df, "HGT")
	ver_raw_const_for_plot <- FormattedTblForPlot(formatted_ver_df, "Ver")

	plot_hor_bars <- PlotBarsRaw(hor_raw_const_for_plot, "HGT")
	plot_ver_bars <- PlotBarsRaw(ver_raw_const_for_plot, "Ver")

	###########################################################################
	### Prepare a plot showing how many groups were reconciled, and within those how many were predicged as HGT or Vertical (RAW) ###

	## Read in the table showing how many groups were reconciled per penalty ##
	reconciled_per_penalty_name	<- paste0("t", paste(penalty_list, collapse = "_t"), ".tsv")
	reconciled_per_penalty_df	<- read.table(paste0(mowg_const_dir, "Groups_reconciled_per_penalty_", reconciled_per_penalty_name), header = TRUE)
	colnames(reconciled_per_penalty_df)[1] <- "Penalty"

	## Merge with the tables showing how many RAW HGT and Vertical predictions were made per penalty ##
	merged_num_reconciled_df <- merge(reconciled_per_penalty_df, hor_raw_num_df, by = "Penalty")
	merged_num_reconciled_df <- merge(merged_num_reconciled_df, ver_raw_num_df, by = "Penalty")
	num_reconciled_molten	 <- melt(merged_num_reconciled_df, id.vars = "Penalty")

	## Plot a point graph ##
	if (i == 1) {
		plot_vert_hgt_total_num <- ggplot(num_reconciled_molten, aes(x = Penalty, y = value, colour = variable)) + geom_line() + geom_point()
	}

	###########################################################################
	### For the consistent dataset at the top penalty for each penalty set, get the relative HGT enrichment by : ###
	### Fraction of HGT set assigned to COG / Fraction of Vertical set assigned to COG ###
	HGT_as_fraction_tbl$Const_hor_ver_prop <- (HGT_as_fraction_tbl$Consistent_HGT / HGT_as_fraction_tbl$Consistent_Ver)

	Const_hor_ver_prop <- HGT_as_fraction_tbl[1:19, c(1, ncol(HGT_as_fraction_tbl))]
	const_prop_name <- paste0("Const_prop_", str_sub(IG_ID, 1, -2), "_", max(penalty_list))
	colnames(Const_hor_ver_prop)[2] <- const_prop_name

	## Add to global list for plotting below ##
	const_prop_list[[length(penalty_list)-1]] <- Const_hor_ver_prop

	plot_consistent_proportion <- ggplot(Const_hor_ver_prop, aes(x = COG, y = Const_hor_ver_prop, group = 1)) + geom_point()

	###########################################################################
	### Read in the long-distance COG annotations at each maximum penalty ###
	maximum_penalty <- max(penalty_list)
	long_dist_COG_filename <- paste0(long_dist_dir, "Narrow/T", maximum_penalty, "_narrow_COG.tsv")

	## Check the relevant file exists ##
	if (file.exists(long_dist_COG_filename) == FALSE) {
		stop(paste0("The long-distance HGT COG annotations for penalty ", maximum_penalty, "does not exist. Exiting..."))
	}

	## Read it in ##
	long_dist_cog_tab <- ReadCogTable(long_dist_COG_filename, "Long_distance_HGT")

	## Export it to global ##
	for_global_ld <- FormatLongDistTables(long_dist_cog_tab)
	global_name_ld <- paste0("Long_distance_cog_", maximum_penalty)
	assign(global_name_ld, for_global_ld)
	
	## Make the values a fraction of the total COGs assigned ##
	long_dist_fract_tab <- FormatLongDistTables(sweep(long_dist_cog_tab, 2, colSums(long_dist_cog_tab), "/") * 100)

	## Make the values a proportional enrichment (i.e. fraction long dist HGT / fraction Vertical) - as with the Constant data above ##
	HGT_as_frac_tbl_ncol <- ncol(HGT_as_fraction_tbl)
	long_dist_as_fract_tbl <- MergeByCOG(HGT_as_fraction_tbl[,c(-(HGT_as_frac_tbl_ncol - 2),-HGT_as_frac_tbl_ncol)], long_dist_fract_tab)
	long_dist_as_fract_tbl$Long_hor_ver_prop <- (long_dist_as_fract_tbl$Long_distance_HGT / long_dist_as_fract_tbl$Consistent_Ver)

	## Extract the proportional enrichment column together with the COGs and add it to the global list ##
	long_dist_hor_ver_prop <- long_dist_as_fract_tbl[1:19, c(1, ncol(long_dist_as_fract_tbl))]
	long_dist_prop_name <- paste0("Long_dist_prop_", str_sub(IG_ID, 1, -2), "_", max(penalty_list))
	colnames(long_dist_hor_ver_prop)[2] <- long_dist_prop_name

	long_dist_prop_list[[length(penalty_list)-1]] <- long_dist_hor_ver_prop


	###########################################################################
	## Read in the functional annotations of those groups that are lost in the consistency check from each penalty ##
	no_more_hor_penalties	<- penalty_list[1:(length(penalty_list) - 1)]
	no_more_hor_path		<- paste0(no_more_hor_dir, paste(no_more_hor_penalties, collapse = "_"))

	if (file.exists(paste0(no_more_hor_path, "/Narrow_COG.tsv")) == FALSE) {
		stop(paste0("The no more HGT COG annotations for penalty range:", no_more_hor_penalties, "do not exist. Exiting..."))
	}

	## Read it in ##
	no_more_hor_cog_tab <- ReadCogTable(paste0(no_more_hor_path, "/Narrow_COG.tsv"), "No_more_HGT")

	global_name <- paste0("RAW_", colnames(no_more_hor_cog_tab[1]), "_", as.character(max(no_more_hor_penalties)))
	no_mor_hor_global <- FormatLongDistTables(no_more_hor_cog_tab)
	colnames(no_mor_hor_global)[2] <- global_name

	assign(global_name, no_mor_hor_global)

	## Make the values a fraction of the total COGs assigned ##
	total_fract_tab <- FormatLongDistTables(sweep(total_cog_tab, 2, colSums(total_cog_tab), "/") * 100)
	no_more_hor_fract_tab	<- FormatLongDistTables(sweep(no_more_hor_cog_tab, 2, colSums(no_more_hor_cog_tab), "/") * 100)

	merged_total_no_more	<- MergeByCOG(total_fract_tab, no_more_hor_fract_tab)
	merged_total_no_more$Normalized	<- (merged_total_no_more[,3] / merged_total_no_more[,2])
	merged_total_no_more[is.na(merged_total_no_more)] <- 0

	colnames(merged_total_no_more)[4] <- paste0("Normal_", colnames(no_more_hor_cog_tab[1]), "_", as.character(max(no_more_hor_penalties)))
	merged_total_no_more <- merged_total_no_more[,-2:-3]

	
	no_more_hor_list[[length(penalty_list)-1]] <- merged_total_no_more



	###########################################################################
	penalty_list <- penalty_list[-length(penalty_list)]
	i = i + 1
}

# total_cog_col_tab <- FormatLongDistTables(total_cog_tab)
# tot_4 <- MergeByCOG(RAW_No_more_HGT_5, Long_distance_cog_4)
# tot_4[is.na(tot_4)] <- 0
# chisq.test(tot_4[,2:3])

total_cog_col_tab <- FormatLongDistTables(total_cog_tab)
tot_3 <- MergeByCOG(RAW_No_more_HGT_3, Long_distance_cog_4)
tot_3[is.na(tot_3)] <- 0
chisq.test(tot_3[,2:3])

tot_4 <- MergeByCOG(RAW_No_more_HGT_4, Long_distance_cog_5)
tot_4[is.na(tot_4)] <- 0
chisq.test(tot_4[,2:3])

tot_5 <- MergeByCOG(RAW_No_more_HGT_5, Long_distance_cog_6)
tot_5[is.na(tot_5)] <- 0
chisq.test(tot_5[,2:3])

tot_6 <- MergeByCOG(RAW_No_more_HGT_6, Long_distance_cog_8)
tot_6[is.na(tot_6)] <- 0
chisq.test(tot_6[,2:3])

xx <- Reduce(MergeByCOG, no_more_hor_list)
xx_melt <- melt(data.frame(as.matrix(xx[1:17,]), row.names = NULL), id.vars = "COG")
xx_melt$value <- as.numeric(xx_melt$value)
xx_melt$COG <- factor(xx_melt$COG, levels = unique(xx_melt$COG[order(-xx_melt$value)]))
ggplot(xx_melt, aes(x = COG, y = value, group = variable)) + geom_point(aes(color = variable)) + geom_line(aes(color = variable))

# # plot_const_relative_enrichment <- PlotRelativeEnrichmentBoxes(const_prop_list, "Consistent")
# # plot_longd_relative_enrichment <- PlotRelativeEnrichmentBoxes(long_dist_prop_list, "Long Distance")



# ###############################
# ## Compare the consistent and long-distance sets ##
# const_prop_df <- Reduce(MergeByCOG, const_prop_list)
# longd_prop_df <- Reduce(MergeByCOG, long_dist_prop_list)

# x <- MergeByCOG(const_prop_df, longd_prop_df)

# y <- (x[,5:7]) - (x[,2:4])
# y <- cbind(x[,1], y)
# names(y)[1] <- "COG"
# y_melt <- melt(data.frame(as.matrix(y), row.names = NULL), id.vars = "COG")
# y_melt$value <- as.numeric(y_melt$value)
# ggplot(y_melt, aes(x = COG, y = value, color = variable)) + geom_point() + scale_y_continuous("Relative HGT: Fraction of all groups assigned to HGT / fraction assigned to Vertical")

# ###############################


###########################################################################
###########################################################################
###########################################################################


































# ###########################################################################
# ## Read in the table corresponding to the refined dataset containing transfers from only OUTSIDE the IG ##



