library(plyr)
library(ggplot2)
library(reshape2)
library(RSvgDevice)


# Get both the stats files. These contain,for each group, the number of total sequences, and the number of geobacilli sequences #
setwd("/users/aesin/Desktop/Geo_analysis/Geo_v_all/2.0/-10")
df_10 <- read.table("-10_group_fastas_stats.tsv")
names(df_10) <- c("Group", "Resolved_at_10", "Geobacs_at_10")

setwd("~/Desktop/Geo_analysis/Geo_v_all/2.0/Stats")
df_comb <- read.table("Combined_group_fastas_stats.tsv", header = TRUE, sep = "\t")
names(df_comb) <- c("Group", "Combined", "Geobacs_in_combined")


df_total <- subset(df_10, select=-Geobacs_at_10)
df_comb_trim <- subset(df_comb, select=Combined)
df_comb_trim <- df_comb_trim[order(-df_comb_trim$Combined),]

df_total$Combined <- df_comb_trim

## First group at which there are fewer than 50 sequences ##
cutoff <- head(df_total[(df_total$Combined == 49), 1], n = 1)

df_total.molten <- melt(df_total, id.vars = "Group", variable.name = "Reconstruction_type", value.name = "Number_of_sequences")


devSVG(file = "/users/aesin/Dropbox/LSR/Report/Figures/Group_size_distribution.svg", width = 15, height = 10)

df_plot <- ggplot(data = df_total.molten, aes(x = Group, y = Number_of_sequences, color = Reconstruction_type)) + 
		geom_line(size = 1) + 
		geom_vline(aes(xintercept = cutoff), color = "red4", linetype = 3, show.legend = T) +
		scale_linetype_manual(name = "Cutoff", labels = c("Minimum_50"), values = c("Minimum_50" = 1)) +
		ggtitle("Something here")
		# theme_classic(base_size = 20)
df_plot

dev.off()


# ####

# # Make a plot to show that the difference in Geobacillus between the Combined and the -10 dataset is not distributed evenly accross the families, but is instead restricted to a small number of groups. These probably represent convergent groups... #

# ####

# df_10_trim<-subset(df_10, select=-Resolved_at_10)
# df_total<-df_10_trim
# df_comb_trim<-subset(df_comb, select=-Combined)
# df_total<-merge(df_10_trim, df_comb_trim, by = "Group")
# Geobac_difference<-sort(df_total$Geobacs_at_10 - df_total$Geobacs_in_combined, decreasing = T)
# Group<-df_total$Group
# df_total_test<-data.frame(Group, Geobac_difference)
# #df_total_test<-df_total_test[order(-df_total_test$Geobac_difference), ]

# #df_total_test<-melt(df_total_test, id.vars = "Group")

# df_plot<-ggplot(data = df_total_test, aes(x = Group, y = Geobac_difference)) + geom_line(color = "red4")
# #df_plot<-df_plot + geom_vline(aes(xintercept = as.numeric(cutoff)))
# df_plot<-df_plot + labs(title = "Geobacillus inclusion in groups", y = "Difference in Geobacillus sequences")
# df_plot<-df_plot + theme_classic(base_size = 20)
# df_plot

# ggsave(filename = "Missing_Geobacillus_distribution.png", plot = df_plot)

# ####

# # Zoom on the above plot #

# ####

# library(plyr)
# library(ggplot2)
# library(reshape2)


# # Get both the stats files. These contain,for each group, the number of total sequences, and the number of geobacilli sequences #
# setwd("/users/aesin/Desktop/Geo_analysis/Geo_v_all/2.0/Stats")
# missing_geo<-read.table("Missing_geobac_families.tsv", sep="\t", header=T)

# missing_geo_number<-subset(missing_geo, select=-Group.number)
# missing_geo_order<-missing_geo_number[order(-missing_geo_number[,1]), ]
# Group<-1:length(missing_geo_order)

# missing_geo_df<-data.frame(Group, missing_geo_order)

# df_plot<-ggplot(data = missing_geo_df, aes(x = Group, y = missing_geo_order)) + geom_line(color = "red4")
# df_plot<-df_plot + labs(title = "Distribution of missing Geobacilli across groups", y = "Number of missing sequences")
# df_plot<-df_plot + theme_classic(base_size = 20)
# df_plot

# ggsave(filename = "Missing_Geobacillus_distribution_ZOOM.png", plot = df_plot)

# ## Plot the density of Geobacillis seqeuences / group ##
# setwd("/users/aesin/desktop/Geo_analysis/Geo_v_all/2.0/Stats")
# group_geo_number <- read.table(file = "Combined_group_fastas_stats.tsv", header = TRUE, sep = "\t")
# num_groups <- nrow(group_geo_number)
# dg <- density(x = group_geo_number$Total.Geobacillus, bw = "sj")
# pdf("Density of Geobacillus sequences across the orthologous groups.pdf")
# plot(dg, main = "Density of Geobacillus sequences across the orthologous groups", xlab = paste0("Number of Geobacillus sequences / group [N = ", num_groups, "]"), col = "firebrick3")
# dev.off()

# setwd("/Users/aesin/Desktop/Geo_analysis/Geo_v_all/2.0/-10")
# group_geo_number <- read.table(file = "-10_group_fastas_stats.tsv", header = FALSE, sep = "\t")
# num_groups <- nrow(group_geo_number)
# dg <- density(x = group_geo_number$V3, bw = "sj")
# pdf("Density of Geobacillus sequences across the orthologous groups.pdf")
# plot(dg, main = "Density of Geobacillus sequences across the orthologous groups", xlab = paste0("Number of Geobacillus sequences / group [N = ", num_groups, "]"), col = "firebrick3")
# dev.off()


