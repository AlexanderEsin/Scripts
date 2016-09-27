roundUp <- function(x,to=100)
{
  to*(x%/%to + as.logical(x%%to))
}

setwd("/users/aesin/desktop/Geo_analysis/Geo_ortholog_nucl/Functional_annotation")
dir.create("Graphs", showWarnings = F)
core_name <- "All_mapped_groups"

## Narrow categories BAR ##
redun_cog_l <- read.table(file = "Narrow_COG_list.txt", sep = "\n")

redun_cog_tab <- table(redun_cog_l)
redun_cot_tab_ordered <- redun_cog_tab[order(-redun_cog_tab)]

pdf(paste0("Graphs/", paste0(core_name, "_functional_COGs.pdf")))
par(cex = 0.8)
barplot(redun_cot_tab_ordered, axes = F, col = "firebrick3", main = "Cluster COG IDs (narrow)")
upper_limit <- roundUp(max(redun_cot_tab_ordered), 100)
axis(2, seq(0, upper_limit, by = 100), xpd = TRUE); #800 by 100 for ALL, 60 by 10 for Non_part_hyp
dev.off()

## Broad categories BAR ##

broad_l <- read.table(file = "Broad_COG_list.txt", sep = "\n")
broad_tbl <- table(broad_l)
broad_tbl_o <- broad_tbl[order(-broad_tbl)]

pdf(paste0("Graphs/", paste0(core_name, "_broad_functional_cats.pdf")))
par(mar = c(10,4,8,2))
barplot(broad_tbl_o, axes = F, col = "firebrick3", cex.names = 0.6, las = 2, main = "Cluster functional categories")
upper_limit <- roundUp(max(broad_tbl_o), 500)
axis(2, seq(0, upper_limit, by = 500), xpd = T); #800 by 100 for ALL, 90 by 10 for Non_part_hyp
dev.off()

## Broad categories PIE ##

pdf(paste0("Graphs/", paste0(core_name, "_broad_functional_cats_PIE.pdf")))
par(mar = c(5,4,4,8), cex = 0.7)
pie(broad_tbl, col=terrain.colors(4))
dev.off()