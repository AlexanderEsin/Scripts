library(gtools)
library(logspline)
library(circular)
library(polyclip)
library(stringr)

scale_range <- function(vector_x) {
	(vector_x - min(vector_x)) / (max(vector_x) - min(vector_x))
}

read_in_position_data <- function(file_name) {
	pos_all_data <- read.table(file = file_name, header = FALSE, sep = "\n")
	pos_all_df <- as.data.frame(mixedsort(pos_all_data$V1))
	colnames(pos_all_df) <- "V1"
	return(pos_all_df)
}

per_species_plot <- function(species, bandwith = 400, set = "long", bg = "all", window = 500, penalty) {

	###############################
	## Position for all datasets ##
	###############################
	master_dir <- paste0("/Users/aesin/Desktop/Geo_analysis/HGT_position/For_circular/T", penalty)

	setwd(paste0(master_dir, "/Per_species_all/Relative_start"))
	pos_all_df		<- read_in_position_data(grep(species, dir(), value = T))
	setwd(paste0(master_dir, "/Per_species_const/Relative_start"))
	pos_const_df 	<- read_in_position_data(grep(species, dir(), value = T))
	setwd(paste0(master_dir, "/Per_species_long/Relative_start"))
	pos_long_df 	<- read_in_position_data(grep(species, dir(), value = T))
	setwd(paste0(master_dir, "/Per_species_vert/Relative_start"))
	pos_vert_df		<- read_in_position_data(grep(species, dir(), value = T))

	## Convert the genome coordinates to circular ##
	circ_all	<- circular(pos_all_df$V1*(2*pi))
	circ_const 	<- circular(pos_const_df$V1*(2*pi))
	circ_long 	<- circular(pos_long_df$V1*(2*pi))
	circ_vert 	<- circular(pos_vert_df$V1*(2*pi))

	## Get density using von Mises ##
	circ_dens_all 	<- density.circular(circ_all, kernel = "vonmises", bw = bandwith)
	circ_dens_const <- density.circular(circ_const, kernal = "vonmises", bw = bandwith)
	circ_dens_long 	<- density.circular(circ_long, kernel = "vonmises", bw = bandwith)
	circ_dens_vert 	<- density.circular(circ_vert, kernel = "vonmises", bw = bandwith)

	## Differentially color the polygons basd on the dataset used ##
	if (set == "long") {
		circ_dens_data <- circ_dens_long
		vector_df <- pos_long_df$V1
		line_col <- rgb(0,1,0,0.5)
	} else if (set == "const") {
		circ_dens_data <- circ_dens_const
		vector_df <- pos_const_df$V1
		line_col <- rgb(0,0.4,0,0.5)
	} else if (set == "vert") {
		circ_dens_data <- circ_dens_vert
		vector_df <- pos_vert_df$V1
		line_col <- rgb(0,0,1,0.5)
	} else {
		message("Set can be one of \"long\", \"const\", or \"vert\"")
		stop()
	}

	## Assign the background dataset ##
	if (bg == "all") {
		circ_dens_bg <- circ_dens_all
	} else if (bg == "const") {
		circ_dens_bg <- circ_dens_const
	} else if (bg == "vert") {
		circ_dens_bg <- circ_dens_vert 
	} else {
		message("Set can be one of \"all\", \"const\", or \"vert\"")
		stop()
	}

	####################################
	## Process the GC enrichment data ##
	####################################
	setwd(paste0("/users/aesin/desktop/Geo_analysis/HGT_position/GC_content/Window_", window))

	gc_for_species_file <- grep(species, dir(), value = T)
	gc_data		<- read.table(gc_for_species_file, header = F, sep = "\t")
	gc_data$V1	<- circular(gc_data$V1 * 2 * pi)
	# Convert the GC-richness to AT-richness #
	gc_data$V3	<- (1 - gc_data$V2)
	# Critical values are those where AT or GC richness > 2 * sd of dataset #
	crit_val_gc	<- 2 * sd(gc_data$V2)
	crit_val_at <- 2 * sd(gc_data$V3)
	# Isolate all windows with AT > critical value #
	high_gc		<- gc_data[which(gc_data$V2 > (mean(gc_data$V2) + crit_val_gc)),]
	high_at		<- gc_data[which(gc_data$V3 > (mean(gc_data$V3) + crit_val_at)),]

	gc_dens		<- density.circular(high_gc$V1, kernel = "vonmises", bw = bandwith)
	at_dens		<- density.circular(high_at$V1, kernel = "vonmises", bw = bandwith)

	#########################
	## Plot rrna positions ##
	#########################
	setwd("/Users/aesin/Desktop/Geo_analysis/HGT_position/Ribosmal_rna_locations/Relative_start")

	rrna_pos	<- NULL
	rrna_pos	<- read_in_position_data(grep(species, dir(), value = T))
	circ_rrna	<- circular(rrna_pos$V1*(2*pi))
	#rrna_dens 	<- density.circular(circ_rrna, kernel = "vonmises", bw = 100)

	#############################################################
	########################### Plots ###########################
	#############################################################
	
	#####################
	## Nucl enrichment ##
	#####################

	# Plot AT enrichment #
	plot_nuc_enrich <- plot(at_dens, tol = 0, ylab = NA, xlab = NA, axes = F, control.circle = circle.control(type = "n"), points.plot = F, lwd = 2, shrink = 2.4, xlim = c(-1, 1), ylim = c(-1, 1), col = "black", zero = (pi/2), rotation = "clock", main = "")
	# Add GC enrichment #
	gc_dens_line 	<- lines(gc_dens, plot.info = plot_nuc_enrich, lwd = 2, shrink = 2.4, col = "chocolate1")

	## 0 line  & axes ##
	lines(x = circular(seq(0, (2 * pi), length.out = 200)), y = rep(0, 200), col = "gray60", plot.info = plot_nuc_enrich, lwd = 1.5)
	axis.circular(at = NULL, labels = c("","", "",""), zero = (pi/2), rotation = "clock", template = "none", tcl = 0.12, tcl.text = -0.4, lwd = 2, tick = T)
	

	#at_line 		<- lines(at_dens, lwd = 1, col = "black", plot.info = plot_for_at_val)
	#points(high_at$V1, plot.info = plot_for_at_val, stack = T, sep = 0.04, pch = 16, cex = 0.8, col = "red3")
	
	points(circ_rrna, plot.info = plot_nuc_enrich, stack = T, sep = 0.1, pch = 16, cex = 2, col = "mediumvioletred", bins = 500)

	#par(new=TRUE)
	#plot_for_rrna <- plot(circ_rrna, tol = 0, ylab = NA, xlab = NA, axes = F, control.circle = circle.control(type = "n"), lwd = 1.5, shrink = 1.3, col = "white", zero = (pi/2), rotation = "clock")
	#points(circ_rrna, plot.info = plot_for_rrna, stack = T, sep = 0.04, pch = 16, cex = 0.8, col = "black")

	# Make a new plot on top for density data #
	par(new=TRUE)
	circ_plot <- plot(circ_dens_all, tol = 0, ylab = NA, xlab = NA, axes = F, control.circle = circle.control(type = "n"), points.plot = F, nosort = T, lwd = 1.5, shrink = 1.5, col = "white", plot.type = "circle", zero = (pi/2), rotation = "clock", main = "")
	title(main = paste0(species, " n = ", length(vector_df), "  T = ", penalty), cex.main = 2.25)

	# Plot the two lines - background and comparison #
	bg_line <- lines(circ_dens_bg, lwd = 1, col = "black", zero = (pi/2), rotation = "clock")
	data_line <- lines(circ_dens_data, lwd= 0, col = "white", zero = (pi/2), rotation = "clock")
	
	# If the background is not "all" gene dataset, plot the "all" dataset as a black line #
	if (bg != "all") {
		all_line <- lines(circ_dens_all, lwd = 1.5, col = "black", zero = (pi/2), rotation = "clock")
	}

	# Only plot a seperate blue line for "vertical" genes if vertical is not being the set being compared to a background #
	if (set != "vert" && bg == "all") {
		vert_line <- lines(circ_dens_vert, lwd = 2, col = "blue", zero = (pi/2), rotation = "clock")
	}

	# Produce polygons #
	lapply(polyclip(A=list("x"=data_line$x, "y"=data_line$y), B=list("x"=bg_line$x, "y"=bg_line$y), op="minus"), polygon, col = line_col, border = line_col)
	lapply(polyclip(B=list("x"=data_line$x, "y"=data_line$y), A=list("x"=bg_line$x, "y"=bg_line$y), op="minus"), polygon, col = rgb(1,0,0,0.5), border = rgb(1,0,0,0.5))

	# Plot axes last #
	axis.circular(at = NULL, labels = c("Origin", "", "Terminus", ""), zero = (pi/2), rotation = "clock", template = "none", tcl = 0.12, tcl.text = -0.4, lwd = 2, tick = F, cex = 2.25)


	pt <- recordPlot()
	dev.off()

	return(pt)
	
}

##########################
## Circular per species ##
##########################

## Get assembled species as list ##
assembled_species_l <- as.character(read.table("/users/aesin/desktop/Geo_analysis/HGT_position/assembled_geo.tsv", header = F, sep = "\n")$V1)

## Window ##
window_size <- 1000
penalty <- "6"

par(mfrow = c(3, 4))
par(mar = c(0, 0, 1, 0), oma = c(0, 0, 1, 0))
par(cex = 0.6)

for (species in assembled_species_l) {
	species_plot <- per_species_plot(species, bandwith = 200, set = "long", bg = "all", window = window_size, penalty = penalty)

	replayPlot(species_plot)
}


# #########################################
# ## For circular density representation ##
# #########################################
# dev.off()
# # All genes #
# setwd("/Users/aesin/Desktop/Geo_analysis/HGT_position/For_circular/Per_species_all/")
# pos_all_df <- read_in_position_data("combined_rel_pos_assembled.tsv")

# # Consistent HGT prediction #
# setwd("/Users/aesin/Desktop/Geo_analysis/HGT_position/For_circular/Per_species_const/")
# pos_const_df <- read_in_position_data("combined_rel_pos_assembled.tsv")

# # Long distance HGT prediction #
# setwd("/Users/aesin/Desktop/Geo_analysis/HGT_position/For_circular/Per_species_long/")
# pos_long_df <- read_in_position_data("combined_rel_pos_assembled.tsv")

# # Vertical prediction #
# setwd("/Users/aesin/Desktop/Geo_analysis/HGT_position/For_circular/Per_species_vert//")
# pos_vert_df <- read_in_position_data("combined_rel_pos_assembled.tsv")


# ## Circularise the data ##
# circ_all <- circular(pos_all_df$V1*(2*pi))
# circ_const <- circular(pos_const_df$V1*(2*pi))
# circ_long <- circular(pos_long_df$V1*(2*pi))
# circ_vert <- circular(pos_vert_df$V1*(2*pi))

# ## Circular density calculation based on vonmises - much faster than wrapped ##
# bandwith <- 600

# circ_dens_all <- density.circular(circ_all, kernel = "vonmises", bw = bandwith)
# circ_dens_const <- density.circular(circ_const, kernal = "vonmises", bw = bandwith)
# circ_dens_long <- density.circular(circ_long, kernel = "vonmises", bw = bandwith)
# circ_dens_vert <- density.circular(circ_vert, kernel = "vonmises", bw = bandwith)


# ## Plot with differentially coloured polygons ##

# plot(circ_dens_all, points.plot = F, axes = F, control.circle = circle.control(type = "n"), nosort = T, lwd = 1.5, col = "white", ylim = c (-1.2, 1.2), plot.type = "circle", zero = (pi/2), rotation = "clock", main = "Gene density average across 12 Geobacillus species")

# axis.circular(at = NULL, labels = c("Origin", "", "Terminus", ""), rotation = "clock", zero = (pi /2), template = "none", tcl = 0.12, tcl.text = 0.2)

# all_line <- lines(circ_dens_all, lwd = 1.5, col = "black", zero = (pi/2), rotation = "clock")
# long_line <- lines(circ_dens_long, lwd= 0, col = "white", zero = (pi/2), rotation = "clock")
# #const_line <- lines(circ_dens_const, lwd= 1.5, col = "black", zero = (pi/2), rotation = "clock")
# #vert_line <- lines(circ_dens_vert, lwd = 0, col = "white", zero = (pi/2), rotation = "clock")

# lapply(polyclip(A=list("x"=long_line$x, "y"=long_line$y), B=list("x"=all_line$x, "y"=all_line$y), op="minus"), polygon, col = rgb(0,1,0,0.5), border = rgb(0,1,0,0.5))

# lapply(polyclip(B=list("x"=long_line$x, "y"=long_line$y), A=list("x"=all_line$x, "y"=all_line$y), op="minus"), polygon, col = rgb(1,0,0,0.5), border = rgb(1,0,0,0.5))

# vert_line <- lines(circ_dens_vert, lwd = 1.5, col = "blue", zero = (pi/2), rotation = "clock")

# legend("topleft", c("Long distance HGT enriched relative to all genes                              ", "Long distance HGT depleted relative to all genes", "Density of vertically inheritted genes"), lty = c(1,1,1), lwd = c(10,10,2), col = c(rgb(0,1,0,0.5), rgb(1,0,0,0.5), "Blue"))
























# #######################################
# ## For linear density representation ##
# #######################################

# # # Logspline #
# # logsp_dens_all <- logspline(pos_all_df)
# # logsp_dens_const <- logspline(pos_const_df)
# # logsp_dens_long<- logspline(pos_long_df)
# # logsp_dens_vert <- logspline(pos_vert_df)

# # # Normal density function #
# # norm_dns_all <- density(pos_all_df[,1])
# # norm_dns_const <- density(pos_const_df[,1])
# # norm_dns_long <- density(pos_long_df[,1])
# # norm_dns_vert <- density(pos_vert_df[,1])

# ### Plot logspline ###
# #plot.logspline(logsp_dens_all, xlim = c(0, 1.2), ylim = c(0, 2.5), lwd = 2, col = "black", xlab = "Genome position", ylab = "Density")
# #plot.logspline(logsp_dens_vert, add = TRUE, col = "blue", lwd = 2)
# #plot.logspline(logsp_dens_long, add = TRUE, col = "red", lwd = 2)

# #legend(1.0,2.0, c("All genes", "Vertical set", "HGT"), lty = c(1,1,1), lwd = c(2,2,2), col = c("Black", "Blue", "Red"))

# # plot.logspline(logsp_dens_const, add = TRUE, col = "red")

# ### Plot normal ###
# #plot(range(norm_dns_all$x, norm_dns_vert$x, norm_dns_long$x), range(norm_dns_all$y, norm_dns_vert$y, norm_dns_long$y), type = "n", xlab = "x", ylab = "Density")
# #lines(norm_dns_all, col = "black")
# #lines(norm_dns_vert, col = "blue")
# #lines(norm_dns_long, col = "red")

