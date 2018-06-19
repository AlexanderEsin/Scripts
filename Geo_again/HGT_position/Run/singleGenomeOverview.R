invisible(sapply(HGTPos.all, source, .GlobalEnv))

require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("tidyverse", "wesanderson", "stringr")

# ------------------------------------------------------------------------------------- #`
# Read in data

message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData			<- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

# GI-specific data
GI_positions_data		<- readRDS(file.path(positionData_path, "AG_GI_positions.rds"))
byType_withGI_data		<- readRDS(file.path(positionData_path, "AG_perTypeGIData.rds"))

# tRNA data
crossSpecies_tRNA_data	<- readRDS(file.path(positionData_path, "AG_tRNA_data.rds"))

# Zone list
zoneBoundaryList	<- readRDS(file.path(positionData_path, "AG_zoneBoundaries.rds"))

message("\rReading in data... done\n")



# -------------------------------------- #
# Output path for figures
singleGenome_path	<- file.path(figureOutput_path, "singleGenomeOverview")
if(!dir.exists(singleGenome_path)) dir.create(singleGenome_path)


# ------------------------------------------------------------------------------------- #

AG_taxids	<- unique(perTypeData$All$allPosData$taxid)

genomeOverview_plotList	<- lapply(AG_taxids, function(speciesTaxid) {

	# Species name
	binomial		<- perTypeData$All$allPosData %>% filter(taxid == speciesTaxid) %>% summarise(binomial = unique(binomial)) %>% pull(binomial)
	message(binomial)
	binomialString	<- str_replace_all(binomial, " ", "_")

	# Process GI boundaries
	GI_boundaries	<- GI_positions_data$rawBoundary_data
	spec_GI_bounds	<- GI_boundaries %>% filter(Taxid == speciesTaxid)

	# Process the HGT and Vertical data
	spec_HGT_data	<- perTypeData$lHGT$`4`$allPosData %>% filter(taxid == speciesTaxid) %>% select(-c(CircStart, CircEnd)) %>% mutate(type = "HGT")
	spec_Ver_data	<- perTypeData$Ver$`3`$allPosData %>% filter(taxid == speciesTaxid) %>% select(-c(CircStart, CircEnd)) %>% mutate(type = "Ver")

	# Combine the HGT and Ver for plotting
	spec_combined	<- bind_rows(list(spec_HGT_data, spec_Ver_data)) %>% mutate(type = factor(type, levels = c("HGT", "Ver")))

	# Extract recent HGTs only
	spec_recent		<- spec_HGT_data %>% filter(Subgroup == TRUE)

	# Process the tRNA data
	spec_tRNA_data	<- crossSpecies_tRNA_data %>% filter(taxid == speciesTaxid)

	# Find the maximum density so we can scale the plot properly downstream
	densityFind		<- ggplot(data = spec_combined, aes(x = relGeneStart, color = type)) + stat_density(geom = "line", position = "identity", n = 2^12, adjust = 1/5)
	densityBuild	<- ggplot_build(densityFind)
	maxDensity		<- max(densityBuild$data[[1]][which(densityBuild$data[[1]]$scaled == 1), "y"])
	density_scale	<- maxDensity / 0.5

	# Plot
	speciesOverview_plot	<- ggplot(data = spec_combined, aes(x = relGeneStart, y = ..density.. / density_scale, color = type)) +
		# X axis
		scale_x_continuous(
			expand = expand_scale(mult = c(0.01, 0.01)),
			name = "Normalized genome position",
			limits = c(0, 1),
			breaks = seq(0, 1, by = 0.5),
			labels = c("Origin", "Terminus", "Origin"),
			minor_breaks = seq(0, 1, by = ((1 / 360) * 30)),
			sec.axis = dup_axis(
				name = NULL,
				breaks = rowMeans(zoneBoundaryList$fullRange[c("zoneMin", "zoneMax")]),
				labels = zoneBoundaryList$fullRange$zoneName)
		) +
		# Y axis
		scale_y_continuous(
			name = "Scaled Gene Density",
			# limits = c(0, 0.5),
			breaks = seq(0, 0.5, by = 0.1),
			minor_breaks = seq(0, 0.5, by = 0.05),
			labels = seq(0, 1, by = 0.2),
			sec.axis = dup_axis(
				name = NULL,
				breaks = c(0.55, 0.65, 0.75, 0.85),
				labels = c("All\nHGTs", "Recent\nHGTs", "Genomic\nIslands", "tRNA\nDensity"))
		) +
		geom_hline(yintercept = c(0, seq(0.5, 0.8, by = 0.1), 0.95), color = axisCol, alpha = 0.8) +
		# Plot the gene densities
		stat_density(
			geom = "line",
			position = "identity",
			n = 2^12,
			adjust = 1/5,
			size = 1) +
		# Line colours
		scale_color_manual(
			name = "Gene Type",
			values = c(dataTypeCols$HGT, dataTypeCols$Ver)
		) +
		# Plot HGT "ticks"
		geom_segment(
			data = spec_HGT_data,
			mapping = aes(x = relGeneStart, xend = relGeneStart, y = 0.51, yend = 0.59),
			color = dataTypeCols$HGT,
			size = 0.2
		) +
		# Plot recent HGT ticks
		geom_segment(
			data = spec_recent,
			mapping = aes(x = relGeneStart, xend = relGeneStart, y = 0.61, yend = 0.69),
			color = dataTypeCols$Recent,
			size = 0.3
		) +
		# Plot the GI boundaries data
		geom_rect(
			data = spec_GI_bounds,
			mapping = aes(xmin = GI_relStart, xmax = GI_relEnd, ymin = 0.71, ymax = 0.79),
			fill = "red",
			color = axisCol,
			inherit.aes = FALSE
		) +
		# Plot the tRNA data
		stat_density(
			data = spec_tRNA_data,
			mapping = aes(x = relStart, y = 0.81 + (..scaled.. / 10)),
			size = 1.5,
			geom = "line",
			n = 2^12,
			adjust = 1/25,
			color = axisCol,
			alpha = 0.9
		) +
		# Plot the zone boundaries
		geom_rect(
			data = zoneBoundaryList$fullRange,
			mapping = aes(xmin = zoneMin, xmax = zoneMax, ymin = -Inf, ymax = Inf),
			fill = zoneBoundaryList$fullRange$zoneCol_alpha,
			alpha = 0.15,
			inherit.aes = FALSE) +
		# Title and themes
		ggtitle(paste0(binomial, " overview")) +
		lightTheme +
		theme(
			axis.ticks = element_blank(),
			legend.background = element_rect(color = axisCol, fill = "white"),
			legend.position = "bottom"
		)

	quartz(width = 24, height = 12)
	print(speciesOverview_plot)
	quartz.save(file = file.path(singleGenome_path, paste0(binomialString, "_overview.pdf")), type = "pdf", dpi = 300)
	invisible(dev.off())

})

# ------------------------------------------------------------------------------------- #
