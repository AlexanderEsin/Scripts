# /// Per-Species position analysis - superceded by the per-Branch analysis ///



# penalty		<- "4"
# bandwith	<- 2000
# binomials	<- unique(perTypeData$All$allPosData$binomial)

# perSpeciesCircStart_list	<- lapply(binomials, function(binomial) {
# 	hgtRelStart		<- perTypeData$lHGT[[penalty]]$allPosData$relGeneStart[which(perTypeData$lHGT[[penalty]]$allPosData$binomial == binomial)]
# 	hgtCircStart	<- circular(hgtRelStart * (2 * pi), zero = pi / 2, rotation = "clock", modulo = "2pi")
# 	return(hgtCircStart)
# })
# names(perSpeciesCircStart_list)	<- binomials



# perSpeciesCircSummary_list <- lapply(1:length(perSpeciesCircStart_list), function(index) {
# 	speciesCircSum	<- summary.circular_AE(perSpeciesCircStart_list[[index]])
# 	speciesName		<- names(perSpeciesCircStart_list)[index]

# 	SpeciesCircSum_df	<- cbind(Species = speciesName, as.data.frame(t(speciesCircSum)), stringsAsFactors = FALSE)
# 	return(SpeciesCircSum_df)
# })
# perSpeciesCircSummary_df	<- bind_rows(perSpeciesCircSummary_list)


# quartz(width = 8, height = 8)
# perSpeciesAvHGT_plot	<- ggplot(perSpeciesCircSummary_df, aes(x = Mean, y = Rho, color = Species, label = Species)) +
# 	guides(color = FALSE) +
# 	scale_x_continuous(labels = c("Origin", "Terminus"), breaks = c(0, pi), limits = c(0, 2 * pi)) +
# 	scale_y_continuous(limits = c(0, 0.25)) +
# 	coord_polar() +
# 	geom_point() +
# 	geom_text_repel(aes(label = ifelse(Rho > 0.1, as.character(Species),'')), hjust = -0.05, vjust = 0.5, angle = -25) +
# 	darkTheme
# perSpeciesAvHGT_plot



# perSpeciesDens_data	<- lapply(binomials, function(binomial) {

# 	# print(binomial)
# 	allCircStart	<- perTypeData$All$allPosData$CircStart[which(perTypeData$All$allPosData$binomial == binomial)]
# 	hgtCircStart	<- perTypeData$lHGT[[penalty]]$allPosData$CircStart[which(perTypeData$lHGT[[penalty]]$allPosData$binomial == binomial)]
# 	vertCircStart	<- perTypeData$Ver[[penalty]]$allPosData$CircStart[which(perTypeData$Ver[[penalty]]$allPosData$binomial == binomial)]
	
# 	# # Calculate circular density (on the start position)
# 	allDensity	<- density.circular(allCircStart, kernel = "vonmises", bw = bandwith)
# 	HGTDensity	<- density.circular(hgtCircStart, kernel = "vonmises", bw = bandwith)
# 	vertDensity	<- density.circular(vertCircStart, kernel = "vonmises", bw = bandwith)

# 	# Number of genes for each density plot
# 	numGenes	<- length(HGTDensity$data)

# 	return(list(dataDensityA = HGTDensity, dataDensityB = vertDensity, bgDensity = allDensity, numGenes = numGenes))
# })
# names(perSpeciesDens_data)	<- binomials

# ## Interesting species from the perSpeciesAvHGT_plot
# binomialsToPlot	<- list("Parageobacillus thermoglucosidasius", "Geobacillus sp. Y4.1MC1", "Geobacillus sp. Y412MC61", "Geobacillus subterraneus")

# quartz(width = 20, height = 6.5)
# par(mfrow = c(1, 4))
# par(mar = c(0, 0, 0, 0))
# invisible(lapply(binomialsToPlot, function(binomial) {
# 	print(binomial)
# 	basicPositionPlot(
# 		dataDensityA = perSpeciesDens_data[[binomial]]$dataDensityA,
# 		dataDensityB = perSpeciesDens_data[[binomial]]$dataDensityB,
# 		bgDensity = perSpeciesDens_data[[binomial]]$bgDensity,
# 		uin = 1.8,
# 		shrink = 1.1,
# 		tcl.offset = 0.8,
# 		titleCex = 1.5,
# 		titleName = paste0(binomial, "\nlHGT Enrichment at Penalty = ", penalty, "\nGenes = ", perSpeciesDens_data[[binomial]]$numGenes)
# 	)
# }))



