#!/usr/bin/env Rscript

# Load master variables and HGT position functions
invisible(sapply(HGTPos.all, source, .GlobalEnv))

# Any libraries explicitly used in the script
require(pacman, warn.conflicts = FALSE, quietly = TRUE)
p_load("tidyverse", "RCircos")


# Replacing native ideogram plot so that we can shade the chromosome band without the need for actual band data #
RCircos.Draw.Chromosome.Ideogram <- function(ideo.pos = NULL, ideo.width = NULL) {
    RCircos.Cyto <- RCircos.Get.Plot.Ideogram()
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    if (is.null(ideo.pos)) 
       ideo.pos <- RCircos.Par$chr.ideo.pos
    if (is.null(ideo.width)) 
       ideo.width <- RCircos.Par$chrom.width
    outerPos <- ideo.pos + ideo.width
    innerPos <- ideo.pos
    chromosomes <- unique(RCircos.Cyto$Chromosome)
    RCircos.Track.Outline(outerPos, innerPos, num.layers = 1, 
        chromosomes, track.colors = rep("white", length(chromosomes)))
    whiteBands <- which(RCircos.Cyto$BandColor == "white")
    # / Changed below: code no longer dependent on there being white bands present for colored band to be colored
    if (!length(whiteBands)) {
        darkBands <- RCircos.Cyto
    } else {
        darkBands <- RCircos.Cyto[-whiteBands, ]
    }

    for (aBand in seq_len(nrow(darkBands))) {
        aColor <- darkBands$BandColor[aBand]
        aStart <- darkBands$StartPoint[aBand]
        aEnd <- darkBands$EndPoint[aBand]
        posX <- c(RCircos.Pos[aStart:aEnd, 1] * outerPos, RCircos.Pos[aEnd:aStart, 
            1] * innerPos)
        posY <- c(RCircos.Pos[aStart:aEnd, 2] * outerPos, RCircos.Pos[aEnd:aStart, 
            2] * innerPos)
        polygon(posX, posY, col = aColor, border = NA)
    }
}
assignInNamespace("RCircos.Draw.Chromosome.Ideogram", RCircos.Draw.Chromosome.Ideogram, "RCircos")

# Replacing native link plot so non-integer linewidth can be used #
RCircos.Link.Plot <- function(link.data = NULL, track.num = NULL, by.chromosome = FALSE, start.pos = NULL, genomic.columns = 3, is.sorted = TRUE, lineWidth = rep(1, nrow(link.data))) {
    if (is.null(link.data)) 
        stop("Link data missing in RCircos.Link.Plot().\n")
    if (by.chromosome != TRUE && by.chromosome != FALSE) 
        stop("Error: by.chromosome must be either TRUE or FALSE.\n")
    ## No reason to have linewidth limited to integers! ##
    # if (length(which(lineWidth < 1)) > 1) 
    #     stop("Line width must be positive integer.")
    RCircos.Par <- RCircos.Get.Plot.Parameters()
    if (is.null(start.pos)) {
        locations <- RCircos.Get.Plot.Boundary(track.num, side = "in", 
            inside.pos = NULL, outside.pos = NULL, FALSE)
        line.start <- locations[1]
    }
    else {
        if (start.pos >= 1) 
            stop("Link line must be inside chromosome ideogram")
        line.start <- RCircos.Par$chr.ideo.pos * start.pos
    }
    if (is.null(genomic.columns) || genomic.columns < 3) 
        stop("Incorrect number of columns for genomic position.\n")
    link.data <- RCircos.Get.Paired.Points.Positions(link.data, 
        genomic.columns, plot.type = "link")
    link.colors <- RCircos.Get.Link.Colors(link.data, genomic.columns, 
        by.chromosome)
    RCircos.Pos <- RCircos.Get.Plot.Positions()
    base.positions <- RCircos.Pos[, 1:2] * line.start
    for (a.link in seq_len(nrow(link.data))) {
        point.one <- as.numeric(link.data$LinkStart[a.link])
        point.two <- as.numeric(link.data$LinkEnd[a.link])
        if (point.one > point.two) {
            point.one <- link.data$LinkEnd[a.link]
            point.two <- link.data$LinkStart[a.link]
        }
        P0 <- as.numeric(base.positions[point.one, ])
        P2 <- as.numeric(base.positions[point.two, ])
        links <- RCircos.Link.Line(P0, P2)
        lines(links$pos.x, links$pos.y, type = "l", lwd = lineWidth[a.link], 
            col = link.colors[a.link])
    }
}
assignInNamespace("RCircos.Link.Plot", RCircos.Link.Plot, "RCircos")


# Make link dataframe for input to circos linkplot. Input is a list data structure output from circos_pairwise_data_list function #
circos_pairwise_link <- function(pairwise_data_list) {
    pairwise_df_list <- lapply(pairwise_data_list, function(group) {
        return(group$Merged_pairwise)
    })
    return(do.call(rbind.data.frame, c(pairwise_df_list, make.row.names = FALSE)))
}

# Plot linked gene positions on bacterial chromosome #
circos_link_plot <- function(link_data, chrom_pad = NULL, title = NULL, linewidth = 0.2, title_cex = 1) {
    bact_ideogram <- data.frame(Chromosome = 1, ChromStart = 0, ChromEnd = 3000000, band = 1, stain = "acen")
    RCircos.Set.Core.Components(cyto.info = bact_ideogram, chr.exclude = NULL, tracks.inside = 1, tracks.outside = 0)

    rcircos.params  <- RCircos.Get.Plot.Parameters()
    rcircos.params$base.per.unit    <- 300
    RCircos.Reset.Plot.Parameters(rcircos.params)
    if (!is.null(chrom_pad)) {
        rcircos.params$chrom.paddings   <- 300
        RCircos.Reset.Plot.Parameters(rcircos.params)
    }
    rcircos.cyto    <- RCircos.Get.Plot.Ideogram()
    rcircos.cyto$BandColor  <- "grey80"
    RCircos.Reset.Plot.Ideogram(rcircos.cyto)

    # RCircos.Set.Plot.Area()
    plot.new()
    par(mai = rep(0, 4))
    plot.window(c(-1.5, 1.5), c(-1.5, 1.5))

    # RCircos.Set.Plot.Area(margins = 0)
    # par(omi = rep(0, 4))
    # warning(rcircos.params$plot.radius)
    RCircos.Draw.Chromosome.Ideogram()
    RCircos.Highligh.Chromosome.Ideogram()
    RCircos.Link.Plot(link_data, track.num = 1, is.sorted = FALSE, lineWidth = rep(linewidth, nrow(link_data)))
    if (!is.null(title)) {
        title(main = title, line = -1, cex.main = title_cex)
    }
}


normalise_to_3M <- function(number) {return(floor(number * 3000000))}


# ------------------------------------------------------------------------------------- #
message("\nReading in data...", appendLF = FALSE)

# General position data
perTypeData         <- readRDS(file.path(positionData_path, "AG_perTypeData.rds"))

message("\rReading in data... done\n")


# ------------------------------------------------------------------------------------- #


geneLink_path    <- file.path(figureOutput_path, "Gene_link")

if (!dir.exists(geneLink_path)) dir.create(geneLink_path, recursive = TRUE)


# ------------------------------------------------------------------------------------- #


x <- mclapply(unique(perTypeData$All$allPosData$orthGroup), function(group) {

    group_data  <- perTypeData$All$allPosData %>% filter(orthGroup == group)

    # Returm if only one gene in group
    if (nrow(group_data) == 1) {
        return(NA)
    }

    # Check that they are 1-to-1 orthologs
    if (length(group_data$binomial) != 23 || length(unique(group_data$binomial)) != 23) {
        return(NA)
    }

    # Get the first species - keystone
    keystone_data   <- group_data[1, ]
    keystone_out    <- data.frame(
        binomial = keystone_data$binomial,
        START = normalise_to_3M(min(c(keystone_data$relGeneStart, keystone_data$relGeneEnd))),
        END = normalise_to_3M(max(c(keystone_data$relGeneStart, keystone_data$relGeneEnd))),
        stringAsFactors = FALSE)

    # Process the rest of the species
    remain_data     <- group_data[-1, ]

    # Fix the dnaA gene positions
    remain_data     %<>% mutate(relGeneStart = case_when(
        relGeneStart == 1 ~ 0,
        TRUE ~ relGeneStart))

    # Linked genes output
    remain_out      <- data.frame(
        Chromosome = rep(1, nrow(remain_data)), 
        chromSTART = rep(keystone_out$START, nrow(remain_data)),
        chromEND = rep(keystone_out$END, nrow(remain_data)),
        Chromosome.1 = rep(1, nrow(remain_data)),
        START.1 = normalise_to_3M(remain_data$relGeneStart),
        END.1 = normalise_to_3M(remain_data$relGeneEnd),
        stringsAsFactors = FALSE
    )

    remain_correct  <- remain_out %>% 
        mutate(chromSTART.1 = case_when(
            START.1 > END.1 ~ END.1,
            TRUE ~ START.1)) %>%
        mutate(chromEND.1 = case_when(
            END.1 < START.1 ~ START.1,
            TRUE ~ END.1)) %>% 
        select(-c(START.1, END.1))

    # Check for distance
    remain_distFilt <- remain_correct %>% filter(abs(chromSTART.1 - keystone_out$START) > 150000)

    return(list(keystone = keystone_out, Merged_pairwise = remain_distFilt))
}, mc.cores = 20)

x <- x[!is.na(x)]

y <- circos_pairwise_link(x)

devSVG(file = file.path(geneLink_path, "1to1_orth_linkPlot_min5percent_dist.svg"), width = 12, height = 12, bg = "white")
circos_link_plot(y)
invisible(dev.off())



# Calculate the delta distance to ori between ortholog pairs
z <- y %>% mutate(startToOri = case_when(
        (chromSTART / 3000000) > 0.5 ~ 1 - (chromSTART / 3000000),
        TRUE ~ chromSTART / 3000000)) %>%
    mutate(start1ToOri = case_when(
        (chromSTART.1 / 3000000) > 0.5 ~ 1 - (chromSTART.1 / 3000000),
         TRUE ~ chromSTART.1 / 3000000)) %>%
    mutate(deltaOriDist = abs(startToOri - start1ToOri))


distOriDist_plot <- ggplot(z, aes(x = deltaOriDist)) +
    scale_x_continuous(limits = c(0, 0.5)) +
    geom_histogram(bins = 100, color = "black", fill = "grey70") +
    lightTheme +
    theme(
        panel.grid.minor.y = element_blank()
    )

quartz(width = 12, heigh = 6)
print(distOriDist_plot)
quartz.save(file = file.path(geneLink_path, "deltaDistOri_histogram.pdf"), type = "pdf")
invisible(dev.off())














