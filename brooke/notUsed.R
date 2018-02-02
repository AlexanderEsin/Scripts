# ### Identifying common motifs ###
# P4Cat 		<- c(rep(1, length(finalUniqSeq_data_list$P4$UniqSeqData)), rep(0, length(randomDNA_bg)))
# P4_motifs	<- findMotif(all.seq = append(finalUniqSeq_data_list$P4$UniqSeqData, randomDNA_bg), category = P4Cat, start.width = 9, mc.cores = 16, both.strand = FALSE, enriched.only = TRUE, max.motif = 10)

# motifPosition_list	<- lapply(names(P4_motifs$motifs), function(motif) {
# 	motifPositions	<- P4_motifs$motifs[[motif]]@match$pos
# 	motifPos_df		<- data.frame(MotifSeq = rep(motif, length(motifPositions)), Position = motifPositions, stringsAsFactors = FALSE)
# 	return(motifPos_df)
# })
# motifPosition_df	<- bind_rows(motifPosition_list)

# ggplot(data = motifPosition_df, aes(x = MotifSeq, y = Position, color = MotifSeq)) + geom_violin()

