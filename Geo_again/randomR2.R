transZoneTotal_length	<- sum(zoneBoundaryList$fullRange$boundaryMax - zoneBoundaryList$fullRange$boundaryMin, na.rm = TRUE)


HTgenesInTransZones_list	<- lapply(1:nrow(zoneBoundaryList$fullRange), function(boundaryIndex) {
	
	boundaryRow	<- zoneBoundaryList$fullRange[boundaryIndex,]
	boundaryVal	<- boundaryRow$boundary
	
	transZoneMin	<- boundaryRow$boundaryMin
	transZoneMax	<- boundaryRow$boundaryMax
	
	if (is.na(boundaryVal)) return(NA)
	
	HGTgenesInTransZone <- subset(perTypeData$lHGT$'4'$allPosData, relGeneStart >= transZoneMin & relGeneEnd >= transZoneMin & relGeneStart <= transZoneMax & relGeneEnd <= transZoneMax & !taxid %in% outlierTaxid, select = -c(CircStart, CircEnd))
	
	HGTgenesInTransZone$transZoneID	<- boundaryRow$boundIndex
	
	return(HGTgenesInTransZone)
})

HTgenesInTransZones_list	<- HTgenesInTransZones_list[!is.na(HTgenesInTransZones_list)]
HTgenesInTransZones_df	<- bind_rows(HTgenesInTransZones_list)



x_tbl		<- dbSendQuery(conn, 'SELECT DISTINCT product FROM t1 WHERE protID = :protid_list')
dbBind(x_tbl, param = list(protid_list = unique(HTgenesInTransZones_df$protID)))
x_df		<- dbFetch(x_tbl)
dbClearResult(x_tbl)


sort(table(x_df$product), decreasing = T)


M genes 4% of all genes
M genes 8% of all genes