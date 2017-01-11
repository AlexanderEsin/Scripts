library(logspline)

setwd("/Users/aesin/Desktop/Deer/Selection_analysis/BGM/BGM_100anc_no_arma/")
posterior_edge_tbl <- read.table(file = "Edge_support.csv", sep = ",", header =T)

names(posterior_edge_tbl)[5] <- "Site1_2_total"
logsp_ppedges <- logspline(posterior_edge_tbl$Site1_2_total)

plot.logspline(logsp_ppedges, xlim = c(0,1), what = "d")
abline(v = 0.655125783, col = "red", lwd =2)