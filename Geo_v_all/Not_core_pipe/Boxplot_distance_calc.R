library(reshape2)
library(ggplot2)

setwd("/users/aesin/desktop/Geo_v_all/-100")
gdata<-read.delim(dir()[18],header=T,fill=T)


data <- melt(gdata,id.vars="Group_no",measure.vars=c(2:ncol(gdata)),na.rm=T)

p<-ggplot(data)+geom_boxplot(aes(x=Group_no,y=value))

p2<-ggplot(data)+geom_boxplot(aes(x=Group_no,y=value),color="forestgreen")+geom_point(data=data[1:nrow(gdata),], aes(x=Group_no,y=value),colour="red",size=2)