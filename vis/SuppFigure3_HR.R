# generate hr plot for 9 psPRSs
rm(list=ls())
library(data.table)

setwd('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/prs_analysis/')

hr <- as.data.frame(fread('hr_cad_updatedfollowup.csv')) # file for HRs

# input file for plot color for each path
col <- as.data.frame(fread('../visualization//path.col'))
color <- col$col
names(color) <- col$pathway

#hr_p <- apply(hr[,2:4],2,function(x){log10(x)})
hr[,2:4] <- apply(hr[,2:4],2,function(x){as.numeric(x)})
hr$pathway <- factor(hr$pathway,
	levels=hr$pathway)

library(ggplot2)
pdf('../visualization/hr_psprs.pdf',width=12)

ggplot()+
  geom_point(data=hr,aes(x=pathway, y=HR,col=pathway),size=1,show.legend=F,alpha=1)+
  xlab('')+
  ylim(1,NA)+
  guides(col=F)+
  geom_errorbar(data=hr,aes(x=`pathway`,ymin=`lower.95`,ymax=`upper.95`,col=`pathway`),width = 0.6,alpha=1)+
  theme_bw()+
  theme(axis.title = element_text(size=10,face = 'bold'),
                  axis.text.x = element_text(size = 12,angle = -45,hjust = 0,face = 'bold'),
                  axis.text.y = element_text(size = 10,face = 'plain'),
                  axis.line.y = element_line(colour=NULL), plot.margin = margin(10, 40, 10, 10)) +
  #ylim(NA,1.65)+
  scale_color_manual(aes(x=pathway),values=c(color,'black'))+
  geom_hline(yintercept=1,linetype=5,col="red")
  #geom_text(data=hr,aes(x=pathway, y=HR,col='black',label=label),size=2)

dev.off()
