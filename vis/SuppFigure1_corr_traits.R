# generate correlation plot for 43 CAD-correlated traits
#source('/gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/code/gnova_43/gnova_43_res.R') # path to gnova_43_res.R
setwd('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/visualization/')
load('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/visualization/data_for_supF1.RData')
library(pheatmap)
library(stringr)
clu <- info$Cluster # pathways for traits

annot <- as.data.frame(clu)
rownames(annot) <- trait # annotation bar in plot
gap <- rep(0,length(unique(clu)))
for(i in 1:length(unique(clu))){gap[i] <- sum(table(clu)[1:i])} # separate traits in different pathways

corr_out[corr_out>1] <- 1
corr_out[corr_out<-1] <- -1

phe <- colnames(corr_out)

### modify names
phe <- str_replace_all(phe,'longetivity','longevity')
phe <- str_replace_all(phe,'Vulvular','Valvular')

colnames(corr_out) <- phe
rownames(corr_out) <- phe

# heatmap plot
pdf('fa1_1.pdf',width=15,height=15)
pheatmap(corr_out,border=TRUE,border_color='grey60',
	 color = colorRampPalette(c("navy",'white',"firebrick3"))(500),
	 clustering_method = "average",cluster_rows=F,cluster_cols=F,
         breaks = seq(-1,1,2/500),legend_breaks = seq(-1,1,2/5),legend_labels = c(seq(-1,1,2/5)),
	 display_numbers = ifelse(corr_p < 0.05/(nrow(corr_p)^2), "*", ""),
         annotation_row=annot,gaps_row=gap,cellwidth=10,cellheight=10,gaps_col=gap)
dev.off()

pdf('fa1_1_p.pdf',width=15,height=15)
pheatmap(corr_p,border=TRUE,border_color='grey60',
	 color = colorRampPalette(c("firebrick3",'white',"navy"))(5000),
	 clustering_method = "average",cluster_rows=F,cluster_cols=F,
         breaks = seq(0,0.002,0.0000004),legend_breaks = seq(0,0.002,0.0004),legend_labels = c(seq(0,0.002,0.0004)),
	 display_numbers = ifelse(corr_p < 0.05/(nrow(corr_p)^2), "*", ""),
         annotation_row=annot,gaps_row=gap,cellwidth=10,cellheight=10,gaps_col=gap)
dev.off()

