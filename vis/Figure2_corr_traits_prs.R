# generate correlation heatmap between 9 psPRSs and phenotypes
setwd('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/prs_analysis/')
library(stringr)
library('pheatmap')
library(data.table)
library(Hmisc)

# input correlation matrix
corr <- as.data.frame(fread('corr_prs_phe_coef.tsv'))
corr_p <-  as.data.frame(fread('corr_prs_phe_p.tsv'))
phe <- corr$V1 # all phenotype names
corr <- corr[,-1]
corr_p <- corr_p[,-1]

# capitalize phenotype names and remve '_'
phe <- capitalize(phe)
phe <- str_replace_all(phe,'_',' ')
phe[c(3,6)] <- toupper(phe[c(3,6)]) # BMI and TDI
rownames(corr) <- phe
rownames(corr_p) <- phe

### For updated version
### need to adjust the sequence of traits
phe_ordered <- c(phe[1:7],phe[56:59],phe[8:55])
corr2 <- corr[match(phe_ordered,rownames(corr)),]
corr_p2 <- corr_p[match(phe_ordered,rownames(corr_p)),]

# phenotype cluster annotations
annot <- c(rep('Demographic',5),rep('Obesity',6),rep('Respiratory',3),rep('Lifestyle',6),
           rep('CVD',3),rep('Family history',4),rep('T2D',2),rep('Lipids',8),
           rep('Protein',2),rep('Sex hormone',3),rep('Urinary system',5),rep('Metabolism',5),
           rep('Liver',5),rep('Immune system',2)) ## biological annotations for phenotypes
Category <- as.data.frame(annot,row.names=phe_ordered)
colnames(Category) <- 'Category'
Category$Category <- factor(Category$Category,
	levels=unique(annot))



# plot heatmap for correlation coeffcients between ps-PRSs and phenotypes
library(grid)
library(gtable)

pdf('../visualization/corr_prs_phe_update.pdf',
	height=9,width=10)

pheatmap(corr2,
	 color = colorRampPalette(c("navy", "white", "firebrick3"))(600),
	 clustering_method = "average", cluster_rows = F, cluster_cols = F, 
	 breaks = seq(-0.06,0.06,0.0002),legend_breaks = seq(-0.06,0.06,0.02),
	 fontsize_row=10,fontsize_col=15,legend_labels = c(seq(-0.06,0.06,0.02)),
	 display_numbers = ifelse(corr_p2 < 0.05/(59*9), "*", ""),
   annotation_row=Category,angle_col=315,
	 cellwidth = 36, cellheight = 9)
# Add custom legend title
grid.text("Beta coefficients", x = 0.74, y = 0.88, rot = 90, gp = gpar(fontsize = 12))

dev.off()


# check the minimal p-value for each PD-PRS
mini <- apply(corr_p2, 2, function(x){rownames(corr_p2)[which.min(x)]})




