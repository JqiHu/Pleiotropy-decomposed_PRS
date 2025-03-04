# generate correlation heatmap between 9 psPRSs
setwd('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/prs_analysis/')
library('corrplot')

# input correlation coefficient matrix and p value matrix
corr <- read.csv('corr_prs_coef.csv')
corr_p <- read.csv('corr_prs_p.csv')
rownames(corr) <- corr$X
corr <- corr[,-1]
rownames(corr_p) <- corr_p$X
corr_p <- corr_p[,-1]

# plot heatmap with upper as color and lower as numeric for correlations between ps-PRSs
col=colorRampPalette(c("navy", "white"))
col2=colorRampPalette(c("white", "firebrick3"))

diag(corr) <- 0 # make diagnal white
corr <- as.matrix(corr)
corr_p <- as.matrix(corr_p)

pdf('../visualization/corr_psprs.pdf',
    height=12,width=12)

corrplot(corr,type='upper',method='color',
         is.corr=F,tl.col='black',tl.cex=0.8,
         col=c(col2(201)),col.lim = c(-0.01,0.01),
         cl.length=3,cl.offset = seq(-0.1,0.1,length.out=3),
         p.mat=corr_p,sig.level=0.05/(9*9),
         insig='label_sig',pch='*',
         pch.col='gold',tl.pos='d',
         addgrid.col='grey')
corrplot(corr_p,add=T,type='lower',number.cex=0.8,
         method='number',diag=FALSE,tl.pos="n", tl.cex=0.2,
         cl.pos="n",col='black',number.digits=4)

dev.off()
