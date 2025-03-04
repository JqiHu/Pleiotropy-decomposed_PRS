# Compare relative change of traits across subgroups
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/analysis')
library(stringr)
library(ggplot2)
library(data.table)

#dat <- fread('simulated_subgroup.txt',data.table=F)
dat <- fread('inter_simulated_subgroup.txt',data.table=F)
dat.high <- dat[dat$high_risk==1,]
psprs <- c('PRS_comp1','PRS_comp2','PRS_comp3','PRS_B_resid')

## function to calculate relative change on one phenotype
## among high-risk individuals
getRelachan <- function(pheno){
  phe <- subset(dat.high,select=c('IID',pheno))
  colnames(phe)[2] <- 'trait'
  out <- c()
  for(i in 1:4){
    # input subgroup data
    iid <- dat.high$IID[dat.high[,paste0('subgroup',i)]==1]

    index <- ifelse(phe$IID %in% iid,1,0) # 1 represents target subgroup
    res <- (mean(phe$trait[index==1]) # proportion of level 1 or mean value in the certain PD-prs subgroup
              -mean(phe$trait[index==0]))/sd(phe$trait)
    out <- rbind(out,cbind(paste0('subgroup',i),res))
  }
  return(out)
}

# Run for each phe
phe <- paste0('Observed_Trait',2:4)
res <- c()
for(i in 1:length(phe)){
  add <- getRelachan(phe[i])
  add <- cbind(phe[i],add)
  res <- rbind(res,add)
}
res <- as.data.frame(res)
colnames(res) <- c('pheno','psPRS','change')
res$change <- as.numeric(res$change)
res$pheno <- str_remove_all(res$pheno,'Observed_')

# ----------------
# Visualization
## color for each psPRS
color <- c('tomato','skyblue','darkgreen','gold')
names(color) <- paste0('subgroup',1:4)
### label for each psPRS
label <- c(paste0('Component ',1:3,'-related Subgroup'),
           'Residual PRS Subgroup')
names(label) <- paste0('subgroup',1:4)

getPlot <- function(data){ # input data for specific phenotype group and group name
  p <-   ggplot(data,aes(x=`pheno`,y=`change`,fill=`psPRS`))+
    geom_bar(stat='identity',position='dodge',width = 0.8,alpha=0.7)+
    scale_fill_manual(values=color,labels=label)+
    #ggtitle(paste0(group,'-related traits'))+
    theme_bw()+
    xlab('')+
    ylab(paste0('Relative change'))+
    scale_y_continuous(
      limits = c(-0.25, 0.8),
      breaks = c(-0.25, 0, 0.25, 0.5, 0.75),
      labels = sprintf("%.2f",c(-0.25, 0, 0.25, 0.50, 0.75))
    ) +
    guides(fill=guide_legend(title="PD-Subgroup"))+
    theme(axis.title = element_text(size=12, face='bold'),
          axis.text.x = element_text(size=12, face = 'bold'),
          axis.text.y = element_text(size = 12,face = 'bold'),
	  legend.text = element_text(size = 14),
          legend.title = element_text(size = 16),
	  plot.margin = unit(c(1, 1, 1, 1), "cm"))+
    coord_flip()
  return(p)
}

p <- getPlot(res)

saveRDS(p,'inter_rela_chan.RDS')
pdf('inter_rela_chan_subgroup.pdf')
print(p)

dev.off()


