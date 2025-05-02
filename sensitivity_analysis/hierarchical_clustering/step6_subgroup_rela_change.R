# Analyze associations between PD-PRSs and phenotypes of interest
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/clustering_sensitivity/hierachical_clustering/data')
library(data.table)
library(stringr)

### import number of clusters
args <- commandArgs(trailingOnly=T)
k <- args[1]
setwd(paste0('./cluster_',k,'/prs_analysis'))

# phenotype information
phe_data <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/pheno_all.tsv') # phenotype for all

### biological clusters for phenotypes
annot <- as.data.frame(c(rep('Demographic',5),rep('Obesity',2),rep('Respiratory',3),rep('Lifestyle',6),
                         rep('CVD',3),rep('Family history',4),rep('T2D',2),rep('Lipids',8),
                         rep('Protein',2),rep('Sex hormone',3),rep('Urinary system',5),rep('Metabolism',5),
                         rep('Liver',5),rep('Immune system',2))) ## biological annotations for phenotypes
rownames(annot) <- colnames(phe_data)[-1]
colnames(annot) <- 'annot'
phe_name <- colnames(phe_data)[-1]

# input high-risk individuals
high_risk <- readRDS('high_risk_ind.RDS')

phe_data <- phe_data[phe_data$eid %in% high_risk$eid,]

# Add fat percentages
update <- readRDS('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/ukbb_all_modifiable.RDS') # original
update <- as.data.frame(update)
update <- update[,-which(duplicated(colnames(update))==T)]
update <- merge(phe_data[,1:2],update,by='eid')

update$arm_fat_percentage <- (update$`23111`+update$`23115`)/2
update$leg_fat_percentage <- (update$`23119`+update$`23123`)/2
update_add <- subset(update,select=c('eid','23099','23127','arm_fat_percentage','leg_fat_percentage'))
colnames(update_add)[2:3] <- c('body_fat_percentage','trunk_fat_percentage')

phe_data <- merge(phe_data,update_add,by='eid',all.x=T)
phe_name <- colnames(phe_data)[-1]


## function to distinguish categorical traits from continuous ones
getIndex <- function(pheno){ # input phenotype name
  if(pheno %in% c('lifestyle_comb','alcohol','physical_activity','income',
                  'smoking','lifestyle_ideal_count','lifestyle_poor_count')){
    return('quali_multi')
  }else{
    res <- ifelse(nrow(unique(na.omit(subset(phe_data,select=pheno))))<3,'quali','quanti')
    return(unique(res))
  }
}

# name of pleiotropy-decomposed regions
files <- list.files('.','subgroup')
files <- files[grep('RDS',files)]
psprs <- str_replace_all(files,'subgroup_','')
psprs <- str_replace_all(psprs,'.RDS','')

## function to calculate relative change on one phenotype
## among top 5% CAD PRS individuals
getRelachan <- function(pheno){
  phe <- subset(phe_data,select=c('eid',pheno))
  colnames(phe)[2] <- 'trait'
  out <- c()
  for(prs in psprs){
    # input subgroup data
    temp <- readRDS(paste0('subgroup_',prs,'.RDS'))
    index <- ifelse(phe$eid %in% temp$eid,1,0) # 1 represents target subgroup
    if(getIndex(pheno)=='quali_multi'){ # only count for highest level
      sub <- length(which(na.omit(phe$trait[index==1])== # for categorical variable, only count for the max level
                                max(na.omit(phe$trait[index==1]))))/length(na.omit(phe$trait[index==1])) # proportion of max-level who are in the certain PD-prs subgroup
      rem <- length(which(na.omit(phe$trait[index==0])==
                                max(na.omit(phe$trait[index==0]))))/length(na.omit(phe$trait[index==0])) # proportion of max-level who are in the the rest PD-prs subgroups
      res <- (sub-rem)/sd(na.omit(phe$trait))
    } else{
      res <- (mean(na.omit(phe$trait[index==1])) # proportion of level 1 or mean value in the certain PD-prs subgroup
                -mean(na.omit(phe$trait[index==0])))/sd(na.omit(phe$trait)) #
    }
    out <- rbind(out,cbind(prs,res))
  }
  return(out)
}

# Phenotype of interest
pheno_name <- rownames(annot)[annot$annot %in% c('Obesity','Respiratory','T2D','Lipids','CVD')]

rela_chan <- data.frame()
for(phe in pheno_name){
  add <- cbind(phe,getRelachan(phe))
  rela_chan <- rbind(rela_chan,add)
}
colnames(rela_chan) <- c('pheno','psPRS','change')

## output relative change
write.table(rela_chan,
            'rela_change.tsv',
           row.names=F,quote=F,sep='\t')

##### Visualize
library(ggplot2)
library(colorspace)
library(Hmisc)
### colors
k <- as.numeric(k)
#colors <- qualitative_hcl(as.numeric(k), palette = "Pastel")
# Generate 15 colors using the HCL function with adjusted parameters:
# Here, l (luminance) is set to 60 (darker than typical pastels)
# and c (chroma) is set to 70 for increased saturation.
colors <- hcl(
  h = seq(0, 360, length.out = k+1)[-(k+1)],  # 15 evenly spaced hues
  c = 70,   # chroma: intensity of the color
  l = 60    # luminance: brightness of the color
)
names(colors) <- psprs

# plots
getPlot <- function(data,group){ # input data for specific phenotype group and group name
  p <- ggplot(data,aes(x=`pheno`,y=`change`,fill=`psPRS`))+
    geom_bar(stat='identity',position='dodge',width = 0.8,alpha=0.7)+
    scale_fill_manual(values=colors)+
    theme_bw()+
    xlab('Quantitative traits')+
    ylab(paste0('Relative change in ',group,'-related traits'))+
    guides(fill=guide_legend(title="PD-PRS"))+
    theme(axis.title = element_text(size=11),
          axis.text.x = element_text(size=11, face = 'plain'),
          axis.text.y = element_text(size = 11,face = 'plain'))+
    coord_flip()
  return(p)
}

annot$annot <- str_replace_all(annot$annot,'CVD','BP')
plot <- list()
for(group in c('BP','Lipids','Obesity',
               'Respiratory','T2D')){
  rela <- rela_chan[rela_chan$pheno %in% rownames(annot)[annot$annot==group],]

  ## formate name for triats
  rela$pheno <- str_replace_all(rela$pheno,'_',' ')
  rela$pheno <- capitalize(rela$pheno)
  
  rela <- as.data.frame(rela)
  rela$change <- as.numeric(rela$change)
  rela$psPRS <- factor(rela$psPRS,levels=psprs)
  
  add <- getPlot(rela,group)
  plot[[group]] <- add
}


# output plots
library(ggpubr)
pdf.options(reset = TRUE, onefile = FALSE)
pdf('rela_change.pdf',
    width=18,height=9)

print(
p <- ggarrange(plot[['BP']]+ggtitle('(A) BP-related traits'),
  plot[['Lipids']] + ggtitle('(B) Lipids-related traits'),
  plot[['Obesity']] + ggtitle('(C) Obesity-related traits'),
  plot[['Respiratory']] + ggtitle('(D) Smoking traits'),
  plot[['T2D']] + ggtitle('(E) T2D-related traits'),
  ncol=3,nrow=2,common.legend=T,
  legend='right',align='hv')
)

dev.off()

saveRDS(p,'plot_rela_change.RDS')

