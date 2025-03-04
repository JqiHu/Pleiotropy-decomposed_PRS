# Organize interaction p-values
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/interaction_new_pipeline/data')
library(data.table)
library(stringr)

# Function to organize interaction results
getBHres <- function(trans,covar){
  tmp <- fread(paste0('interaction_',trans,'_adj_',covar,'.csv'),data.table=F)
  # clean trait names
  tmp$trait <- str_remove_all(tmp$trait,'trans_')
  ## Extract the last component
  tmp$trait <- sub("_([^_]+)$", "", tmp$trait)
  ##### format trait names
  tmp$trait[tmp$trait=='fid_1160'] <- 'sleep_duration'
  tmp$trait[tmp$trait=='fid_1269'] <- 'exposure_to_tobacco_smoke_at_home'
  tmp$trait[tmp$trait=='fid_1279'] <- 'exposure_to_tobacco_smoke_outside_home'
  tmp$trait[tmp$trait=='fid_189'] <- 'TDI'
  tmp$trait[tmp$trait=='fid_23099'] <- 'body_fat_percentage'
  tmp$trait[tmp$trait=='fid_23127'] <- 'trunk_fat_percentage'
  tmp$trait[tmp$trait=='fid_3062'] <- 'FVC'
  tmp$trait[tmp$trait=='fid_3063'] <- 'FEV1'
  tmp$trait[tmp$trait=='fid_3064'] <- 'PEF'
  tmp$trait[tmp$trait=='fid_4079'] <- 'DBP'
  tmp$trait[tmp$trait=='fid_4080'] <- 'SBP'

  ### adjusted p-values by BH
  tmp$pval.adj <- p.adjust(abs(tmp$`p-value`),method='BH')
  tmp$pval.adj <- ifelse(tmp$`p-value`<0,tmp$pval.adj*(-1),tmp$pval.adj)

  ##### significant pairs
  tmp$sig <- 0
  tmp$sig[abs(tmp$pval.adj)<0.05] <- 1
  tmp <- tmp[order(abs(tmp$pval.adj),decreasing=F),]

  return(tmp)
}

# Run for all results
###### Baseline covariates
out <- getBHres('INT','baseline')
write.csv(out,'interaction_INT_baseline_BH.csv',row.names=F)

###### Residuals
out <- getBHres('Residual','all_all')
write.csv(out,'interaction_Residual_adj_all_all_BH.csv',row.names=F)
