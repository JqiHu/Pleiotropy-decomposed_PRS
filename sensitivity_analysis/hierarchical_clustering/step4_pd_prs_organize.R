# Organize PD-PRSs into one table
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/clustering_sensitivity/hierachical_clustering/data')
library(data.table)
library(stringr)

### import number of clusters
args <- commandArgs(trailingOnly=T)
k <- args[1]
setwd(paste0('./cluster_',k,'/prs_calculation'))

dat <- c()
for(i in 1:k){
  if(k==15 & i==13){
    next
  }
  tmp <- fread(paste0('Cluster',i,'.profile'),
	       select=c('IID','SCORE'))
  if(nrow(tmp)==0){
    next 
  }else{
    colnames(tmp) <- c('eid',paste0('PD_PRS_cluster',i))
    if(i==1){
      dat <- tmp
    }else{
      dat <- merge(dat,tmp,by='eid')
    }
  }
}

write.table(dat,'../PD_PRS.txt',row.names=F,sep='\t',quote=F)
