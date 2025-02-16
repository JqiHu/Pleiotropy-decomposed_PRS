# Organize PRS and phenotype data
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/analysis')
library(data.table)
library(stringr)

# Read PRSs
files <- list.files('../prs/','profile')
prefix <- str_remove_all(files,'simu_PRS_')
prefix <- str_remove_all(prefix,'.profile')
for(i in 1:length(files)){
  tmp <- fread(paste0('../prs/',files[i]),
	       data.table=F,select=c('IID','SCORESUM'))
  colnames(tmp)[2] <- paste0('PRS_',prefix[i])
  if(i==1){
    dat.prs <- tmp
  }else{
    dat.prs <- merge(dat.prs,tmp,by='IID')
  }
}
## standardize PRSs
dat.prs <- as.data.frame(dat.prs)
dat.prs[,-1] <- apply(dat.prs[,-1],2,scale)

# phenotypes
#dat.phe <- fread('../../../simu_data/data/simu_pheno/simulated_phenotype.txt')
dat.phe <- fread('../../../simu_data_interaction/data/simu_pheno_withInter.txt')

dat <- merge(dat.prs,dat.phe,by='IID')

write.table(dat,'simulated_phe_inter_decomposed_prs.txt',
	    row.names=F,sep='\t',quote=F)
