# organize PD-PRSs calcualted
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/pd_prs/data/pd_prs')
library(data.table)
library(stringr)

files <- list.files('.','profile')
prefix <- str_remove_all(files,'.profile')

dat <- c()
for(i in 1:length(prefix)){
  tmp <- fread(files[i],data.table=F,select=c('IID','SCORE'))
  colnames(tmp) <- c('eid',prefix[i])
  if(i==1){
    dat <- tmp
  }else{
    dat <- merge(dat,tmp,by='eid')
  }
}

# standardize PRS
dat <- as.data.frame(dat)
dat[,-1] <- apply(dat[,-1],2,scale)

write.table(dat,'../simulated_PD_PRS_all.txt',row.names=F,sep='\t',quote=F)
