# Identify subgroups
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/analysis')
library(data.table)

# read data
#dat <- fread('simulated_phe_decomposed_prs.txt',data.table=F)
dat <- fread('simulated_phe_inter_decomposed_prs.txt',data.table=F)


### Identify subjects at high genetic risk for Trait 1
thresh=0.95
dat$high_risk <- ifelse(dat$PRS_B_disease >
                          quantile(dat$PRS_B_disease,
                                   thresh), # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                        1,0)
print(sum(dat$high_risk))
##### identify 4 subgroups
prs_name <- c('PRS_comp1','PRS_comp2','PRS_comp3','PRS_B_resid')
for(i in 1:4){
  prs <- unlist(subset(dat,select=prs_name[i]))
  cutoff <- quantile(prs,thresh)

  subgroup <- rep(NA,nrow(dat))
  subgroup[prs>cutoff & dat$high_risk==1] <- 1
  subgroup[prs<=cutoff & dat$high_risk==1] <- 0

  print(table(subgroup))
  dat <- cbind.data.frame(dat,subgroup)
  colnames(dat)[ncol(dat)] <- paste0('subgroup',i)
}

write.table(dat,'inter_simulated_subgroup.txt',
            row.names=F,sep='\t',quote=F)


