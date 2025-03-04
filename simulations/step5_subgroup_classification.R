# divide highrisk subjects into subjects
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data_interaction/data')
library(data.table)

# read into PD-PRSs
dat.prs <- fread('../../pd_prs/data/simulated_PD_PRS_all.txt',data.table=F)
# read into simulated phenotypes
dat.phe <- fread('simu_pheno_withInter.txt',data.table=F)

dat <- merge(dat.phe,dat.prs,by.x='IID',by.y='eid')

### Identify subjects at high genetic risk for Trait 1
thresh=0.95
dat$high_risk <- ifelse(dat$simu_PRS_Observed_Beta_Trait1_all >
			  quantile(dat$simu_PRS_Observed_Beta_Trait1_all,
				   thresh), # !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
			1,0)
print(sum(dat$high_risk))
##### identify 4 subgroups
for(i in 1:4){
  col <- paste0('simu_PRS_Observed_Beta_Trait1_',i)
  prs <- unlist(subset(dat,select=col))
  cutoff <- quantile(prs,thresh)

  subgroup <- rep(NA,nrow(dat))
  subgroup[prs>cutoff & dat$high_risk==1] <- 1
  subgroup[prs<=cutoff & dat$high_risk==1] <- 0

  dat <- cbind.data.frame(dat,subgroup)
  colnames(dat)[ncol(dat)] <- paste0('subgroup',i)
}

write.table(dat,'simulated_subgroup.txt',
	    row.names=F,sep='\t',quote=F)

# Check mean values of phenotypes
prs <- paste0('simu_PRS_Beta_Trait1_',1:4)
phe <- paste0('Observed_Trait',1:4)

out <- c()
for(p in prs){
  for(pheno in phe){
    formu <- formula(paste0(pheno,'~',p))
    fit <- lm(formu,data=dat[dat$high_risk==1,])
    coef <- summary(fit)$coefficients
    add <- coef[2,]
    add <- c(p,pheno,add)

    out <- rbind(out,add)
  }
}
print(out)

