# evaluate prediction preformance of PD-PRSs for CAD

library('data.table')
library('stringr')
library('survival')

setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data') # path to phenotype and PD-PRS files

# PRS for CAD
prs_cad <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS/cad_prs.tsv')
prs_cad <- prs_cad[,c(2,6)] # subtract iid and score
colnames(prs_cad)[2] <- 'CAD_prs'
prs_cad$CAD_prs <- scale(prs_cad$CAD_prs,
                center=T, scale=T) # scale

ps_prs <- fread('./psPRS_V2_absolute/ps_prs.tsv') # input PD-prs files
data <- merge(prs_cad,ps_prs,by='IID')

pheno_cad <- fread('./phenotype/pheno_cad.tsv') # input CAD 
time_cad <- fread('./phenotype/ukbb_CAD.tsv',select=c('eid','age_end')) # update follow-up time
data <- merge(pheno_cad,data,
	by.x='eid',by.y='IID',
	all.y=T) # merge CAD and psPRS
data <- merge(data,time_cad,
	      by='eid',all.x=T)
## covariates
cov <- fread('./phenotype/pheno_all.tsv',select=c('eid','sex'))
pc <- fread('/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_pheno/MR_time/adjust/adjust.csv') # principle components
pc <- subset(pc,select=c('eid',
                  paste0('PC',1:4)))
cov <- merge(cov,pc,by='eid')

data <- merge(data,cov,by='eid',all.x=T)

ps_names <- colnames(ps_prs)[-1] # name for pathways
ps_names <- c(ps_names,'CAD_prs') # add the overall CAD PRS

# separate models to predict risk for CAD
getHR <- function(prs){ # select one prs one time
  sur <- Surv(data$age_end,data$CAD)
  res_cox <- summary(coxph(sur ~ unlist(subset(data,select=prs))+
			   data$age_recruit+data$sex+
			   data$PC1+data$PC2+data$PC3+data$PC4))
  res <- cbind.data.frame(prs,
	   t(res_cox$conf.int[1,]),res_cox$coefficients[1,5])
  return(res)
}

hr <- data.frame()
for(i in 1:length(ps_names)){
  hr <- rbind.data.frame(hr,getHR(ps_names[i]))
}
hr <- hr[,-3] # remove exp(-beta)
colnames(hr) <- c('pathway','HR','lower.95','upper.95','p_val')
hr <- hr[order(hr$`HR`,decreasing=T),] # rank by HR
rownames(hr) <- NULL

write.csv(hr,
	  './prs_analysis_V2_absolute/hr_cad_updatedfollowup.csv',
	  row.names=F)
