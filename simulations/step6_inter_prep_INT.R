# Transform simulated phenotype by adjusting PD-PRSs
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data_interaction/data')
library(data.table)

# read PRS and phenotypes
dat <- fread('simulated_subgroup.txt',data.table=F)

# Function to adjust specific phenotype
getAdj <- function(phe,prs){
  formu <- formula(paste0(phe,'~',prs))
  fit <- lm(formu,data=dat)
  residuals <- fit$residuals
  trans_resid <- qnorm((rank(residuals,na.last="keep")-0.5)/sum(!is.na(residuals)))

  return(trans_resid)
}

# Transform Trait2 to Trait4
for(i in 2:4){
  phe <- paste0('Observed_Trait',i)
  prs <- paste0('simu_PRS_Observed_Beta_Trait1_',i-1)

  add <- getAdj(phe,prs)
  dat <- cbind.data.frame(dat,add)
  colnames(dat)[ncol(dat)] <- paste0('Transformed_observed_Trait',i)
}

print(head(dat))

write.table(dat,'simulated_prs_transPhe.txt',
	    row.names=F,sep='\t',quote=F)
