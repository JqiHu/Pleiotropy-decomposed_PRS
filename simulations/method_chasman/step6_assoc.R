# Associations between PRSs and phenotypes
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/analysis')
library(data.table)

# read data
#dat <- fread('simulated_phe_decomposed_prs.txt',data.table=F)
dat <- fread('simulated_phe_inter_decomposed_prs.txt',data.table=F)
var.exp <- colnames(dat)[grep('PRS',colnames(dat))]
var.outcome <- paste0('Observed_Trait',1:4)

out <- c()
for(e in var.exp){
  for(o in var.outcome){
    formu <- formula(paste0(o,'~',e))
    fit <- lm(formu,data=dat)
    coef <- summary(fit)$coefficients
    add <- coef[2,]
    add <- c(e,o,add)
    out <- rbind(out,add)
  }
}
colnames(out)[1:2] <- c('PRS','Outcome')
print(out)

write.csv(out,'inter_assoc_prs.csv',row.names=F)
