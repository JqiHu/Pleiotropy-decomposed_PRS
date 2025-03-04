# Test for interactions 
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data_interaction/data')
library(data.table)
library(lmtest)

# read into data
dat <- fread('simulated_prs_transPhe.txt')
## INT-transformed traits
trait.exp <- paste0('Transformed_observed_Trait',2:4)
## PD-PRSs
prs <- paste0('simu_PRS_Observed_Beta_Trait1_',c(1:4,'all'))


# Function to get interaction
getInterP <- function(t,p){
  # model 1 with interaction
  formu <- formula(paste0('Observed_Trait1~',t,'*scale(',p,')'))
  fit <- lm(formu,data=dat)
  coef <- summary(fit)$coefficients
  
  res <- coef[4,]
  res <- c(t,p,res)
  return(res)
}

# Run for each trait and prs
out <- c()
for(t in trait.exp){
  for(p in prs){
    add <- getInterP(t,p)
    out <- rbind(out,add)
  }
}
colnames(out)[1:2] <- c('exposure','prs')

### FDR control
out <- as.data.frame(out)
out$p.bh <- p.adjust(out$`Pr(>|t|)`,method='BH')

print(out)

write.csv(out,'inter_test.csv',row.names=F)
