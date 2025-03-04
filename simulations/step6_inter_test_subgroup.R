# Test the interaction between specific pairs of interest
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data_interaction/data')
library(data.table)

# read data
dat <- fread('simulated_subgroup.txt',data.table=F)
#### read transformed phenotype values
dat.phe <- fread('simulated_prs_transPhe.txt',select=c('IID',paste0('Transformed_observed_Trait',2:4)))
dat <- merge(dat,dat.phe,by='IID')

# Function to get interactions
getInter <- function(subgroup,trait){
  tmp <- subset(dat,select=c(subgroup,trait,'Observed_Trait1'))
  colnames(tmp)[1:2] <- c('group','trait')
  tmp <- tmp[complete.cases(tmp),]

  # model for interaction
  fit <- lm('Observed_Trait1~group*trait',data=tmp)
  coef <- summary(fit)$coefficients
  add <- coef[4,]

  res <- c(subgroup,trait,add)
  return(res)
}

# Pairs of interest 
inter.res.sig <- data.frame(subgroup=c('subgroup1','subgroup2','subgroup3'),
			    trait=c('Transformed_observed_Trait2','Transformed_observed_Trait3','Transformed_observed_Trait4'))
out <- c()
for(i in 1:nrow(inter.res.sig)){
  add <- getInter(inter.res.sig$subgroup[i],
		  inter.res.sig$trait[i])
  out <- rbind(out,add)
}
print(out)

write.csv(out,'subgroup_inter_test_or.csv',row.names=F)
