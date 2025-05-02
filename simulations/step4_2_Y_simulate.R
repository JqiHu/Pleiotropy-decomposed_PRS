# Simulate phenotypes
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data_interaction/data/')
library(data.table)
library(stringr)

# Import phe1
dat.phe1 <- c()
num_trait=4
for(i in 1:num_trait){
  tmp2 <- fread(paste0('../../simu_data/data/simu_pheno/phe1/simu_UKB_sum_PRS_Observed_Beta_Trait',i,'.profile'),
                select=c('FID','IID','SCORESUM'))
  colnames(tmp2)[3] <- paste0('Observed_Trait',i)

  add <- tmp2
  if(i==1){
    dat.phe1 <- add
  }else{
    dat.phe1 <- merge(dat.phe1,add,by=c('FID','IID'))
  }
}
##### standardize phe1
dat.phe1 <- as.data.frame(dat.phe1)
dat.phe1[,-c(1:2)] <- apply(dat.phe1[,-c(1:2)],2,scale)

###### First, simulate Trait 2 to Trait 4
#### Official simulation across 8 phenotypes
set.seed(123)  # For reproducibility
dat.final <- dat.phe1[,1:2]
for(i in 2:4){
  y_genetic <- unlist(subset(dat.phe1,select=paste0('Observed_Trait',i)))
  # Step 1: Simulate Y2 as random noise
  n_individuals <- nrow(dat.phe1)  # Number of individuals
  y_noise <- rnorm(n = n_individuals, mean = 0, sd = 1)  # Random noise Y2

  # Step 2: Standardize genetic and noise components
  y_genetic <- scale(y_genetic)
  y_noise <- scale(y_noise)

  # Step 3: Define fixed beta coefficients
  beta1 <- 2  # Genetic effect
  beta2 <- 1  # Noise effect
  beta3 <- 0.3  # Interaction effect

  # Step 4: Compute interaction term
  interaction_term <- y_genetic * y_noise

  # Step 5: Compute final phenotype
  y_final <- beta1 * y_genetic + beta2 * y_noise + beta3 * interaction_term
  y_final <- scale(y_final)

  # Combine into a data frame
  dat.final <- cbind.data.frame(dat.final,y_noise,y_final)
  colnames(dat.final)[c(ncol(dat.final)-1,ncol(dat.final))] <- c(paste0('Noise_Trait',i),paste0('Observed_Trait',i))
}

##### align IIDs
dat.phe1 <- dat.phe1[match(dat.final$IID,dat.phe1$IID),]
#### Simulate Trait 1 by adding trait 2 to trait 4
y_1 <- scale(dat.phe1$Observed_Trait1)
y_noise <- scale(rnorm(n = n_individuals, mean = 0, sd = 1)) # random noise

###### add PD-PRSs
pdprs <- fread('../../pd_prs/data/simulated_PD_PRS_all.txt',data.table=F)
pdprs <- pdprs[match(dat.final$IID,pdprs$eid),]
#### extract pd-prss
y_1.2 <- scale(pdprs$simu_PRS_Observed_Beta_Trait1_1)
y_1.3 <- scale(pdprs$simu_PRS_Observed_Beta_Trait1_2)
y_1.4 <- scale(pdprs$simu_PRS_Observed_Beta_Trait1_3)

beta1 <- 2 # genetic
beta2 <- 1 # Trait 2
beta3 <- 1 # Trait 3
beta4 <- 1 # Trait 4
beta5 <- 1 # noise
beta6 <- 0.3 # interaction with Trait 2 
beta7 <- 0.3 # interaction with Trait 3
beta8 <- 0.3 # interaction with Trait 4
beta9 <- 0.3 # interaction with noise

y_final <- beta1*y_1 + 
	beta2*scale(dat.final$Observed_Trait2) + 
	beta3*scale(dat.final$Observed_Trait3) +
	beta4*scale(dat.final$Observed_Trait4) +
	beta5*y_noise +
	beta6*y_1.2*scale(dat.final$Observed_Trait2) +
	beta7*y_1.3*scale(dat.final$Observed_Trait3) +
	beta8*y_1.4*scale(dat.final$Observed_Trait4) +
	beta9*y_1*y_noise

dat.final <- cbind.data.frame(dat.final,y_final)
colnames(dat.final)[ncol(dat.final)] <- 'Observed_Trait1'

write.table(dat.final,'simu_pheno_withInter.txt',
	    row.names=F,sep='\t',quote=F)

