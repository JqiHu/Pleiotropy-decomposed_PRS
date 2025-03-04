# Scale beta for twice!!!!!!
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/sumstats')
library(data.table)

# read data
dat0 <- fread('simu_causalfortrait1.txt')
### MAF for dat
maf <- fread('../ld_prune/1000G_mac5eur_hapmap3_chr22_pruned.frq')

dat1 <- merge(dat0,maf,
	      by.x=c('chromosome','marker.ID','allele1','allele2'),
	      by.y=c('CHR','SNP','A1','A2'))# all matched!

# Phenotype data
dat.phe <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_pheno/simulated_UKBphenotype.txt')

# 1. Scaling by phenotype SE
## variances across 4 traits
variance <- apply(subset(dat.phe,select=c(paste0('Observed_Trait',1:4))),
		  2,var) # all 1 bc of standardization
#### SE for phenotype
N.gwas <- 50000 # simulated sample size for GWAS
dat1$se.phe <- sqrt(1/(N.gwas*2*dat1$MAF*(1-dat1$MAF)))

# 2. Scaling by SNP sd.
dat1$snp.sd <- 2*dat1$MAF*(1-dat1$MAF)

#### Deal with beta across 4 traits
for(i in 1:4){
  # Simulated beta and scaling factors
  tmp <- subset(dat1,select=c(paste0('Observed_Beta_Trait',i),'se.phe','snp.sd'))
  
  # Scale
  tmp.scale <- apply(tmp,1,function(x){
		    x[1]/x[2]*x[3]
		  })
  add <- tmp.scale

  dat1 <- cbind.data.frame(dat1,add)
  colnames(dat1)[ncol(dat1)] <- paste0('Scaled_observed_beta_trait',i)
}

write.table(dat1,'scaled_simu_causalfortrait1.txt',
	    row.names=F,sep='\t',quote=F)
