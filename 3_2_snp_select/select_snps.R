# form SNPs for 43 traits 

args = commandArgs(trailingOnly=TRUE)
corr_S = args[1] # path to correlation file
factor = args[2] # phenotype name

library('data.table')
library('stringr')

# input PRS for CAD
beta <- fread('/gpfs/gibbs/pi/zhao/zhao-data/yy496/pathway/ps_prs/data/CAD_AnnoPred.txt')
beta$chrom <- str_replace(beta$chrom,'chrom_','')
beta$chom <- as.numeric(beta$chrom)
beta$pos <- as.numeric(beta$pos)

# input correlation file
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/supergnova/res') # path to SUPERGNOVA results
corr <- fread(corr_S)
corr <- as.data.frame(corr)
corr$start <- as.numeric(corr$start)
corr$end <- as.numeric(corr$end)
corr$p <- as.numeric(corr$p)
corr <- corr[order(corr$p,decreasing=F),]

# extract SNPs in phenotype-CAD correlated regions
p_corr <- corr[which(corr$p <= 0.05),] # regions where correlation p<=.05
i <- nrow(p_corr) # number of significant regions

res <- c() # combination of row number and region p value
index <- nrow(p_corr)
for (j in 1:index){
  index_pos <- unique(which(beta$pos < p_corr$end[j] # row number for SNPs in this region
                            &beta$pos > p_corr$start[j]
                            &beta$chrom==p_corr$chr[j]))
  index_p <- rep(p_corr$p[j],length(index_pos)) # p values for SNPs in this region
  index_rho <- rep(p_corr$rho[j],length(index_pos))
  res <- rbind(res,cbind(index_pos,index_p,index_rho))
}

res <- res[order(res[,2]),] # rank SNPs from lowest p to highest
res_unique <- res[which(duplicated(res[,1])==F),] # remove duplicated SNPs
num_snp <- nrow(res_unique) # number of SNPs

beta_gather <- beta[res_unique[,1],] # extract betas for these SNPs
# combine beta with p-values and genetic covariance 
beta_gather <- cbind(beta_gather,res_unique[,2],res_unique[,3])

# check uniqueness
print(length(unique(beta_gather$sid)))
print(nrow(beta_gather))

# output file
write.table(beta_gather,
		paste0('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_beta_selection_V2/',
			factor,'_beta_',num_snp,
			'_sigregion_',i,'.txt'),
			quote=F,row.names=F)
