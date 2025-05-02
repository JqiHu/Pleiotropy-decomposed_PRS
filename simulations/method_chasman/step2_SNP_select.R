# Select independent SNPs in causal regions for trait 1
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/sumstats')
library(data.table)

# Read simulated GWAS
simu <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_beta/simulated_beta_by_region.txt')

# Read pruned SNPs
snp <- fread('../ld_prune/1000G_mac5eur_hapmap3_chr22.prune.in',data.table=F,header=F)

# extract independent SNPs
simu.ind <- simu[simu$marker.ID %in% snp$V1,]

##### Keep SNPs in causal regions for trait1!!!!!
res <- simu.ind[simu.ind$Causal_Trait1==1,]
print(dim(res))

write.table(res,'simu_causalfortrait1.txt',
	    row.names=F,sep='\t',quote=F)
