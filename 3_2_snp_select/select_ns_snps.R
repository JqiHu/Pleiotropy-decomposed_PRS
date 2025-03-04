# extract non-specific SNPs
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_beta_selection_final_V2_absolute/') # path to SNPs that are assigned to one of the pleiotropy-decomposed regions
library('data.table')
library(stringr)

files <- list.files('.')
beta <- fread('/gpfs/gibbs/pi/zhao/zhao-data/yy496/pathway/ps_prs/data/CAD_AnnoPred.txt')

# input specific SNPs and find non-specific ones
snp_ns <- setdiff(beta$sid,fread(files[1])$sid)
for(i in 2:length(files)){
  snp_ns <- setdiff(snp_ns,fread(files[i])$sid)
  print(length(snp_ns))
}
beta_out <- beta[beta$sid %in% snp_ns,] 
beta_out$chrom <- str_replace_all(beta_out$chrom,'chrom_','')

print(nrow(beta_out))
# output 
write.table(beta_out,
	    'ns_snp.txt',
            col.names=T,row.names=F,sep=' ',quote=F)

