# generate job for PRScs

library(stringr)

path <- '/ysm-gpfs/pi/zhao-data/jh2875/ps_PRS/clean/data/trait_prs/gwas/' # path to formated sum stats for 8 traits
files <- list.files(path,'*sumstats',
		full.names=F)
size <- c(795640,210967,898130,135638,
	106736,130347,132236,74053) # sample size for GWASs
trait <- c('BMI','HDL','T2D',
	'Mothers_age_at_death','Cerebravas','hypertension',
	'PBC','eversmk') # name for trait

# generate job for PRScs
out <- c()
for(j in 1:length(files)){
  for(i in 1:22){
    add <- paste0('python /gpfs/gibbs/pi/zhao/jh2875/PRScs/PRScs.py ',
		'--ref_dir=/ysm-gpfs/pi/zhao-data/yy496/PRScs/eur_1kg/ldblk_1kg_eur ',
		'--bim_prefix=/ysm-gpfs/pi/zhao-data/yy496/ukbb_v3/ukbb3_neale/ukbb3_imp_nealqc/tmp_Ukb_imp_v3_neal ',
		'--chrom=',i,
		' --sst_file=',path,files[j],
		' --n_gwas=',size[j],
		' --out_dir=/ysm-gpfs/pi/zhao-data/jh2875/ps_PRS/clean/data/trait_prs/beta/',trait[j])
    out <- c(out,add)  
  }
}

# output
write.table(out,
	'/gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/code/prs_analysis/trait_prs/beta.job',
	row.names=F,quote=F,col.names=F)
