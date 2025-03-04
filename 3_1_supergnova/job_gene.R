# Rscript to generate SUPERGNOVA job

library(data.table)
library(stringr)

trait <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/trait_43/trait_43_info.csv') # input trait information
path <- fread('/ysm-gpfs/pi/zhao-data/jh2875/corr_CAD/information_update.csv') # input file containing paths to curated file
path43 <- unlist(path[path$Code %in% trait$Code,4]) # subtract path to curated file for 43 traits
pathcad <- "/gpfs/ysm/pi/zhao-data/jh2875/munge_gwas/CARDIoGRAMplusC4D_1000G_CAD_2015.sumstats.gz" # path to CAD sum stats
pathsu="/gpfs/gibbs/pi/zhao/jh2875/SUPERGNOVA/" # path to SUPERGNOVA
trait <- trait[order(trait$Code),]
path <- path[order(path$Code),]
name <- paste0(trait$Group,'_',
		trait$Code,'_CAD.txt') # output file name

# generate job for each sum stats
job <- c()
for(i in 1:length(path43)){
  add <- paste0('cd ',pathsu,';',
		'~/anaconda3/bin/python3; ',
		'python supergnova.py ',
		pathcad,' ',path43[i],' ',
		'--bfile /ysm-gpfs/pi/zhao/yy496/PRS/GWAS_1000G_mac5eur_mapping/1000G_mac5eur',
		' --partition data/partition/eur_chr@.bed ',
		'--out /gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/data/supergnova/res/',name[i])
  job <- c(job,add)
}

write.table(job,
	'/gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/code/supergnova/su43.job',
	row.names=F,col.names=F,quote=F)
