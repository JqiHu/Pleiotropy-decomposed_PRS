# R script to generate job file for GNOVA
# calculate correlations between 2,128 traits with CAD

library(data.table)
### 2,128 traits saved in two different folders
path1 <- '/gpfs/ysm/pi/zhao-data/ag2589/VD/More_trait/Munged/'
disease1 <- list.files(path1,"*gz") # name of munged file
disease1_name <- sub('(.gz)$', '', disease) # name of munged traits

path2 <- '/ysm-gpfs/pi/zhao-data/jh2875/munge_gwas/'
disease2 <- list.files(path2,'*.sumstats.gz') # name of munged file
disease2_name <- c(disease2,sub('(.gz)$', '',a)) # name of munged traits

# generate job file for GNOVA
my <- c()
for (i in 1:length(disease1)){
  add <- paste0("~/anaconda2/bin/python /gpfs/gibbs/pi/zhao/jh2875/GNOVA/gnova.py  ",
		path1, disease1[i],
		" /gpfs/ysm/pi/zhao-data/jh2875/munge_gwas/CARDIoGRAMplusC4D_1000G_CAD_2015.sumstats.gz ",
		"--bfile /ysm-gpfs/pi/zhao/yy496/PRS/GWAS_1000G_mac5eur_mapping/1000G_mac5eur ",
		"--out /ysm-gpfs/pi/zhao-data/jh2875/ps_PRS/clean/data/gnova_2128/corr_cad",disease1_name[i],"_CAD.txt")
  my <- c(my,add)
}
for (i in 1:length(disease2)){
  add <- paste0("~/anaconda2/bin/python /gpfs/gibbs/pi/zhao/jh2875/GNOVA/gnova.py  ",
                path2, disease2[i],
                " /gpfs/ysm/pi/zhao-data/jh2875/munge_gwas/CARDIoGRAMplusC4D_1000G_CAD_2015.sumstats.gz ",
                "--bfile /ysm-gpfs/pi/zhao/yy496/PRS/GWAS_1000G_mac5eur_mapping/1000G_mac5eur ",
                "--out /ysm-gpfs/pi/zhao-data/jh2875/ps_PRS/clean/data/gnova_2128/corr_cad",disease2_name[i],"_CAD.txt")
  my <- c(my,add)
}

# output job file
write.table(my,
	   '/gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/code/gnova_2128/all_cad_gnova.job',
	   quote=F,row.names=F,col.names=F)
