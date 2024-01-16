# generate gnova job file for 43 CAD-correlated traits

library(data.table)
library(stringr)

trait <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/trait_43/trait_43_info.csv') # input trait information
path <- fread('/ysm-gpfs/pi/zhao-data/jh2875/corr_CAD/information_update.csv') # input file containing paths to curated file
path43 <- path[path$Code %in% trait$Code,c(4,6)] # subtract path to curated file for 43 traits

# organize path files
# make path43 and trait the same order
trait <- as.data.frame(trait)
path43$Code <- factor(path43$Code,level=trait$Code)
path43 <- path43[order(path43$Code),]

# function to get GNOVA jobs for one trait without overlapping
getJob <- function(num){
  res <- c()
  for(i in 1:num){
    res2 <- paste0("~/anaconda2/bin/python /gpfs/gibbs/pi/zhao/jh2875/GNOVA/gnova.py ", 
		path43$Curated_file[i], " ", 
		path43$Curated_file[num], 
		" --bfile /ysm-gpfs/pi/zhao/yy496/PRS/GWAS_1000G_mac5eur_mapping/1000G_mac5eur",
		" --out  /ysm-gpfs/pi/zhao-data/jh2875/clean/data/gnova_43/corr/",
		trait$Group[i],'_',trait$Code[i],"_",path43$Code[num],'.txt')
    res <- c(res,res2)
  }
  return(res)
}

job <- c()
for(i in 1:43){
  add <- getJob(i)
  job <- c(job,add)
}

# output job file for gnova
write.table(job,
	   '/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/code/gnova_43/gnova_43.job',
	   quote=F,row.names=F,col.names=F)
