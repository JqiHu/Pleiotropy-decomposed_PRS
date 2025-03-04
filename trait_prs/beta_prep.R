# temparory file to format gwas sum stats
setwd('/ysm-gpfs/pi/zhao-data/jh2875/ps_PRS/clean/data/trait_prs/gwas/')

file <- list.files('.')
for(gwas in file){
#data <- fread(paste0('~/scratch60/trait_PRS/data/',gwas),header=T)
data <- fread(gwas,header=T)
#colnames(data) <- str_replace(colnames(data),'Z','BETA') 
data$P <- 2*(pnorm(abs(data$BETA),lower.tail=F))
#data <- subset(data,select=c('SNP','A1','A2','BETA','P'))
write.table(data,paste0('/ysm-gpfs/pi/zhao-data/jh2875/ps_PRS/clean/data/trait_prs/gwas/',gwas),row.names=F,col.names=T,quote=F)
}
