# aggregate correlation results

library('data.table')
library('stringr')

setwd('/gpfs/ysm/pi/zhao-data/jh2875/ps_PRS/clean/data/gnova_2128/corr_cad') # path to gnova output
a <- list.files('./results','_CAD.txt')
res <- data.frame()
for(i in 1:length(a)){
  # input gnova output and combine together
  res <- rbind.data.frame(res,fread(paste0('./results/',a[i]),header=T))
}

# remove suffix for trait names
a <- str_replace(a,'_CAD.txt','')
a <- str_replace(a,'.txt','')
a <- str_replace(a,'.tbl','')
a <- str_replace(a,'.csv','')

res$trait <- a

# FDR control
res <- res[order(res$pvalue,decreasing=F),]
threshold <- (1:nrow(res))*0.05/nrow(res)
res <- cbind(res,threshold)
fdr_reject <- ifelse(res$pvalue < res$threshold,1,0)

# output significant results
res_fdr <- res[which(fdr_reject==1),]
write.table(res_fdr,
	'../gnova_res_fdr.txt',
	row.names=F,col.names=T,sep=' ',quote=F)

