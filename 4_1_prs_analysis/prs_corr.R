# evaluate correlations between PD-PRSs

library('data.table')
library('stringr')

setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS_V2_absolute/') # path to PD-PRSs

ps_prs <- fread('ps_prs.tsv')
ps_names <- c('Basic','BP','CVD','Immune',
	      'Lipids','Others','Obesity','Respir','T2D') # name for pathways
colnames(ps_prs)[-1] <- ps_names

# calculate correlations for certain PD-PRS with other PD-PRSs
getcor_prs <- function(path,mark){ # mark = estimate for coefficients or p.value for p values
  res <- c()
  for (i in 1:(ncol(ps_prs)-1)){
    temp <- subset(ps_prs,select=c(path,ps_names[i]))
    temp <- as.matrix(temp)
    add <- as.numeric(cor.test(temp[,1],temp[,2])[mark]) # calculate correlations between path and ps_name[i] 
    res <- c(res, add)
  }
  return(res)
}

corr_res <- getcor_prs(ps_names[1],'estimate') # coefficient matrix
corr_p <- getcor_prs(ps_names[1],'p.value') # p value matrix
for (i in 2:length(ps_names)) {
  corr_res <- cbind.data.frame(corr_res,
	        getcor_prs(ps_names[i],'estimate'))
  corr_p <- cbind.data.frame(corr_p,
	      getcor_prs(ps_names[i],'p.value'))
}
row.names(corr_res) <- ps_names
colnames(corr_res) <- ps_names
row.names(corr_p) <- ps_names
colnames(corr_p) <- ps_names

# output correlations
write.csv(corr_res,'../prs_analysis_V2_absolute/corr_prs_coef.csv')
write.csv(corr_p,'../prs_analysis_V2_absolute/corr_prs_p.csv')
