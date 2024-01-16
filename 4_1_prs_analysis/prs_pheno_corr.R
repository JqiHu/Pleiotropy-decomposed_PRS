# evaluate correlations between PD-PRSs and phenotypes

library('data.table')
library('stringr')

setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS_V2_absolute/') # path to PD-PRSs

ps_prs <- fread('ps_prs.tsv')
ps_names <- c('Basic Condition','BP','CVD','Immune System',
              'Lipids','Others','Obesity','Respiratory system','T2D') # name for pathways
colnames(ps_prs)[-1] <- ps_names 

# phenotype file
pheno <- fread('../phenotype/pheno_all.tsv')
phe_name <- colnames(pheno)[-1]

# Add fat percentages
update <- readRDS('../phenotype/ukbb_all_modifiable.RDS') # original
update <- as.data.frame(update)
update <- update[,-which(duplicated(colnames(update))==T)]
update <- merge(pheno[,1:2],update,by='eid')
### Use the mean value of left and right arm/leg
update$arm_fat_percentage <- (update$`23111`+update$`23115`)/2
update$leg_fat_percentage <- (update$`23119`+update$`23123`)/2
update_add <- subset(update,select=c('eid','23099','23127','arm_fat_percentage','leg_fat_percentage'))
colnames(update_add)[2:3] <- c('body_fat_percentage','trunk_fat_percentage')

pheno <- merge(pheno,update_add,by='eid',all.x=T)
phe_name <- colnames(pheno)[-1]

# merge prs and phenotype
data <- merge(pheno,ps_prs,
	by.x='eid',by.y='IID',
	all.y=T)

# calculate correlations for certain PD-PRS with phenotypes
getcor_prs <- function(path,mark){ # mark = estimate for coefficients or p.value for p values 
  res <- c()
  for (i in 1:(ncol(pheno)-1)){
    temp <- subset(data,select=c(path,phe_name[i]))
    temp <- as.matrix(temp)
    add <-  as.numeric(cor.test(temp[,1],temp[,2])[mark]) # calculate correlations between path and phe_name[i]
    res <- c(res, add)
  }
  return(res)
}


corr_res <- getcor_prs(ps_names[1],'estimate') # coefficient matrix
corr_p <- getcor_prs(ps_names[1],'p.value') # p value matrix
for(i in 2:length(ps_names)){
  corr_res <- cbind.data.frame(corr_res,
		getcor_prs(ps_names[i],'estimate'))
  corr_p <- cbind.data.frame(corr_p,
              getcor_prs(ps_names[i],'p.value'))
}
rownames(corr_res) <- phe_name # each row is one phenotype
colnames(corr_res) <- ps_names # each column is one psPRS
rownames(corr_p) <- phe_name
colnames(corr_p) <- ps_names

write.table(corr_res,'../prs_analysis_V2_absolute/corr_prs_phe_coef.tsv',
                row.names=T,col.names=T,quote=F,sep='\t')
write.table(corr_p,'../prs_analysis_V2_absolute/corr_prs_phe_p.tsv',
                row.names=T,col.names=T,quote=F,sep='\t')
