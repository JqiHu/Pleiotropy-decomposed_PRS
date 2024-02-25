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
### binary traits
phe_bi <- c('sex','smoking_current','smoking_ever',
	    'diet_status','hypertation','father_heart_disease','mother_heart_disease','sibling_heart_disease',
	    'family_disease_history','DIA','cholesterol_lowering_medication')
phe_cont <- setdiff(phe_name,phe_bi)

# merge prs and phenotype
data <- merge(pheno,ps_prs,
	by.x='eid',by.y='IID',
	all.y=T)


# calculate correlations for certain PD-PRS with phenotypes
getcor_prs <- function(path,mark){ # mark = estimate for coefficients or p.value for p values 
  res <- c()
  ### for continuous traits
  for (i in 1:length(phe_cont)){
    temp <- subset(data,select=c(path,phe_cont[i]))
    temp <- as.matrix(temp)
    add <-  as.numeric(cor.test(temp[,1],temp[,2])[mark]) # calculate correlations between path and phe_name[i]
    res <- c(res, add)
  }
  ### for binary traits
  for (i in 1:length(phe_bi)){
    temp <- subset(data,select=c(path,phe_bi[i]))
    colnames(temp) <- c('prs','phe')
    temp$phe <- factor(temp$phe,levels=unique(sort(temp$phe)),labels=0:1)
    print(table(temp$phe))
    fit <- glm(phe ~ prs,data=temp,family='binomial')
    coef <- summary(fit)$coefficients
    add <- ifelse(mark=='estimate',coef[2,1],coef[2,4])
    res <- c(res,add)
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
rownames(corr_res) <- c(phe_cont,phe_bi) # each row is one phenotype
colnames(corr_res) <- c(ps_names) # each column is one psPRS
rownames(corr_p) <- c(phe_cont,phe_bi)
colnames(corr_p) <- c(ps_names)

# change sequence of row and column names
seq <- match(phe_name,c(phe_cont,phe_bi))
corr_res2 <- corr_res[seq,]
corr_p2 <- corr_p[seq,]

write.table(corr_res2,'../prs_analysis_V2_absolute/corr_prs_phe_coef.tsv',
                row.names=T,col.names=T,quote=F,sep='\t')
write.table(corr_p2,'../prs_analysis_V2_absolute/corr_prs_phe_p.tsv',
                row.names=T,col.names=T,quote=F,sep='\t')
