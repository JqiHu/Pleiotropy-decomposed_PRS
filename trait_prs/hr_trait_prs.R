# generate HR table for trait PRSs on CAD

library('survival')
library('data.table')
library('stringr')

path <- '/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/trait_prs/score'
prs_name <- list.files(path,'profile',full.names=F) # list prs files 
prs <- paste0(path,'/',prs_name) # path to read prs files
prs_name <- str_replace_all(prs_name,'.profile','')

# CAD phenotypes
cad <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/pheno_cad.tsv')
time_cad <- fread('./phenotype/ukbb_CAD.tsv',select=c('eid','age_end')) # update follow-up time
data <- merge(cad,time_cad,by='eid')

## covariates
cov <- fread('./phenotype/pheno_all.tsv',select=c('eid','sex'))
pc <- fread('/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_pheno/MR_time/adjust/adjust.csv') # principle components
pc <- subset(pc,select=c('eid',
                  paste0('PC',1:4)))
cov <- merge(cov,pc,by='eid')

data <- merge(data,cov,by='eid',all.x=T)

# separate models to predict risk for CAD
getHR <- function(prs_path){ # select one prs one time
  res <- subset(fread(prs_path),
		select=c('IID','SCORE')) # keep iid and score
  colnames(res) <- c('eid','score')
  temp <- merge(data,res,by='eid')
  temp$score <- scale(temp$score,
		  center = TRUE, scale = TRUE)
  attach(temp)

  # cox ph model
  sur <- Surv(temp$`age_end`,temp$CAD)
  res_cox <- summary(coxph(sur ~ score + age_recruit + sex + 
			   PC1 + PC2 + PC3 + PC4),
		data=temp)

  # output HR, 95% CI and p value
  res <- cbind(t(res_cox$conf.int[1,]),res_cox$coefficients[1,5])
  return(res)
}

hr <- data.frame()
for(i in 1:length(prs)){
  add <- cbind(prs_name[i],getHR(prs[i]))
  hr <- rbind.data.frame(hr,add)
}
hr <- hr[,-3] # remove -exp(beta)
colnames(hr) <- c('pathway','HR','lower.95','upper.95','p_val')
hr <- hr[order(hr$`HR`,decreasing=T),] # rank by HR
rownames(hr) <- NULL

setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/res/table')
write.csv(hr,
	  'hr_cad_traitprs.csv',
          row.names=F)
