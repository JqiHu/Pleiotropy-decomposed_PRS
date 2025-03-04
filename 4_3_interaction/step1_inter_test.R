# calculate hr for interactions 
# between psPRSs and transformed selected traits
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/interaction_new_pipeline/data')
library(survival)
library(lmtest)
library(stringr)

# import args
args <- commandArgs(trailingOnly=T)
### transformation
trans <- args[1] # INT or Residual
### covariates
covar.group <- args[2] # simple or all

data <- readRDS(paste0('psPRS_phe',trans,'.RDS'))

## input phenotypes needed to be adjusted
data2 <- readRDS('psPRS_phe.RDS')
data2.sub <- subset(data2,select=c('eid',
				   'bmi','smoking','cholesterol_lowering_medication','antihtn'))

data <- merge(data,data2.sub,by='eid')
data$smoking <- factor(data$smoking)

# Binary variables
binary_list <- c("insomnia_2","insomnia_3",
                 "cur_tobacco_smk_1","cur_tobacco_smk_2",
                 "past_tobacco_smk_1","past_tobacco_smk_2","past_tobacco_smk_3",
                 'smoking_current',"smoking_ever")

# calculate interactions in all subjects
getInter <- function(trait,prs){ # calculate interactions between trait and specific psPRSs
    trait2 <- str_remove_all(trait,'trans_')
    trait2 <- str_remove_all(trait2,paste0('_',prs))

    if(covar.group=='all'){
      covar.base <- c('age_std','age_std^2','sex',paste0('PC',1:4),'age_std*sex','age_std^2*sex')
   
      if(trait2 %in% c('Apolipoprotein_A','Apolipoprotein_B','Cholesterol',
                    'HDL_cholesterol','LDL_direct','Lipoprotein_A','Triglycerides')){
        covar <- c(covar.base,'bmi','smoking','cholesterol_lowering_medication')
      }else if(trait2 %in% c('fid_4079','fid_4080')){
        covar <- c(covar.base,'bmi','smoking','antihtn')
      }else if(trait2 %in% binary_list[-c(1:2)]){
        covar <- c(covar.base,'bmi')
      }else{
        covar <- c(covar.base,'bmi','smoking')
      }
     }else{
      covar <- c('age_recruit','sex',paste0('PC',1:4))  
     }

    formu <- as.formula(paste0('Surv(follow.up,CAD) ~ ',
			       prs,'*',trait,'+',
			       paste(covar,collapse='+')))
    formu2 <- as.formula(paste0('Surv(follow.up,CAD) ~ ',
                               prs,'+',trait,'+',
                               paste(covar,collapse='+')))
    cox <- coxph(formu,data=data)
    cox2 <- coxph(formu2,data=data)
    ## Likelihood ratio test
    test <- lrtest(cox,cox2)
    add <- test$`Pr(>Chisq)`[2]

    ## add direction
    beta <- cox$coefficients[9]
    add <- ifelse(beta<0,-add,add)

    add <- c(prs,trait,add)
    return(add)
}

phenotype <- colnames(data)[grep('trans_',colnames(data))] 
### get the corresponding psPRS
# Extract the last component
last_segment <- sapply(strsplit(phenotype, "_"), function(x) tail(x, 1))

inter_out <- data.frame()
for(i in 1:length(phenotype)){
  inter_out <- rbind(inter_out,
		     getInter(phenotype[i],last_segment[i]))
}
colnames(inter_out) <- c('psPRS','trait','p-value')


# output
write.csv(inter_out,
	  paste0('interaction_',trans,'_adj_all_',covar.group,'.csv'),
  	  row.names=T)

