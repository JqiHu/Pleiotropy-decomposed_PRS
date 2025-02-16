# calculate hr for interactions 
# between psPRSs and transformed selected traits
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/interaction_new_pipeline/data')
library(survival)
library(lmtest)
library(stringr)

# import args
data <- readRDS('psPRS_pheINT_baselineCov.RDS')

# Binary variables
binary_list <- c("insomnia_2","insomnia_3",
                 "cur_tobacco_smk_1","cur_tobacco_smk_2",
                 "past_tobacco_smk_1","past_tobacco_smk_2","past_tobacco_smk_3",
                 'smoking_current',"smoking_ever")

# calculate interactions in all subjects
getInter <- function(trait,prs){ # calculate interactions between trait and specific psPRSs
    trait2 <- str_remove_all(trait,'trans_')
    trait2 <- str_remove_all(trait2,paste0('_',prs))

    covar <- c('age_recruit','sex',paste0('PC',1:4))  

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
write.csv(inter_out,'interaction_INT_adj_baseline.csv',row.names=T)

