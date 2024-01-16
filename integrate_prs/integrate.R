# combine 9 PD-PRS together
# and compare it to overall CAD PRS
library('data.table')
library('stringr')
library('survival')
library('pROC')

# CAD phenotypes
cad <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/pheno_cad.tsv')
# PD-PRS files
prs <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS_V2_absolute/ps_prs.tsv')
# PRS for CAD
cad_prs <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS/cad_prs.tsv',
		 select=c('IID','SCORESUM'))

# Merge prs and CAD
prs <- merge(cad_prs,prs,by='IID')
colnames(prs)[2] <- 'cad_prs'
#### Standardize CAD PRS
prs$cad_prs <- scale(prs$cad_prs,
	scale=T,center=T)
ps_names <- c('Basic','BP','CVD','Immune','Lipids','Others',
	'Obesity','Respir','T2D')
colnames(prs)[3:11] <- ps_names

data <- merge(cad,prs,by.x='eid',by.y='IID',
	all.y=T)

# Follow-up time for CAD
time_cad <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/ukbb_CAD.tsv',select=c('eid','age_end')) # update follow-up time
data <- merge(data,time_cad,
              by='eid',all.x=T)

#integrate all PD-PRSs into overall PRS for CAD
prs_cox <- as.formula(paste0('Surv(age_end,CAD) ~ ',
	     paste(ps_names,collapse='+')))
beta <- coef(coxph(prs_cox,data=data,id=eid))

## weighted sum of PD-PRSs
data$allprs <- 0
for(i in 1:length(ps_names)){ 
  data$allprs <- data$allprs+
	beta[i]*unlist(subset(data,select=ps_names[i]))
}
data$allprs <- scale(data$allprs,
	scale=T,center=T)

#compare performance of prediction
getHR <- function(name){ # name of column in data
  #formu <- as.formula(paste0('Surv(age_recruit,CAD) ~ ',name))
  formu <- as.formula(paste0('Surv(age_end,CAD) ~ ',name))
  model <- summary(coxph(formu,data=data,id=eid))
  res <- c(name,model$conf.int,
	   model$coefficients[,5])
  return(res)
}

hr <- rbind(getHR('allprs'),getHR('cad_prs'))
hr <- hr[,-3] # remove exp(-coef)

###top 10% vs. bottom 90%
getHR <- function(name){ # name of column in data
  prs <- unlist(subset(data,select=name))
  quant <- quantile(prs,0.9) # top 10%
  data$label <- ifelse(prs>quant,1,0)
  formu <- as.formula('Surv(age_end,CAD) ~ label')
  sum <- summary(coxph(formu,data=data))
  res <- c(name,
	sum$conf.int,sum$coefficients[,5])
  return(res)
}
hr2 <- rbind(getHR('allprs'),getHR('cad_prs'))
hr2 <- hr2[,-3] # remove exp(-coef)

# combine two results
hr_out <- cbind(hr,hr2)
colnames(hr_out) <- rep(c('prs','HR','lower','upper','pval'),2)

write.csv(hr_out,
	  '/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/integrate_prs_overlap/hr_integrateprs_updatefollowup.csv',
	  row.names=F)

## compare prediction performance
getAUC <- function(name){ # name of column in data
  set.seed(123)
  sample <- sample(c(TRUE, FALSE), nrow(data), replace=TRUE, prob=c(0.9,0.1))
  train <- data[sample,]
  test <- data[!sample,]
  formu <- as.formula(paste0('Surv(age_recruit,CAD) ~ ',name))
  fit <- coxph(formu,data=train)
  prediction <- predict(fit,test)
  auc_ci <- ci.auc(test$CAD,prediction)  
  res <- c(auc_ci)
  return(res)
}

auc_out <- rbind(getAUC('allprs'),getAUC('cad_prs'))
auc_out <- cbind(c('iPRS','CAD_PRS'),auc_out)
colnames(auc_out) <- c('PRS','lower','AUC','upper')

write.csv(auc_out,
	  '/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/integrate_prs_overlap/auc_integrateprs_updatefollowup.csv',
	  row.names=F)






