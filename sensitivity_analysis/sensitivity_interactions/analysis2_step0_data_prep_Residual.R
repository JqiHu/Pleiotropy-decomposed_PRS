# calculate hr for interactions 
# between psPRSs and transformed selected traits
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/interaction_new_pipeline/data')
library(survival)
library(lmtest)

data <- readRDS('psPRS_phe.RDS')
#### Confounding controlled
#### And Rank-based inverse normal transformed 
##### Function to do this for specific continuous traits
##### and psPRS pair
getINT_cont <- function(trait,psprs){
  temp <- data[is.na(data[trait])==F,]
  temp <- temp[is.na(temp$smoking)==F & is.na(temp$cholesterol_lowering_medication)==F &
	       is.na(temp$bmi)==F,]
  if(trait %in% c('Apolipoprotein_A','Apolipoprotein_B','Cholesterol',
		  'HDL_cholesterol','LDL_direct','Lipoprotein_A','Triglycerides')){
    formu <- formula(paste0(trait,'~ age_std+age_std^2+sex+PC1+PC2+PC3+PC4+bmi+
			  smoking+age_std*sex+age_std^2*sex+cholesterol_lowering_medication+',psprs))
  }else if(trait %in% c('fid_4079','fid_4080')) {
    formu <- formula(paste0(trait,'~ age_std+age_std^2+sex+PC1+PC2+PC3+PC4+bmi+
                          smoking+age_std*sex+age_std^2*sex+antihtn+',psprs))
  }else {
    formu <- formula(paste0(trait,'~ age_std+age_std^2+sex+PC1+PC2+PC3+PC4+bmi+
                          smoking+age_std*sex+age_std^2*sex+',psprs))
  }
  fit <- lm(formu,data=temp)
  residuals <- fit$residuals
  out <- cbind(temp$eid,residuals)
  colnames(out) <- c('eid',paste0('trans_',trait,'_',psprs))

  print(paste0('trans_',trait,'_',psprs))
  return(out)
}

# All continuous traits
lipids_list <- c( 'Apolipoprotein_A','Apolipoprotein_B','Cholesterol',
                 'HDL_cholesterol','LDL_direct','Lipoprotein_A','Triglycerides',
                 'fid_1160','fid_1269','fid_1279','alcohol_intake',
                 'fid_3062','fid_3063','fid_3064','fid_4079','fid_4080',
                 'fid_23099','leg_fat_percentage','arm_fat_percentage','fid_23127',
                 'fid_189',"qualification",'Vitamin_D')
# 9 psPRSs + overall PRS
psprs <- c('Basic','BP','CVD','Immune','Lipids','Others',
	   'Obesity','Respir','T2D','overall')

for(i in 1:length(lipids_list)){
  for(j in 1:length(psprs)){
    add <- getINT_cont(lipids_list[i],psprs[j])
    data <- merge(data,add,by='eid',all.x=T)
  }
}

##### Function to do this for specific binary traits
binary_list <- c("insomnia_2","insomnia_3",
		 "cur_tobacco_smk_1","cur_tobacco_smk_2",
		 "past_tobacco_smk_1","past_tobacco_smk_2","past_tobacco_smk_3",
		 'smoking_current',"smoking_ever")
getINT_bi <- function(trait,psprs){
  temp <- data[is.na(data[trait])==F,]
  temp <- temp[is.na(temp$smoking)==F & is.na(temp$cholesterol_lowering_medication)==F &
               is.na(temp$bmi)==F,]
  if(trait %in% c("insomnia_2","insomnia_3")){ 
    formu <- formula(paste0(trait,'~ age_std+age_std^2+sex+PC1+PC2+PC3+PC4+bmi+
                          smoking+age_std*sex+age_std^2*sex+',psprs))
  }else{
    formu <- formula(paste0(trait,'~ age_std+age_std^2+sex+PC1+PC2+PC3+PC4+bmi+
                          age_std*sex+age_std^2*sex+',psprs))
  }
  fit <- glm(formu,data=temp,family='binomial')
  residuals <- residuals(fit, type = "working")

  out <- cbind(temp$eid,residuals)
  colnames(out) <- c('eid',paste0('trans_',trait,'_',psprs))

  print(paste0('trans_',trait,'_',psprs))
  return(out)
}
for(i in 1:length(binary_list)){
  for(j in 1:length(psprs)){
    add <- getINT_bi(binary_list[i],
		     psprs[j])
    data <- merge(data,add,by='eid',all.x=T)
  }
}
data2 <- data[,-c(4:34,51:60)]

saveRDS(data2,'psPRS_pheResidual.RDS')

