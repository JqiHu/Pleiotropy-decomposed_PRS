# create a new categorical variable
# stratified by PRS and subgroups
library(stringr)
library(data.table)
library(survival)

# Import threshold for subgroup
args <- commandArgs(trailingOnly=T)
thresh <- as.character(args[1])
if(thresh=='5%'){
  ### 5%
  setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_analysis_V2_absolute/')
}else if(thresh=='10%'){
  ### 10%
  setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/subgroup_sensitivity/data/0.1')
}else{
  ### 1%
  setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/subgroup_sensitivity/data/0.01')
}


# input data
data <- readRDS('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/interaction_new_pipeline/data/psPRS_phe.RDS')

##### -------------------------------- #####
### Dichotomize specific phenotypes
#### sleep duration
data$sleep_bi <- 0
data$sleep_bi[data$`fid_1160`<6] <- 1
#### FVC, FEV1, PEF
data$FVC_bi[is.na(data$`fid_3062`)==F] <- 0
data$FVC_bi[is.na(data$`fid_3062`)==F &
            data$`fid_3062`>median(data$`fid_3062`,na.rm=T)] <- 1
data$FEV1_bi[is.na(data$`fid_3063`)==F] <- 0
data$FEV1_bi[is.na(data$`fid_3063`)==F &
            data$`fid_3063`>median(data$`fid_3063`,na.rm=T)] <- 1
data$PEF_bi[is.na(data$`fid_3064`)==F] <- 0
data$PEF_bi[is.na(data$`fid_3064`)==F &
            data$`fid_3064`>median(data$`fid_3064`,na.rm=T)] <- 1
#### Triglyceride
data$tg_bi[is.na(data$Triglycerides)==F] <- 0
data$tg_bi[is.na(data$Triglycerides)==F &
           data$Triglycerides>2.3] <- 1
data$tg_bi[data$cholesterol_lowering_medication==1] <- 1
#### TC 
data$tc_bi[is.na(data$Cholesterol)==F] <- 0
data$tc_bi[is.na(data$Cholesterol)==F &
           data$Cholesterol>5.17] <- 1
data$tc_bi[data$cholesterol_lowering_medication==1] <- 1
#### LDL
data$ldl_bi[is.na(data$LDL_direct)==F] <- 0
data$ldl_bi[is.na(data$LDL_direct)==F &
              data$LDL_direct>5.6] <- 1
data$ldl_bi[data$cholesterol_lowering_medication==1] <- 1
#### vitamin D
data$VD_bi[is.na(data$Vitamin_D)==F] <- 0
data$VD_bi[is.na(data$Vitamin_D)==F &
           data$Vitamin_D<30] <- 1
#### BP
data$dbp_bi[is.na(data$fid_4079)==F] <- 0
data$dbp_bi[is.na(data$fid_4079)==F &
            data$fid_4079>=90] <- 1
data$dbp_bi[data$antihtn==1] <- 1
####### SBP
data$sbp_bi[is.na(data$fid_4080)==F] <- 0
data$sbp_bi[is.na(data$fid_4080)==F &
            data$fid_4080>=140] <- 1
data$sbp_bi[data$antihtn==1] <- 1

# input high-risk individuals
highrisk <- readRDS('high_risk_ind.RDS')

# name of pathways
files <- list.files('.','subgroup')
files <- files[grep('RDS',files)]
psprs <- str_replace_all(files,'subgroup_','')
psprs <- str_replace_all(psprs,'.RDS','')

# Binary variables
binary_list <- c("insomnia_2","insomnia_3",
                 "cur_tobacco_smk_1","cur_tobacco_smk_2",
                 "past_tobacco_smk_1","past_tobacco_smk_2","past_tobacco_smk_3",
                 'smoking_current',"smoking_ever")

# calculate interactions for specific traits and specific psPRS
getInter <- function(trait,prs_name,covar.index){ 
    if(covar.index=='all'){
      covar.base <- c('age_std','age_std^2','sex',paste0('PC',1:4),'age_std*sex','age_std^2*sex')
      
      if(trait %in% c('Apolipoprotein_A','Apolipoprotein_B','Cholesterol',
                    'HDL_cholesterol','LDL_direct','Lipoprotein_A','Triglycerides','ldl_bi','tc_bi','tg_bi')){
        covar <- c(covar.base,'bmi','smoking','cholesterol_lowering_medication')
      }else if(trait %in% c('fid_4079','fid_4080','dbp_bi','sbp_bi')){
        covar <- c(covar.base,'bmi','smoking','antihtn')
      }else if(trait %in% binary_list[-c(1:2)]){
        covar <- c(covar.base,'bmi')
      }else{
        covar <- c(covar.base,'bmi','smoking')
      }
     }else{
      covar <- c('age_recruit','sex',paste0('PC',1:4))
     }

    temp <- subset(data,select=c('eid','follow.up','CAD',trait,
                                 setdiff(colnames(data),c('eid','follow.up','CAD',trait))))
    colnames(temp)[4] <- 'Phenotype'
    temp <- temp[is.na(temp$Phenotype)==F,]

    # high-risk subjects
    temp <- temp[temp$eid %in% highrisk$IID,]

    # read subgroup data
    i=grep(prs_name,psprs)
    subgroup <- readRDS(files[i])
    temp$index <- 0
    temp$index[temp$eid %in% subgroup$IID] <- 1
    print(table(temp$index))
    
    # create a new variable
    temp$index[temp$Phenotype==1] <- temp$index[temp$Phenotype==1]+2 
    temp$index <- factor(temp$index,levels=0:3) #0: remains and phe=0[REF]; 1: subgroup and phe=0; 2: remains and phe=1; 3: subgroup and phe=1

    # Count of CAD in each group
    table <- aggregate(temp$CAD,list(temp$index),table)

    # Fit Cox
    formu <- formula(paste0('Surv(follow.up,CAD) ~ index+',paste(covar,collapse='+')))
    cox <- summary(coxph(formu,data=temp))
    out <- cox$coefficients[1:3,c(1:3,5)]
    out <- rbind(c(0,1,0,0),out)
    add <- cbind(table,out)

    return(add)
}


### Manually select traits of intersct (significant ones in step1!!!!!!)
phenotype <- c('FVC_bi','FEV1_bi','sleep_bi',
	       'smoking_current','cur_tobacco_smk_1','cur_tobacco_smk_2',
	       'tg_bi','ldl_bi','smoking_current',
	       'FVC_bi','FEV1_bi',
	       'cur_tobacco_smk_1',
               'cur_tobacco_smk_1','smoking_current')
prs_names <- c(rep('BP',3),
               rep('Immune',3),
               rep('Lipids',3),
               rep('Obesity',2),
               'Others',
               rep('Respiratory_system',2))

res <- list()
for(i in 1:length(phenotype)){
  res[[i]] <- getInter(phenotype[i],prs_names[i],'all')
}
names(res) <- paste0(phenotype,'*',prs_names)

saveRDS(res,paste0('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/interaction_new_pipeline/data/inter_select_HR_vis_',thresh,'_covarAll.RDS'))

res <- list()
for(i in 1:length(phenotype)){
  res[[i]] <- getInter(phenotype[i],prs_names[i],'simple')
}
names(res) <- paste0(phenotype,'*',prs_names)

setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/interaction_new_pipeline/data')
saveRDS(res,paste0('inter_select_HR_vis_',thresh,'_covarSimple.RDS'))

