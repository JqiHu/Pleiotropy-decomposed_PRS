# create a new categorical variable
# stratified by PRS and subgroups
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_analysis_V2_absolute/')
library(stringr)
library(data.table)
library(survival)

# input data
phe_all <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/pheno_all.tsv')
time <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/phe_interaction.tsv',
              selec=c('eid','follow.up'))
phe <- readRDS('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/ukbb_modifiable_final.RDS')
colnames(phe)[2:76] <- paste0('fid_',colnames(phe)[2:76])
cad <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/pheno_cad.tsv') # cad 
pc <- fread('/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_pheno/MR_time/adjust/adjust.csv') # principle components
pc <- subset(pc,select=c('eid',
                  paste0('PC',1:4)))
cad <- cad[,-2] # remove age

phe_sub <- subset(phe_all,select=c('eid','sex','age_recruit','cholesterol_lowering_medication',
                                   'Apolipoprotein_A','Apolipoprotein_B',
                                   'Cholesterol','HDL_cholesterol','LDL_direct',
                                   'Lipoprotein_A','Triglycerides',
                                   'bmi','smoking_current','smoking_ever','smoking',
                                   'Vitamin_D'))
phe_sub2 <- subset(phe,select=c('eid',paste0('fid_',c(1160,1200,1239,1249,1269,1279,1558,
                                              3062,3063,3064,4079,4080,23099,23127,
                                              189,6138)),
                                'leg_fat_percentage','arm_fat_percentage'))

phe_sub$smoking <- factor(phe_sub$smoking,levels=c(0:2))
pheno <- merge(phe_sub,phe_sub2,by='eid')
pheno <- merge(pheno,time,by='eid')
pheno <- merge(pheno,cad,by='eid')
pheno <- merge(pheno,pc,by='eid')

psprs <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS_V2_absolute//ps_prs.tsv')
colnames(psprs)[c(2,5,9)] <- c('Basic','Immune','Respir')
overallprs <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS/cad_prs.tsv')
overallprs <- subset(overallprs,select=c('IID','SCORESUM'))
colnames(overallprs)[2] <- 'overall'
overallprs$overall <- scale(overallprs$overall,
        center=T,scale=T)
psprs <- merge(psprs,overallprs,by='IID')

ps_names <- colnames(psprs)[-1]

data <- merge(pheno,psprs,
        by.x='eid',by.y='IID',
        all.y=T)

#### Anti-hypertension medicine
anti_htn <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/phe_interaction_antihtn.tsv')
data <- merge(data,anti_htn,by='eid',all.x=T)

### Deal with lipids phenotypes
data <- as.data.frame(data)
#### Adjust medication
data$LDL_direct[data$cholesterol_lowering_medication==1 & 
                is.na(data$LDL_direct)==F] <- data$LDL_direct[data$cholesterol_lowering_medication==1 & is.na(data$LDL_direct)==F]/0.7
data$Cholesterol[data$cholesterol_lowering_medication==1 &
                    is.na(data$Cholesterol)==F] <- data$Cholesterol[data$cholesterol_lowering_medication==1 &  is.na(data$Cholesterol)==F]/0.8

### Recode some of the variables
#### Insomnia
data$insomnia_2[data$`fid_1200`==1] <- 0
data$insomnia_2[data$`fid_1200`==2] <- 1
data$insomnia_3[data$`fid_1200`==1] <- 0
data$insomnia_3[data$`fid_1200`==3] <- 1
#### current tobacco smoking
data$cur_tobacco_smk_1[data$`fid_1239`==0] <- 0
data$cur_tobacco_smk_1[data$`fid_1239`==1] <- 1
data$cur_tobacco_smk_2[data$`fid_1239`==0] <- 0
data$cur_tobacco_smk_2[data$`fid_1239`==2] <- 1
#### past tobacco smoking
data$past_tobacco_smk_1[data$`fid_1249`==4] <- 0
data$past_tobacco_smk_1[data$`fid_1249`==1] <- 1
data$past_tobacco_smk_2[data$`fid_1249`==4] <- 0
data$past_tobacco_smk_2[data$`fid_1249`==2] <- 1
data$past_tobacco_smk_3[data$`fid_1249`==4] <- 0
data$past_tobacco_smk_3[data$`fid_1249`==3] <- 1
#### alcohol intake
data$alcohol_intake <- as.numeric(as.character(data$`fid_1558`))
#### qualification/education level
data$qualification <- as.numeric(as.character(data$`fid_6138`))

#### Remove outliers +- 5 sd for continuous variables
lipids_list <- c( 'Apolipoprotein_A','Apolipoprotein_B','Cholesterol',
                 'HDL_cholesterol','LDL_direct','Lipoprotein_A','Triglycerides',
                 'fid_1160','fid_1269','fid_1279','alcohol_intake',
                 'fid_3062','fid_3063','fid_3064','fid_4079','fid_4080',
                 'fid_23099','leg_fat_percentage','arm_fat_percentage','fid_23127',
                 'fid_189',"qualification",'Vitamin_D')
for(i in 1:length(lipids_list)){
  temp <- unlist(data[lipids_list[i]])
  mean_value <- mean(temp,na.rm=T)
  sd_value <- sd(temp,na.rm=T)
  out_row <- which(temp>(mean_value+5*sd_value))
  out_row <- c(out_row,which(temp<(mean_value-5*sd_value)))
  data[lipids_list[i]][out_row,] <- NA
  print(length(out_row))
}

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
high_risk <- readRDS('high_risk_ind.RDS')

# name of pathways
files <- list.files('.','subgroup')
files <- files[grep('RDS',files)]
#files <- files[-c(1,2)]
psprs <- str_replace_all(files,'subgroup_','')
psprs <- str_replace_all(psprs,'.RDS','')

# calculate interactions for specific traits and specific psPRS
getInter <- function(trait,prs_name){ 
    temp <- subset(data,select=c('eid','follow.up','CAD','age_recruit',
				 'sex','PC1','PC2','PC3','PC4',
				 trait))
    colnames(temp)[10] <- 'Phenotype' 
    temp <- temp[is.na(temp$Phenotype)==F,]

    # high-risk subjects
    highrisk <- readRDS('high_risk_ind.RDS')
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
    formu <- formula('Surv(follow.up,CAD) ~ index+age_recruit+sex+PC1+PC2+PC3+PC4')
    cox <- summary(coxph(formu,data=temp))
    out <- cox$coefficients[1:3,c(1:3,5)]
    out <- rbind(c(0,1,0,0),out)
    add <- cbind(table,out)

    return(add)
}


phenotype <- c('tc_bi','ldl_bi','tg_bi','sleep_bi',
               'FVC_bi','FEV1_bi','PEF_bi',
               'dbp_bi','sbp_bi',
               'cur_tobacco_smk_1','cur_tobacco_smk_1',
               'smoking_current','smoking_current')
prs_names <- c('Lipids','Lipids','Lipids','BP',
               'BP','BP','BP',
               'BP','BP',
               'Respiratory_syste','T2D',
               'Respiratory_syste','T2D')
res <- list()
for(i in 1:length(phenotype)){
  res[[i]] <- getInter(phenotype[i],prs_names[i])
}
names(res) <- paste0(phenotype,'*',prs_names)

saveRDS(res,'inter_tranas_select_update_vis.RDS')

