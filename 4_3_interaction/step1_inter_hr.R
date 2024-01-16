# calculate hr for interactions 
# between psPRSs and transformed selected traits
library(survival)
library(lmtest)
library(data.table)

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

##### Function to do this for specific binary traits
binary_list <- c("insomnia_2","insomnia_3",
		 "cur_tobacco_smk_1","cur_tobacco_smk_2",
		 "past_tobacco_smk_1","past_tobacco_smk_2","past_tobacco_smk_3",
		 'smoking_current',"smoking_ever")

# calculate interactions in all subjects
getInter <- function(trait){ # calculate interactions between trait and 9 psPRSs
  res <- c()
  for(i in 1:length(ps_names)){
    formu <- as.formula(paste0('Surv(follow.up,CAD) ~ ',
		ps_names[i],'*',trait,'+ age_recruit + sex + PC1 + PC2 + PC3 + PC4',
		collapse=" "))
    formu2 <- as.formula(paste0('Surv(follow.up,CAD) ~ ',
                ps_names[i],'+',trait,'+ age_recruit + sex + PC1 + PC2 + PC3 + PC4',
                collapse=" "))
    cox <- coxph(formu,data=data)
    cox2 <- coxph(formu2,data=data)
    ## Likelihood ratio test
    test <- lrtest(cox,cox2)
    add <- test$`Pr(>Chisq)`[2]
    res <- c(res,add)  
  }
  return(res)
}

phenotype <- c(lipids_list,binary_list)
inter_out <- data.frame()
for(i in 1:length(phenotype)){
  inter_out <- rbind(inter_out,
	getInter(phenotype[i]))
  print(phenotype[i])
}
inter_out <- cbind.data.frame(phenotype,inter_out)
colnames(inter_out) <- c('trait',ps_names)


# output
write.csv(inter_out,
	  '/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_analysis_V2_absolute//interaction_selected.csv',
  	  row.names=T)

