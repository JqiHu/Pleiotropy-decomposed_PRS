# compare level of specific traits
# across genetic subgroups
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_analysis_V2_absolute/')
library('data.table')
library('stringr')

# phenotype information
phe_data <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/pheno_all.tsv') # phenotype for all
phe_name <- colnames(phe_data)[-1]

# input high-risk individuals
high_risk <- readRDS('high_risk_ind.RDS')

phe_data <- phe_data[phe_data$eid %in% high_risk$IID,]

# Add fat percentages
#update <- readRDS('../phenotype/ukbb_modifiable_final.RDS') # imputed
update <- readRDS('../phenotype/ukbb_all_modifiable.RDS') # original
update <- as.data.frame(update)
update <- update[,-which(duplicated(colnames(update))==T)]
update <- merge(phe_data[,1:2],update,by='eid')

update$arm_fat_percentage <- (update$`23111`+update$`23115`)/2
update$leg_fat_percentage <- (update$`23119`+update$`23123`)/2
update_add <- subset(update,select=c('eid','23099','23127','arm_fat_percentage','leg_fat_percentage'))
colnames(update_add)[2:3] <- c('body_fat_percentage','trunk_fat_percentage')

phe_data <- merge(phe_data,update_add,by='eid',all.x=T)
phe_name <- colnames(phe_data)[-1]


## function to distinguish categorical traits from continuous ones
getIndex <- function(pheno){ # input phenotype name
  if(pheno %in% c('lifestyle_comb','alcohol','physical_activity','income',
                  'smoking','lifestyle_ideal_count','lifestyle_poor_count')){
    return('quali_multi')
  }else{
    res <- ifelse(nrow(unique(na.omit(subset(data,select=pheno))))<3,'quali','quanti')
    return(unique(res))
  }
}

# name of pathways
files <- list.files('.','subgroup')
files <- files[grep('RDS',files)]
psprs <- str_replace_all(files,'subgroup_','')
psprs <- str_replace_all(psprs,'.RDS','')

# study subjects
cad <- fread('../phenotype/pheno_cad.tsv')
iid <- cad$eid[is.na(cad$CAD)==F]

# Covariates
pc <- fread('/gpfs/gibbs/pi/zhao/zhao-data/yy496/ukbb_pheno/MR_time/adjust/adjust.csv') # principle components
pc <- subset(pc,select=c('eid',
                  paste0('PC',1:4)))
data <- merge(pc,phe_data,by='eid')

# Create a new variable for smoking
# Compare ever to never
data$smoking1[data$smoking==1] <- 1
data$smoking1[data$smoking==0] <- 0

data$smoking2[data$smoking==2] <- 1
data$smoking2[data$smoking==0] <- 0

# function to check the distribution 
# of a trait comparing subgroup and the remaining
getTest <- function(phe,estimate){ # estimate can be either p-val or beta
  temp <- subset(data,select=c('eid',paste0('PC',1:4),'sex','age_recruit',phe))
  colnames(temp)[8] <- 'Phenotype'
  # compare 9 subgroups separately
  out <- c()
  for(i in 1:9){
    subgroup <- readRDS(files[i])
    temp$index <- 0
    temp$index[temp$eid %in% subgroup$IID] <- 1 
    print(table(temp$index))
    formu <- formula('Phenotype ~ index+age_recruit+sex+PC1+PC2+PC3+PC4')
    if(getIndex(phe)=='quanti'){fit <- lm(formu,temp)}
    if(getIndex(phe)=='quali'){fit <- glm(formu,temp,family='binomial')} 
    res <- summary(fit)
    add <- ifelse(estimate=='p-val',res$coefficients[2,4],
		  res$coefficients[2,1])
    out <- c(out,add)
  }
  return(out)
}

# phenotype of interest
phenos <- c('bmi','waist_circumference','body_fat_percentage','trunk_fat_percentage','arm_fat_percentage','leg_fat_percentage',
	    'SBP','DBP','hypertation',
	    'DIA','Glucose',
	    'smoking1','smoking2','smoking_current','smoking_ever',
	    'Apolipoprotein_A',"Apolipoprotein_B","Cholesterol",'HDL_cholesterol',"LDL_direct","Lipoprotein_A","Triglycerides",
	    "Creatinine_enzymatic_in_urine","Urea","Creatinine","Cystatin_C","Urate")

res <- c()
for(i in 1:length(phenos)){
  add <- getTest(phenos[i],'p-val')
  res <- rbind(res,add)
}
res <- cbind.data.frame(phenos,res)
colnames(res) <- c('triat',psprs)

write.csv(res,'trait_compare_pval.csv',row.names=F)



