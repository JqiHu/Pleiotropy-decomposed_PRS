# ARR for subgroups
library('data.table')
library('stringr')
rm(list=ls())

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
high_risk <- readRDS('high_risk_ind.RDS')

phe_data <- data[data$eid %in% high_risk$IID,]
data <- phe_data # high-risk subjects only

# name of pathways
files <- list.files('.','subgroup')
files <- files[grep('RDS',files)]
psprs <- str_replace_all(files,'subgroup_','')
psprs <- str_replace_all(psprs,'.RDS','')

# ARR
## and define binary traits based on median value
getARR <- function(pheno,prs){ # get ARR for 1 pheno in 1 subgroups 
  temp <- subset(data,select=c('eid','CAD',pheno))
  colnames(temp)[3] <- 'phe'
  print(table(temp$phe))
  # for specific subgroup
  i=grep(prs,psprs)
    subgroup <- readRDS(files[i])
#    remain <- readRDS('remain_ind.RDS')
#    temp$sub[temp$eid %in% remain] <- 0
    temp$sub <- 0
    temp$sub[temp$eid %in% subgroup$IID] <- 1
    print(table(temp$sub))
    temp <- na.omit(temp)
    # number of CAD patients/total number of sub=1 and phe=1
    add1 <- c(psprs[i],1,sum(temp$CAD[temp$phe==1&temp$sub==1]),sum(temp$phe==1&temp$sub==1))
    # number of CAD patients/total number of sub=1 and phe=0
    add2 <- c(psprs[i],0,sum(temp$CAD[temp$phe==0&temp$sub==1]),sum(temp$phe==0&temp$sub==1))
    # number of CAD patients/total number of sub=0 and phe=1
    add3 <- c(paste0(psprs[i],'_remain'),1,sum(temp$CAD[temp$phe==1&temp$sub==0]),sum(temp$phe==1&temp$sub==0))
    # number of CAD patients/total number of sub=0 and phe=0
    add4 <- c(paste0(psprs[i],'_remain'),0,sum(temp$CAD[temp$phe==0&temp$sub==0]),sum(temp$phe==0&temp$sub==0))

    add <- rbind.data.frame(add1,add2,add3,add4)
    colnames(add) <- c('subgroup','phenotype','num_cad','num_tot')
    add[,-c(1:2)] <- apply(add[-c(1:2)],2,as.numeric)
    add$ar <- add$num_cad/add$num_tot # absolute risk = incidence
    # risk difference across different phenotype groups
    arr <- add$ar[add$phenotype==1]-add$ar[add$phenotype==0]
    se <- sqrt(add$ar[add$phenotype==1]*(1-add$ar[add$phenotype==1])/add$num_tot[add$phenotype==1]+
                add$ar[add$phenotype==0]*(1-add$ar[add$phenotype==0])/add$num_tot[add$phenotype==0])
    add$arr <- arr[c(1,1,2,2)]
    add$se <- se[c(1,1,2,2)]

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


arr_out <- c()
for(i in 1:length(phenotype)){
  arr_add <- cbind(phenotype[i],getARR(phenotype[i],prs_names[i]))
  arr_out <- rbind(arr_out,arr_add)
}
colnames(arr_out) <- c('trait','subgroup','phenotype',
  'num_cad','num_tot','ar','arr','se')

# output arr table
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/interaction_new_pipeline/data')
write.table(arr_out,paste0('sub_arr_trans_selected_',thresh,'.tsv'),
  row.names=F,col.names=T,quote=F,sep='\t')
