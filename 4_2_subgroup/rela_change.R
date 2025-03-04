# relative change on traits
# between subgroups
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_analysis_V2_absolute/')
library('data.table')
library('stringr')
rm(list=ls())

# phenotype information
phe_data <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/pheno_all.tsv') # phenotype for all

### biological clusters for phenotypes
annot <- as.data.frame(c(rep('Demographic',5),rep('Obesity',2),rep('Respiratory',3),rep('Lifestyle',6),
                         rep('CVD',3),rep('Family history',4),rep('T2D',2),rep('Lipids',8),
                         rep('Protein',2),rep('Sex hormone',3),rep('Urinary system',5),rep('Metabolism',5),
                         rep('Liver',5),rep('Immune system',2))) ## biological annotations for phenotypes
rownames(annot) <- colnames(phe_data)[-1]
phe_name <- colnames(phe_data)[-1]

# input high-risk individuals
high_risk <- readRDS('high_risk_ind.RDS')

phe_data <- phe_data[phe_data$eid %in% high_risk$IID,]

# Add fat percentages
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
    res <- ifelse(nrow(unique(na.omit(subset(phe_data,select=pheno))))<3,'quali','quanti')
    return(unique(res))
  }
}

# name of pleiotropy-decomposed regions
files <- list.files('.','subgroup')
files <- files[grep('RDS',files)]
files <- files[-1]
psprs <- str_replace_all(files,'subgroup_','') 
psprs <- str_replace_all(psprs,'.RDS','')

## function to calculate relative change on one phenotype
## among top 5% CAD PRS individuals
getRelachan <- function(pheno){
  phe <- subset(phe_data,select=c('eid',pheno))
  colnames(phe)[2] <- 'trait'
  out <- c()
  for(prs in psprs){
    # input subgroup data
    temp <- readRDS(paste0('subgroup_',prs,'.RDS'))
    index <- ifelse(phe$eid %in% temp$IID,1,0) # 1 represents target subgroup
    if(getIndex(pheno)=='quali_multi'){ # only count for highest level
      sub <- length(which(na.omit(phe$trait[index==1])== # for categorical variable, only count for the max level
				max(na.omit(phe$trait[index==1]))))/length(na.omit(phe$trait[index==1])) # proportion of max-level who are in the certain PD-prs subgroup
      rem <- length(which(na.omit(phe$trait[index==0])==
				max(na.omit(phe$trait[index==0]))))/length(na.omit(phe$trait[index==0])) # proportion of max-level who are in the the rest PD-prs subgroups
      res <- (sub-rem)/sd(na.omit(phe$trait))   
    } else{
      res <- (mean(na.omit(phe$trait[index==1])) # proportion of level 1 or mean value in the certain PD-prs subgroup
		-mean(na.omit(phe$trait[index==0])))/sd(na.omit(phe$trait)) # 
    }
    out <- rbind(out,cbind(prs,res))
  }
  return(out)
}

## relative change
#target <- c('Lipids','Obesity','Urinary system') # target phenotype clusters
#pheno_name <- phe_name[which(annot[,1] %in% target)]
pheno_name <- phe_name

rela_chan <- data.frame()
for(phe in pheno_name){
  add <- cbind(phe,getRelachan(phe)) 
  rela_chan <- rbind(rela_chan,add)
}
colnames(rela_chan) <- c('pheno','psPRS','change')

## output relative change
write.table(rela_chan,
	    'rela_change.tsv',
	   row.names=F,quote=F,sep='\t')




