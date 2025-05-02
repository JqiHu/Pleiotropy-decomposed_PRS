# classify top 5% CAD PRS into 9 groups
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_analysis_V2_absolute/')
library('data.table')
library('stringr')
rm(list=ls())

## input PD-PRSs and CAD PRS
psprs <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS_V2_absolute//ps_prs.tsv')
ps_names <- c('Basic Condition','BP','CVD','Immune System',
              'Lipids','Others','Obesity','Respiratory system','T2D') # name for the pleiotropy-decomposed regions
colnames(psprs)[-1] <- ps_names # !!!! sequence!!!!!

prs_cad <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS/cad_prs.tsv')
prs_cad <- prs_cad[,c(2,6)] # subtract iid and score
colnames(prs_cad)[2] <- 'CAD'
# prs_cad$CAD <- scale(prs_cad$CAD,
#                center=T, scale=T) # scale

prs <- merge(prs_cad,psprs,by='IID') # merge prs

## top 5% of PD-PRSs in general population
cutoffs <- apply(prs[,-c(1,2)],2,
		 function(x){quantile(x,0.95)})

## identify top 5% CAD PRS
prs_high <- prs[which(prs$CAD > 
                        quantile(prs$CAD,0.95)),] # extract top 5% risk subjects for CAD
prs_high <- apply(prs_high,2,as.numeric)
prs_high <- as.data.frame(prs_high)

saveRDS(prs_high,'high_risk_ind.RDS')

# assign subjects into subgroups by 
# selecting the subjects in top 5% PD-PRSs among individuals with 5% top CAD PRS
subgroup <- c()
for(i in 1:length(ps_names)){
  temp <- subset(prs_high,select=c('IID',ps_names[i])) 
  colnames(temp)[2] <- 'psprs'
  add_ind <- temp$IID[temp$psprs>=cutoffs[i]]
  # output subgroup with PRS
  out <- prs_high[prs_high$IID %in% add_ind,]
  print(nrow(out))
  saveRDS(out,paste0('subgroup_',
		     str_replace_all(ps_names[i],' ','_'),'.RDS'))
  out <- cbind.data.frame(ps_names[i],out)
  subgroup <- rbind(subgroup,out)
}
colnames(subgroup)[1] <- 'subgroup'

## calculate proportions of each PD-PRS on overall
sum_psprs <- apply(abs(prs_high[,3:11]),1,sum) # sum of 9 psPRSs
mean_psprs <- apply(prs_high[,3:11],2,mean) # mean of each psPRS

stack <- cbind(ps_names,abs(mean_psprs)/sum(abs(mean_psprs))
                        *mean(prs_high$CAD)) # mean proportions
# absolute value of CAD PRS explained by PD-PRS
prop <- apply(abs(prs_high[,3:11]),2,
                function(x){(abs(x)/sum_psprs)*prs_high$CAD}) 

# Combine with subgroup information
wide <- cbind.data.frame(subgroup$subgroup,
			 prs_high$IID[match(subgroup$IID,prs_high$IID)],
			 prop[match(subgroup$IID,prs_high$IID),])
colnames(wide)[1:2] <- c('subgroup','IID')
#wide <- cbind.data.frame(prs_high$IID,prop,
#			 subgroup$subgroup[match(prs_high$IID,subgroup$IID)])
#colnames(wide)[c(1,11)] <- c('IID','subgroup')

long <- melt(setDT(wide),id.vars=c('IID','subgroup'),
                variable.name='psprs') # long format for plot 

write.table(long,
	    'psprs_prop_long.tsv',
            row.names=F,quote=F,sep='\t')
