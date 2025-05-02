# Analyze associations between PD-PRSs and phenotypes of interest
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/clustering_sensitivity/hierachical_clustering/data')
library(data.table)
library(stringr)

### import number of clusters
args <- commandArgs(trailingOnly=T)
k <- args[1]
setwd(paste0('./cluster_',k,'/prs_analysis'))

# PD-PRS
prs <- fread('../PD_PRS.txt',data.table=F)
ps_names <- colnames(prs)[-1]
### CAD prs
prs_cad <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS/cad_prs.tsv')
prs_cad <- prs_cad[,c(2,6)] # subtract iid and score
colnames(prs_cad) <- c('eid','PRS_CAD')

prs <- merge(prs_cad,prs,by='eid') # merge prs

#### extract study subjects
study_subject <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/phenotype/pheno_all.tsv',
		       select=c('eid'))
prs.sub <- prs[prs$eid %in% study_subject$eid,]
##### standardize PRS
prs.sub <- as.data.frame(prs.sub)
prs.sub[,-1] <- apply(prs.sub[,-1],2,scale)

## top 5% of PD-PRSs in general population
cutoffs <- apply(prs.sub[,-c(1,2)],2,
                 function(x){quantile(x,0.95)})

## identify top 5% CAD PRS
prs_high <- prs.sub[which(prs.sub$PRS_CAD >
                        quantile(prs.sub$PRS_CAD,0.95)),] # extract top 5% risk subjects for CAD
prs_high <- apply(prs_high,2,as.numeric)
prs_high <- as.data.frame(prs_high)

saveRDS(prs_high,'high_risk_ind.RDS')

# assign subjects into subgroups by
# selecting the subjects in top 5% PD-PRSs among individuals with 5% top CAD PRS
subgroup <- c()
for(i in 1:length(ps_names)){
  temp <- subset(prs_high,select=c('eid',ps_names[i]))
  colnames(temp)[2] <- 'psprs'
  add_ind <- temp$eid[temp$psprs>=cutoffs[i]]
  # output subgroup with PRS
  out <- prs_high[prs_high$eid %in% add_ind,]
  print(nrow(out))
  saveRDS(out,paste0('subgroup_',
                     str_replace_all(ps_names[i],' ','_'),'.RDS'))
  out <- cbind.data.frame(ps_names[i],out)
  subgroup <- rbind(subgroup,out)
}
colnames(subgroup)[1] <- 'subgroup'

## calculate proportions of each PD-PRS on overall
num_col <- ncol(prs_high)
prs_high <- as.data.frame(prs_high)
sum_psprs <- apply(abs(prs_high[,3:num_col]),1,sum) # sum of 9 psPRSs
mean_psprs <- apply(prs_high[,3:num_col],2,mean) # mean of each psPRS

stack <- cbind(ps_names,abs(mean_psprs)/sum(abs(mean_psprs))
                        *mean(prs_high$PRS_CAD)) # mean proportions
# absolute value of CAD PRS explained by PD-PRS
prop <- apply(abs(prs_high[,3:num_col]),2,
                function(x){(abs(x)/sum_psprs)*prs_high$PRS_CAD})

# Combine with subgroup information
wide <- cbind.data.frame(subgroup$subgroup,
                         prs_high$eid[match(subgroup$eid,prs_high$eid)],
                         prop[match(subgroup$eid,prs_high$eid),])
colnames(wide)[1:2] <- c('subgroup','eid')

long <- melt(setDT(wide),id.vars=c('eid','subgroup'),
                variable.name='psprs') # long format for plot

write.table(long,
            'psprs_prop_long.tsv',
            row.names=F,quote=F,sep='\t')

