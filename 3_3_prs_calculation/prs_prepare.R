# aggregate 9 PD-PRSs
library(data.table)
library(stringr)

setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/psPRS_V2_absolute/')# path to PD-PRSs

prs_file <- list.files('.','profile') # file for PD-PRSs
name <- c('Basic_condition','BP','CVD',
	  'Immune_system','Lipids','Others',
	  'Obesity','Respiratory_system','T2D') # name for PD-PRS file

# combine all PD-PRSs in to one file
ps_prs <- subset(fread(prs_file[1]),select=c('IID','SCORE'))
colnames(ps_prs)[2] <- name[1] # change column name into the name of pleiotropy-decomposed regions
for (i in 2:length(prs_file)) {
  prs <- subset(fread(prs_file[i]),select=c('IID','SCORE'))
  colnames(prs)[2] <- name[i]
  ps_prs <- merge(ps_prs,prs,by = 'IID')
}

# Input CAD study subjects
cad <- fread('../phenotype/pheno_cad.tsv')

ps_prs <- ps_prs[ps_prs$IID %in% cad$eid[is.na(cad$CAD)==F],]
ps_prs_out <- cbind.data.frame(ps_prs[,1],
	        apply(ps_prs[,-1],2,
		  function(x){scale(x,scale=T,center=T)})) # scale psPRSs


# output prs file
write.table(ps_prs_out,
	    'ps_prs.tsv',
	   row.names=F,col.names=T,quote=F,sep='\t')
