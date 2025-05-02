# Combine selected SNPs to PD clusters defined 
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/clustering_sensitivity/hierachical_clustering/data')
library(data.table)
library(stringr)

### import number of clusters
args <- commandArgs(trailingOnly=T)
k <- args[1]
setwd(paste0('./cluster_',k,'/prs_beta_selection'))

## info of clustering
cluster <- fread('../cluster_info.csv',data.table=F)
cluster$prefix <- str_replace_all(cluster$Feature,' ','_')
######## change Inflammation to inflammation
cluster$prefix <- str_replace_all(cluster$prefix,'Inflamm','inflamm')
### files of selected SNPs for each trait
rootpath='/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_beta_selection_V2'
file <- list.files(rootpath)
###### get prefix
file.list <- strsplit(file,'_beta_')
file.prefix <- unlist(lapply(file.list,function(x){x[1]}))
######## match file names
cluster$file <- file[match(cluster$prefix,file.prefix)]

# combine SNPs for the same path
res <- list() # list for all correlation results
path <- unique(cluster$Cluster)
for(k in 1:length(path)){
  index <- which(cluster$Cluster==path[k]) # row numbers for certain pleiotropy-decomposed region
  corr <- c() # file for beta of one pleiotropy-decomposed region
  for(i in index){
    temp <- fread(paste0(rootpath,'/',cluster$file[i])) # beta file with region correlation for specific trait
    colnames(temp)[10] <- 'index_p_unique' # this column is the p-value from SUPERGNOVA of local genetic covariance
    colnames(temp)[11] <- 'index_rho_unique' # genetic covariance
    corr <- rbind(corr,temp)
  }
  corr <- corr[order(abs(corr$index_rho_unique),decreasing=T),] # rank SNPs by decreasing the absolute value of genetic covariance from SUPERGNOVA
  # remove duplicated info with larger pval/smaller rho
  corr2 <- corr[which(duplicated(corr$sid)==F),]
  res[[k]] <- corr2
  print(nrow(res[[k]]))
}

# construct functions to deal with overlapping SNPs
rmOverlap <- function(beta,num){ # beta is the file to be compared with others and num is the index in `path`
  n <- 1:length(path)
  n <- n[-num] # path number except for the compared one
  for(i in n){
    overlap <- merge(subset(beta,select=c('sid','index_rho_unique')), # merge by SNP rsid to find overlapping
                     subset(res[[i]],select=c('sid','index_rho_unique')),by='sid')
    # remove SNPs from less significant/smaller rho path
    overlap$rm <- ifelse(abs(overlap$index_rho_unique.x)<abs(overlap$index_rho_unique.y),1,2)
    rm_sid <- overlap[overlap$rm==1,]$sid
    keep <- setdiff(beta$sid,rm_sid) # SNPs to keep
    beta <- beta[beta$sid %in% keep,]
  }
  return(beta)
}

# compare one pleiotropy-decomposed region to others one by one
sum <- c() # number of unique SNPs in each pathway
for(i in 1:length(path)){
  out <- rmOverlap(res[[i]],i)
  write.table(out,paste0('Cluster',path[i],'.txt'),
                row.names=F,col.names=T,quote=F)
  sum <- c(sum,nrow(out))
}
sum <- cbind.data.frame(path,sum)
colnames(sum) <- c('pathway','number_of_SNP')
sum <- sum[order(sum$number_of_SNP,decreasing = T),]

write.csv(sum,'table_num_SNP_cluster.csv',row.names=F)

