# deal with overlapping SNPs between pleiotropy-decomposed regions.
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_beta_selection_V2')
library(stringr)
library(data.table)

trait43 <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/trait_43/trait_43_info.csv') # information for 43 traits
trait43$id <- paste(trait43$Group,trait43$Code,sep='_') # for matching
trait43 <- trait43[order(trait43$id),]
file <- list.files('.')
file <- sort(file)
file <- file[c(1:28,31:32,29:30,33:43)] # match order to trait43

# combine SNPs for the same path
res <- list() # list for all correlation results
path <- unique(trait43$Cluster)
for(k in 1:length(path)){
  index <- which(trait43$Cluster==path[k]) # row numbers for certain pleiotropy-decomposed region
  corr <- c() # file for beta of one pleiotropy-decomposed region
  print(path[k])
  print(file[index])
  for(i in 1:length(index)){
    temp <- fread(file[index[i]]) # beta file with region correlation for specific trait
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
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/prs_beta_selection_final_V2_absolute/') # path for output
sum <- c() # number of unique SNPs in each pathway
for(i in 1:length(path)){
  out <- rmOverlap(res[[i]],i)
  write.table(out,paste0(str_replace_all(path[i],' ','_'),'_',nrow(out),'.txt'),
		row.names=F,col.names=T,quote=F)
  sum <- c(sum,nrow(out))
}
sum <- cbind.data.frame(path,sum)
colnames(sum) <- c('pathway','number_of_SNP')
sum <- sum[order(sum$number_of_SNP,decreasing = T),]
write.csv(sum,
	'/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/res/table/ta1_1_V2_absolute.csv',
	row.names=F)

