# Generate SNP subsets based on simulated local genetic correlations
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/pd_prs/data/snp_subset')
library(data.table)

# Simulated SNP betas
gwas <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_beta/simulated_beta_by_region.txt')

# Simulated local genetic covariance
covar <- readRDS('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_beta/simu_local_covar.RDS')
# Information on simulated causal regions and correlated regions
snp_info <- fread('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_beta/snp_info.txt')

#### Partition SNPs into three subgroups and 1 NS-group
#### Group 1 for trait 2
group1 <- snp_info[snp_info$Correlated_Trait2==T,]
#### Group 2 for trait 3
group2 <- snp_info[snp_info$Correlated_Trait3==T,]
#### Group 3 for trait 4
group3 <- snp_info[snp_info$Correlated_Trait4==T,]
#### Group 4 for non-specific SNPs
group4 <- snp_info[snp_info$marker.ID %in% c(group1$marker.ID,group2$marker.ID,group3$marker.ID)==F,]

##### Check overlap
### Function
getOverlapG <- function(t1,t2){ # t1 t2 are the trait numbers
  g1 <- get(paste0('group',t1-1))
  g2 <- get(paste0('group',t2-1))
  overlap <- merge(g1,g2,by='marker.ID')
  if(nrow(overlap)==0){
    return(NULL)
  }
  overlap.region <- unique(overlap$block.x)

  overlap.assign <- c()
  for(r in overlap.region){
    snps.region <- overlap$marker.ID[overlap$block.x==r] 
    local.covar.matrix <- covar[[r]]
    # Convert to matrix if necessary
#    if (class(local.covar.matrix) != "matrix") {
#      local.covar.matrix <- as.matrix(local.covar.matrix)
#    }
    local.covar <- local.covar.matrix[1,c(t1,t2)]
    
    ### assign SNPs to group with larger covar
    larger.group <- c(t1,t2)[which.max(abs(local.covar))]
    add <- cbind(snps.region,larger.group)
    overlap.assign <- rbind(overlap.assign,add)
  }

  colnames(overlap.assign) <- c('snp_overlap','trait')
  return(overlap.assign)
}
### get all overlap SNPs
overlaps <- c()
for(i in 2:3){
  for(j in (i+1):4){
    add <- getOverlapG(i,j)
    overlaps <- rbind(overlaps,add)
  }
}
overlaps <- as.data.frame(overlaps)
overlaps$trait <- as.numeric(overlaps$trait)
### some variants shared by multiple groups
overlaps <- overlaps[order(overlaps$snp_overlap,overlaps$trait,decreasing=F),]
overlaps <- overlaps[duplicated(overlaps$snp_overlap)==F,]
### Remove overlaps from current list
for(t in 2:4){
  tmp <- get(paste0('group',t-1))
  snp.rm <- overlaps$snp_overlap[overlaps$trait!=t]
  snp.rm <- intersect(tmp$marker.ID,snp.rm)

  if(length(snp.rm)==0){
    out <- tmp
  }else{
    out <- tmp[-which(tmp$marker.ID %in% snp.rm),]
  }
  assign(paste0('group',t-1),out)
  print(dim(out))
}

# output all SNP subsets
gwas <- as.data.frame(gwas)
for(i in 1:4){
  tmp <- get(paste0('group',i))
  out <- gwas[gwas$marker.ID %in% tmp$marker.ID,]
  print(dim(out))
  write.table(out,
	      paste0('simulated_beta_group',i,'.txt'),
	      row.names=F,sep='\t',quote=F)
}
