# Simulate beta as GWASs
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/simu_data/data/simu_beta')
library(data.table)
library(MASS)
library(bigsnpr)
library(Matrix)
#remotes::install_github("QuantGen/BEDMatrix")
library(BEDMatrix)

# Read Genotype data from 1kg hapmap3 chr22
plink_file <- "../hapmap3/1000G_mac5eur_hapmap3_chr22_forUKB"
##### Check if already exist files
index <- list.files('../hapmap3','forUKB.bk')
if(length(index)>0){
  unlink('../hapmap3/1000G_mac5eur_hapmap3_chr22_forUKB.bk')
  unlink('../hapmap3/1000G_mac5eur_hapmap3_chr22_forUKB.rds')
}

# Load PLINK files using bigsnpr
snp_data <- snp_attach(snp_readBed(paste0(plink_file, ".bed")))
# Access genotype matrix (individuals x SNPs)
genotypes <- snp_data$genotypes  # A "FBM.code256" object
snp_info <- snp_data$map         # SNP metadata: chromosome, position, etc.
genotypes.col <- snp_info$marker.ID

# Load LD block definitions
ref <- fread('./LDBlock/eur_chr22.bed',data.table=F)
print(nrow(ref)) # 47 regions

# Assign SNPs to each region
snp_info2 <- c()
for(i in 1:nrow(ref)){
  tmp <- snp_info[snp_info$physical.pos<=ref$stop[i] &
		  snp_info$physical.pos>=ref$start[i],]
  add <- cbind.data.frame(tmp,i)
  colnames(add)[ncol(add)] <- 'block'
  snp_info2 <- rbind(snp_info2,add) 
}

# Some missing
missing_snps <- snp_info[which(snp_info$marker.ID %in% snp_info2$marker.ID==F),]
missing_snps$block <- 0
# Assign missing SNPs to the nearest block
for (i in seq_len(nrow(missing_snps))) {
  # distance to blocks
  dist1 <- abs(missing_snps$physical.pos[i]-ref$start)
  dist2 <- abs(missing_snps$physical.pos[i]-ref$stop)
  distances <- apply(cbind.data.frame(dist1,dist2),1,min)
  nearest_block <- which.min(distances)

  # Assign SNP to the nearest block
  missing_snps$block[i] <- nearest_block
}
snp_info <- rbind(snp_info2,missing_snps)


# -------------------
# Set up parameters
set.seed(123)             # For reproducibility
num_traits <- 4
num_snps <- nrow(snp_info)
n_regions <- max(snp_info$block)
num_causal_blocks <- ceiling(n_regions*0.2) # in total, 47 blocks. Assume 20% causal

# Assign causal blocks for each trait
# Step 1: Ensure at least one overlapping block between Trait 1 and Traits 2-4
#overlap_blocks <- sample(1:n_regions, 3, replace = FALSE)  # Select 3 overlapping blocks

# Step 2: Assign remaining causal blocks for each trait
#causal_blocks <- list()
#for (trait in 1:num_traits) {
#  if (trait == 1) {
#    # Trait 1: Ensure overlapping blocks are included
#    causal_blocks[[trait]] <- c(overlap_blocks, sample(setdiff(1:n_regions, overlap_blocks), 
#                                                       num_causal_blocks - 3, replace = FALSE))
#  } else {
#    # Traits 2-4: Ensure each overlaps with Trait 1
#    causal_blocks[[trait]] <- c(overlap_blocks[trait-1],  # Overlap with Trait 1
#                                sample(setdiff(1:n_regions, overlap_blocks[trait-1]), 
#                                       num_causal_blocks - 1, replace = FALSE))
#  }
#}

## Randomly assign blocks
##### and restrict the number of shared causal blocks to be 3!!!!!
causal_blocks <- list()
n_overlap <- 3 # number of shared causal regions
for(trait in 1:num_traits) {
  if(trait==1){
    causal_blocks[[trait]] <- sample(1:n_regions,num_causal_blocks,replace=F)
  }else{
    overlap <- causal_blocks[[1]][seq(trait,10,3)]
    #overlap <- sample(causal_blocks[[1]],n_overlap)
    causal_blocks[[trait]] <- c(overlap,
				sample(setdiff(1:n_regions,causal_blocks[[1]]),
				       num_causal_blocks - n_overlap,replace=F))
  }
   
}

### check overlap
print(intersect(causal_blocks[[1]],causal_blocks[[2]]))
print(intersect(causal_blocks[[1]],causal_blocks[[3]]))
print(intersect(causal_blocks[[1]],causal_blocks[[4]]))

# Mark causal SNPs
##### Allow overlapping
for (trait in 1:num_traits) {
  add <- rep(0,nrow(snp_info))
  add[snp_info$block %in% causal_blocks[[trait]]] <- 1

  snp_info <- cbind.data.frame(snp_info,add)
  colnames(snp_info)[ncol(snp_info)] <- paste0("Causal_Trait", trait)
}
print(table(snp_info$Causal_Trait1,snp_info$Causal_Trait2))
print(table(snp_info$Causal_Trait1,snp_info$Causal_Trait3))
print(table(snp_info$Causal_Trait1,snp_info$Causal_Trait4))

# Heritabilities across simulated diseases/traits
### CAD: https://www.nature.com/articles/s41591-022-01891-3#Sec2
h2 <- c(0.244,0.2,0.2,0.2)
causal_mean_effect <- rep(1,4) # Non-zero mean for causal SNPs
###### Local heritabilities
######## Assumption: causal regions should be assigned with larger local heritability
causal_weight <- 0.8  # 80% of heritability is from causal regions
non_causal_weight <- 1 - causal_weight

scaled_local_h2 <- matrix(0, nrow = n_regions, ncol = num_traits)

for (trait in 1:num_traits) {
  # Separate causal and non-causal regions
  causal_regions <- which(1:n_regions %in% causal_blocks[[trait]])
  non_causal_regions <- setdiff(1:n_regions, causal_regions)

  # Randomly assign heritability within causal and non-causal regions
  raw_causal_h2 <- runif(length(causal_regions), 0.02, 0.1)
  raw_non_causal_h2 <- runif(length(non_causal_regions), 0.001, 0.02)

  # Scale heritabilities to match overall heritability
  total_causal_h2 <- sum(raw_causal_h2)
  total_non_causal_h2 <- sum(raw_non_causal_h2)

  scaled_causal_h2 <- raw_causal_h2 * (h2[trait] * causal_weight / total_causal_h2)
  scaled_non_causal_h2 <- raw_non_causal_h2 * (h2[trait] * non_causal_weight / total_non_causal_h2)

  # Assign scaled heritabilities back to regions
  scaled_local_h2[causal_regions, trait] <- scaled_causal_h2
  scaled_local_h2[non_causal_regions, trait] <- scaled_non_causal_h2
}

# Confirm the sum of local heritabilities equals the overall heritability
total_local_h2 <- colSums(scaled_local_h2)
print(total_local_h2)  # Should match `h2` approximately

###### local genetic correlations
# Assign correlated regions for each trait except Trait 1!!!!!
n_regions <- max(snp_info$block)
num_corr_region <- 10
correlated_regions <- lapply(2:num_traits, function(trait) {
  causal_bs <- intersect(causal_blocks[[1]],
		    	 causal_blocks[[trait]]) # Shared causal blocks for trait 1 and the target trait;
  causal_b2 <- setdiff(causal_blocks[[trait]],
		       causal_bs) # nonshared
  non_causal_blocks <- setdiff(1:n_regions, c(causal_bs,causal_b2))

  all_blocks <- c(causal_bs,
		  sample(causal_b2,size=2),
		  sample(non_causal_blocks,size=5))
  unique(all_blocks)
})
names(correlated_regions) <- 2:num_traits

# Mark correlated regions for each trait with TRAIT 1!!!
for (trait in 2:num_traits) {
  add <- rep(F,nrow(snp_info))
  add[snp_info$block %in% unlist(correlated_regions[as.character(trait)])] <- T
  snp_info <- cbind.data.frame(snp_info,add)
  colnames(snp_info)[ncol(snp_info)] <- paste0("Correlated_Trait", trait)
}

##### check simulated correlated regions and causal regions
simulated_data <- snp_info
print(table(simulated_data$Causal_Trait2,simulated_data$Correlated_Trait2))
print(table(simulated_data$Causal_Trait3,simulated_data$Correlated_Trait3))
print(table(simulated_data$Causal_Trait4,simulated_data$Correlated_Trait4))


write.table(snp_info,'snp_info.txt',row.names=F,sep='\t',quote=F)

# Initialize correlation matrices for each region
local_correlations <- list()
snp_info <- as.data.frame(snp_info)
for (block in unique(snp_info$block)) {
  tmp <- snp_info[snp_info$block == block, ]

  # Check which traits are correlated with the region
  is_correlated <- subset(tmp,select=paste0('Correlated_Trait',2:4))

  # Initialize the correlation matrix
  r_g <- matrix(0, ncol = num_traits, nrow = num_traits)  # Start with all zeros

  # Assign correlations for Trait 1
  for (j in 2:num_traits) {
    if (any(is_correlated[, j - 1])) {
      r_g[1, j] <- runif(1, 0.1, 0.5)  # Positive correlation with Trait 1
      r_g[j, 1] <- r_g[1, j]         # Symmetric assignment
    } else {
      r_g[1, j] <- 0  # No correlation for non-correlated traits
      r_g[j, 1] <- 0
    }
  }

  # Self-correlation for all traits
  diag(r_g) <- 1

  # No inter-correlations for Traits 2-4
  r_g[2:num_traits, 2:num_traits] <- diag(num_traits - 1)  # Only diagonal elements remain

  # Store correlation matrix
  local_correlations[[as.character(block)]] <- r_g
}


# Save the local correlation results
saveRDS(local_correlations,'simu_local_corr.RDS')

# Calculate genetic covariance matrices for each region
local_covariances <- list()
for (region in 1:n_regions) {
  h2_diag <- diag(sqrt(scaled_local_h2[region, ]))  # Diagonal matrix of sqrt(local h2)
  local_covariances[[region]] <- h2_diag %*% local_correlations[[region]] %*% h2_diag
}
saveRDS(local_covariances,'simu_local_covar.RDS')

# Initialize storage
simulated_data <- list()

# Simulate beta-values for each block
for (block in unique(snp_info$block)) {
  # Subset SNPs in this block
  num.col <- which(snp_info$block==block)
  block_snps <- snp_info[num.col,]
  block_snps <- as.data.frame(block_snps)
  ### map SNPs to genotype
  num.col2 <- match(block_snps$marker.ID,genotypes.col) 
  block_genotypes <- genotypes[, num.col2]

  # Compute LD matrix for the block
  ld_matrix <- cor(block_genotypes)

  # Get genetic covariance matrix for the region
  genetic_cov <- local_covariances[[block]]

  # Define combined covariance matrix
  combined_cov_matrix <- kronecker(ld_matrix, genetic_cov)
  combined_cov_matrix <- Matrix::nearPD(combined_cov_matrix)$mat

  # Define mean vector for \(\beta\)-values
  num_snps <- nrow(block_snps)
  mu <- rep(0, num_snps * num_traits)  # Default mean = 0
  for (trait in 1:num_traits) {
    causal_indices <- which(block_snps[, paste0('Causal_Trait', trait)] == 1)
    if (length(causal_indices) > 0) {
      # Assign larger mean effect for causal SNPs
      mu[(causal_indices - 1) * num_traits + trait] <- causal_mean_effect[trait] 
    }
  }

  # Simulate beta-values
  beta_block <- mvrnorm(
    n = 1,                            # Single realization
    mu = mu,                          # Adjusted mean vector
    Sigma = combined_cov_matrix       # Combined covariance matrix
  )

  # Reshape beta-values into a matrix (SNPs x Traits)
  beta_matrix <- matrix(beta_block, ncol = num_traits, byrow = T)

  # Store \(\beta\)-values
  block_snps <- as.data.table(block_snps)
  block_snps[, paste0("Beta_Trait", 1:num_traits) := as.data.table(beta_matrix)]
  simulated_data[[block]] <- block_snps

}

# Combine results
simulated_data <- rbindlist(simulated_data)

##### match sequence
new.order <- match(genotypes.col,simulated_data$marker.ID)
simulated.data <- simulated_data[new.order,]

# Check MAF
stats <- bigstatsr::big_colstats(genotypes)
# Extract MAF from allele frequencies
maf <- stats$sum / (2 * nrow(genotypes))  # Total allele count divided by total observations

# Add noise to observed betas
sample_size <- 50000
for (trait in 1:num_traits) {
  simulated.data[, paste0("Observed_Beta_Trait", trait) :=
    get(paste0("Beta_Trait", trait)) + rnorm(.N, mean = 0, sd = sqrt(1 / (2 * sample_size * maf)))]
}

write.table(simulated.data,'simulated_beta_by_region.txt',
	    row.names=F,sep='\t',quote=F)
