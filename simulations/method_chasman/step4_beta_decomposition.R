# SVD on GWAS matrix and get residuals
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/beta_component')
library(data.table)

# Read GWAS
dat <- fread('../sumstats/scaled_simu_causalfortrait1.txt',data.table=F)

# Matrix B_pheno
# with phenotypes in rows and 
# SNPs in columns
b_pheno <- as.matrix(t(dat[,28:30]))

# Check means of Beta across phenotypes
print(apply(b_pheno,2,mean)) # clearly not close to 0....

#### SVD on Bpheno
svd_res <- svd(b_pheno)
# Components of SVD
U <- svd_res$u           # Left singular vectors
# Extend V 
V_partial <- svd_res$v       # Right singular vectors 
V_full <- qr.Q(qr(cbind(V_partial, matrix(0, nrow = nrow(V_partial), ncol = nrow(V_partial)-ncol(V_partial)))))  # Pad to full dimension

# Expanded D to match 
D_expanded <- matrix(0, nrow = nrow(b_pheno), ncol = ncol(b_pheno))
diag(D_expanded) <- svd_res$d

# Reconstruct the original matrix
b_phe_re <- U %*% D_expanded %*% t(V_full) # arbitrary directions

prop <- svd_res$d^2/sum(svd_res$d^2)
n_comp <- 3 # cumulative variance > 80%

V_r <- V_full[,1:n_comp]

#### switch direction!!!!!!!!!!
V_r <- (-1)*V_r

###### Get residuals 
B_disease <- as.matrix(dat[,27])
tmp <- cbind.data.frame(B_disease,V_r)
colnames(tmp)[2:4] <- paste0('comp',1:3)

fit <- lm(B_disease~comp1+comp2+comp3,data=tmp)
B_resid <- fit$residuals

#### Out put weights for PRS
out <- cbind.data.frame(dat[,1:6],tmp,B_resid)

write.table(out,'weights_for_PRS.txt',
	    row.names=F,sep='\t',quote=F)
