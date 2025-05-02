# aggregate GNOVA correlation results for 43 CAD-correlated traits

library(data.table)
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean//data/gnova_43/corr') # path to correlation files
a <- list.files('.') # all files

info <- fread('../../trait_43/trait_43_info.csv') # information for 43 traits
code <- info$Code
trait <- paste(info$Group,code,sep='_') # generate unique name for traits
code <- factor(code,level=code)
trait <- factor(trait,level=trait)


# split file name into trait1, trait1's code, trait2's code
library(stringr)
name <- str_replace(a,'.txt','')

group <- data.frame() # first column represents code for trait1 and second code for trait2
for(i in 1:length(a)){
  b <- str_split(name[i],'_')[[1]]
  group <- rbind(group,rev(b)[c(2,1)])
  group <- apply(group,2,as.character)
}

# function to input correlation file for a certain trait based on its code 
getCorr <- function(cod){ # enter code for trait
  cod <- as.character(cod)
  file <- list.files('.',pattern=paste0('_',cod)) # all the files including the trait
#  if(cod!=2092){file <- file[-grep(260,file)]}
  g <- group[which(a %in% file),] # get codes for traits involved in each file 
  tr <- c() # array of the code for second traits
  for(j in 1:nrow(g)){
    t <- ifelse(g[j,1]==cod,g[j,2],g[j,1])
    tr <- c(tr,t)
  }

  corr <- data.frame() # input correlation files
  for(i in 1:length(file)){
    corr <- rbind(corr,fread(file[i]))
  }
  corr$tr1 <- cod # code for trait1
  corr$tr2 <- tr # code for trait2
  corr$name <- str_replace(file,'.txt','') # file name
  return(corr)
}

# generate correlation file
corr_out <- c()
corr_p <- c()
for(i in 1:length(code)){
  cor <- getCorr(code[i]) # all correlation files for code[i]
  cor <- cor[cor$tr2 %in% code,] # remain 43 traits
  cor$tr2 <- factor(cor$tr2,level=code) 
  cor <- cor[order(cor$tr2),] 
  corr_out <- cbind(corr_out,cor$corr) # correlation coefficients matrix
  corr_p <- cbind(corr_p,cor$pvalue) # correlation p value matrix
}

# format matrix for plots
corr_out <- as.data.frame(corr_out)
corr_p <- as.data.frame(corr_p)
trait <- str_replace_all(trait,'_',' ')
rownames(corr_out) <- trait
colnames(corr_out) <- trait
rownames(corr_p) <- trait
colnames(corr_p) <- trait
