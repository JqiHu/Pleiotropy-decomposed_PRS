# Check the patterns of subgroups
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/subgroup_V2_V5')
library(stringr)

subgroup <- list.files('.','subgroup')
subgroup <- subgroup[grep('RDS',subgroup)]
subgroup <- subgroup[-c(1,2)]
name <- str_replace_all(subgroup,'subgroup_','')
name <- str_replace_all(name,'.RDS','')
all <- readRDS('high_risk_ind.RDS')

# add all subgroup tables together
sum_sub <- c()
for(i in 1:length(subgroup)){
  temp <- readRDS(subgroup[i])
  temp2 <- cbind(name[i],temp$IID)
  colnames(temp2) <- c('Subgroup','IID')
  sum_sub <- rbind(sum_sub,temp2)
}
sum_sub <- as.data.frame(sum_sub)
dup <- table(sum_sub$IID)  

remain <- setdiff(all$IID,sum_sub$IID)
print(length(remain)) # 3173 not subgrouped

# Check the patterns of overlapping
pattern <- table(sum_sub$IID,sum_sub$Subgroup)
pattern2 <- as.data.frame.matrix(pattern)
## Create an indicator of patterns
index <- apply(pattern2,1,function(x){paste(colnames(pattern2)[x==1],collapse='_')})
print(length(table(index))) # 250 combinations
count <- table(index)
print(length(count[count>10])) # 114
print(length(count[count>100])) # 38

out <- cbind.data.frame(names(index),index)
colnames(out) <- c('eid','index')
write.csv(out,'pattern.csv',row.names=F)
## sum up the number of groups individauls in
num_subgroup <- apply(pattern2,1,sum)
print(table(num_subgroup))
