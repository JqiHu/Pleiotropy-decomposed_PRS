# generate plot of subgroups
rm(list=ls())
library('data.table')
library('stringr')
library('ggplot2')
setwd('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/subgroup/')

## data contains IID, subgroup, psprs proportions in long format
long <- fread('psprs_prop_long.tsv')
long <- as.data.frame(long)
path <- unique(sort(long$psprs)) # psprs and subgroup name in alphabetical order

## data contains CAD PRS
cad_prs <- fread('../../update_overlap/Explorations/subgroup_version2/cad_prs.tsv')
iid <- cad_prs$IID[order(cad_prs$SCORESUM,decreasing=T)] # rank by CAD PRS from large to small

long$IID <- factor(long$IID,levels=iid)
long$psprs <- factor(long$psprs, levels=path)

#### Calculate mean contribution of psPRSs comparing subgroup vs. remaining
## value is the absolute value of CAD PRS explained by PD-PRS
###### step 1. transform absolute values to ratio
long$index <- paste0(long$IID,'_',long$subgroup)
sum_prs <- aggregate(long$value,list(long$index),sum)
long$sum_prs <- sum_prs$x[match(long$index,sum_prs$Group.1)]
long$ratio <- long$value/long$sum_prs
mean_contri <- c()
for(i in 1:9){
  subgroup <- long[long$subgroup==path[i],]
  remain <- long[-which(long$subgroup==path[i]),]
  
  ## absolute values
  add1 <- aggregate(subgroup$ratio,list(subgroup$psprs),mean)
  add2 <- aggregate(remain$ratio,list(remain$psprs),mean)
  ## SD
  add_sd1 <- aggregate(subgroup$ratio,list(subgroup$psprs),sd)
  add_sd2 <- aggregate(remain$ratio,list(remain$psprs),sd)
  
  add <- round(rbind(add1$x,add2$x)*100,1)
  add_sd <- round(rbind(add_sd1$x,add_sd2$x),2)
  
  add_out <- rbind(paste0(add,'% (',add_sd,')')[seq(1,18,2)],
                   paste0(add,'% (',add_sd,')')[seq(2,18,2)])
  
  rownames(add_out) <- c(path[i],paste0(path[i],'_Remain'))
  mean_contri <- rbind(mean_contri,add_out)
}
colnames(mean_contri) <- path
mean_contri <- cbind.data.frame(rownames(mean_contri),mean_contri)
colnames(mean_contri)[1] <- 'Subgroup'
library(writexl)
write_xlsx(mean_contri,'mean_psPRS_contribution_to_subgroup.xlsx')

set.seed(116)
library('Cairo')
library(ggpubr)

## color for each psPRS
col <- fread('../visualization//path.col')
color <- col$col[c(1:5,7,6,8:9)]
names(color) <- path

p <- list() # list for plots
for(i in 1:length(path)){
  temp <- long[which(long$subgroup==path[i]),] # extract individuals in this subgroup
  temp <- temp[order(temp$IID),] # rank by CAD PRS
  temp$value <- as.numeric(temp$value)
  n <- length(unique(temp$IID)) # number of individuals
  print(n)

  p[[i]] <- ggplot(temp,aes(x=IID,y=value,fill=psprs))+
         geom_bar(stat="identity",width = 0.6, position = 'stack',alpha=1)+
         ylab(paste0(path[i],' subgroup'))+
         xlab(paste0('N = ',n))+ # n is the number of individuals* 9 psprs
         scale_fill_manual(values=color)+
         scale_y_continuous(expand = c(0, 0),limits=c(0,0.75))+
         guides(fill=guide_legend(title="PD-PRS"))+
         theme_classic()+
         theme(axis.title = element_text(size=12,face = 'bold'),
               axis.text.x = element_blank(),
               axis.text.y = element_text(size = 8,face = 'plain'),
               axis.line.y = element_line(colour=NULL))
}

# output plot and arrange them as 3*3
CairoPNG('../visualization//subgroup_stack.png',
	bg='transparent',width=600*3,height=600*3,res=72*3)

ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
          ncol=3,nrow=3,common.legend=T,legend='right')
dev.off()

