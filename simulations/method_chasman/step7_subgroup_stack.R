# Visualize the composition of PD-PRSs for subgroups
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/simulation/chasman/data/analysis')
library(data.table)
library(ggplot2)

# read data
#dat <- fread('simulated_subgroup.txt',data.table=F)
dat <- fread('inter_simulated_subgroup.txt',data.table=F)
ps_names <- c('PRS_comp1','PRS_comp2','PRS_comp3','PRS_B_resid')
subgroup_names <- paste0('subgroup',1:4)
# High-risk subjects
dat.high <- dat[dat$high_risk==1,]

subgroup <- c()
for(i in 1:length(ps_names)){
  # Individuals in the subgroup
  out <- dat.high[dat.high[,subgroup_names[i]]==1,]
  print(nrow(out))
  if(nrow(out)==0){
    next
  }
  out <- cbind.data.frame(ps_names[i],out)
  subgroup <- rbind(subgroup,out)
}
colnames(subgroup)[1] <- 'subgroup'
#### add remaining subjects
add_remain <- dat.high[which(dat.high$IID %in% subgroup$IID==F),]
if(nrow(add_remain)>0){
  add_remain$subgroup <- 'remain'
  subgroup <- rbind(subgroup,add_remain)
}

## calculate proportions of each PD-PRS on overall
sum_psprs <- apply(abs(dat.high[,3:6]),1,sum) # sum of 4 PD-PRSs
# absolute value of overall PRS explained by PD-PRS
prop <- apply(abs(dat.high[,3:6]),2,
                function(x){(abs(x)/sum_psprs)*dat.high$PRS_B_disease})

# Combine with subgroup information
wide <- cbind.data.frame(subgroup$subgroup,
                         dat.high$IID[match(subgroup$IID,dat.high$IID)],
                         prop[match(subgroup$IID,dat.high$IID),])
colnames(wide)[1:2] <- c('subgroup','IID')

long <- melt(setDT(wide),id.vars=c('IID','subgroup'),
                variable.name='psprs') # long format for plot
## rank by CAD PRS
iid <- dat.high$IID[order(dat.high$PRS_B_disease,decreasing=T)]
long$IID <- factor(long$IID,levels=unique(iid))

# ------------------------------------
#### Calculate mean contribution of psPRSs comparing subgroup vs. remaining
###### step 1. transform absolute values to ratio
long$index <- paste0(long$IID,'_',long$subgroup)
sum_prs <- aggregate(long$value,list(long$index),sum)
long$sum_prs <- sum_prs$x[match(long$index,sum_prs$Group.1)]
long$ratio <- long$value/long$sum_prs
mean_contri <- c()
for(i in 1:4){
  subgroup <- long[long$subgroup==ps_names[i],]
  if(nrow(subgroup)==0){next}
  remain <- long[-which(long$subgroup==ps_names[i]),]

  ## absolute values
  add1 <- aggregate(subgroup$ratio,list(subgroup$psprs),mean)
  add2 <- aggregate(remain$ratio,list(remain$psprs),mean)
  ## SD
  add_sd1 <- aggregate(subgroup$ratio,list(subgroup$psprs),sd)
  add_sd2 <- aggregate(remain$ratio,list(remain$psprs),sd)

  add <- round(rbind(add1$x,add2$x)*100,1)
  add_sd <- round(rbind(add_sd1$x,add_sd2$x),2)

  add_out <- rbind(paste0(add,'% (',add_sd,')')[seq(1,8,2)],
                   paste0(add,'% (',add_sd,')')[seq(2,8,2)])

  rownames(add_out) <- c(paste0('subgroup',i),paste0(paste0('subgroup',i),'_Remain'))
  mean_contri <- rbind(mean_contri,add_out)
}
mean_contri <- cbind.data.frame(rownames(mean_contri),mean_contri)
colnames(mean_contri)[1] <- 'Subgroup'
library(writexl)
write_xlsx(mean_contri,'inter_mean_psPRS_contribution_to_subgroup.xlsx')

set.seed(116)
library('Cairo')
library(ggpubr)

## color for each psPRS
color <- c('tomato','skyblue','darkgreen','gold')
names(color) <- ps_names
### label for each psPRS
label <- c(paste0('Component ',1:3,' PRS'),
	   'Residual PRS')
names(label) <- ps_names

p <- list() # list for plots
for(i in 1:4){
  temp <- long[which(long$subgroup==ps_names[i]),] # extract individuals in this subgroup
  temp <- temp[order(temp$IID),] # rank by CAD PRS
  temp$value <- abs(as.numeric(temp$value))
  n <- length(unique(temp$IID)) # number of individuals
  print(n)

  temp$psprs <- factor(temp$psprs,levels=c(paste0('PRS_comp',1:3),'PRS_B_resid'))

  ylab.name <- ifelse(i==4,'Residual PRS Subgroup',
		      paste0('Component ',i,'-related Subgroup'))

  p[[i]] <- ggplot(temp,aes(x=IID,y=value,fill=psprs))+
         geom_bar(stat="identity",width = 0.6, position = 'stack',alpha=1)+
         ylab(ylab.name)+
         xlab(paste0('N = ',n))+ # n is the number of individuals
         scale_fill_manual(values=color,labels=label)+
         scale_y_continuous(expand = c(0, 0),limits=c(0,5.2))+
         guides(fill=guide_legend(title="PD-PRS"))+
         theme_classic()+
         theme(axis.title = element_text(size=12,face = 'bold'),
               axis.text.x = element_blank(),
               axis.text.y = element_text(size = 8,face = 'plain'),
               axis.line.y = element_line(colour=NULL),
               legend.text = element_text(size = 14),
               legend.title = element_text(size = 16),
	       plot.margin = unit(c(0, 0, 0, 0), "cm"))
}

# output plot and arrange them as 2*2
CairoPNG('inter_subgroup_stack.png',
	 bg='white',width=600*3,height=600*3,res=72*3)

p.out <- ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],
               ncol=2,nrow=2,common.legend=T,legend='right')
print(p.out)

dev.off()

saveRDS(p,'inter_subgroup_stack.RDS')
