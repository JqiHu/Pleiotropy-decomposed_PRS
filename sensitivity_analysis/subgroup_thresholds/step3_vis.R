# Supplemental Figure 5
library('data.table')
library('stringr')
rm(list=ls())

### Generate Figure objects
plot <- list()
for(tresholds in c('0.01','0.1')){
  setwd(paste0('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/subgroup_sensitivity/data/',thresholds))
  ## data contains IID, subgroup, psprs proportions in long format
  long <- fread('psprs_prop_long.tsv')
  long <- as.data.frame(long)
  path <- unique(sort(long$psprs)) # psprs and subgroup name in alphabetical order

  ## data contains CAD PRS
  cad_prs <- fread('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_overlap/Explorations/subgroup_version2/cad_prs.tsv')
  iid <- cad_prs$IID[order(cad_prs$SCORESUM,decreasing=T)] # rank by CAD PRS from large to small

  long$IID <- factor(long$IID,levels=iid)
  long$psprs <- factor(long$psprs, levels=path)

  #### Calculate mean contribution of psPRSs comparing subgroup vs. remaining
  ###### step 1. transform absolute values to ratio
  long$index <- paste0(long$IID,'_',long$subgroup)
  sum_prs <- aggregate(long$value,list(long$index),sum)
  long$sum_prs <- sum_prs$x[match(long$index,sum_prs$Group.1)]
  long$ratio <- long$value/long$sum_prs
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
  }

  set.seed(116)
  library(ggpubr)

  ## color for each psPRS
  col <- fread('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/visualization//path.col')
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
  p.out <- ggarrange(p[[1]],p[[2]],p[[3]],p[[4]],p[[5]],p[[6]],p[[7]],p[[8]],p[[9]],
                     ncol=3,nrow=3,common.legend=T,legend='right')
  plot[[paste0('stack_',thresholds)]] <- p.out



  # Relative change
  getPlot <- function(data,group){ # input data for specific phenotype group and group name
    p <-   ggplot(data,aes(x=`pheno`,,y=`change`,fill=`psPRS`))+
      geom_bar(stat='identity',position='dodge',width = 0.8,alpha=0.7)+
      scale_fill_manual(values=color)+
      #ggtitle(paste0(group,'-related traits'))+
      theme_bw()+
      xlab('')+
      ylab(paste0('Relative change in ',group,'-related traits'))+
      guides(fill=guide_legend(title="PD-Subgroup"))+
      theme(axis.title = element_text(size=11),
            axis.text.x = element_text(size=11, face = 'plain'),
            axis.text.y = element_text(size = 11,face = 'plain'))+
      coord_flip()
    return(p)
  }

  p.list <- list()
  for(group in c('BP','Lipids','Obesity-obe',
                 'Respiratory','T2D')){
    rela <- read_xlsx('relative_change_all.xlsx',sheet = group)

    ## formate name for triats
    rela$pheno <- str_replace_all(rela$pheno,'_',' ')
    rela$pheno[rela$pheno=='hypertation'] <- 'hypertension'
    rela$pheno[rela$pheno=='bmi'] <- 'BMI'
    rela$pheno <- capitalize(rela$pheno)
    path <- sort(unique(rela$psPRS)) # names for psPRS in alphabetical order
    names(color) <- path

    rela <- as.data.frame(rela)
    rela$psPRS <- factor(rela$psPRS,levels=path)
    rela <- rela[is.na(rela$psPRS)==F,]

    add <- getPlot(rela,group)
    p.list[[group]] <- add
  }

  # output plots
  library(ggpubr)
  p.out <- ggarrange(p.list[['BP']]+ggtitle('(A) BP-related traits'),
                     p.list[['Lipids']] + ggtitle('(B) Lipids-related traits'),
                     p.list[['Obesity-obe']] + ggtitle('(C) Obesity-related traits'),
                     p.list[['Respiratory']] + ggtitle('(D) Smoking-related traits'),
                     p.list[['T2D']] + ggtitle('(E) T2D-related traits'),
                     ncol=3,nrow=2,common.legend=T,
                     legend='right',align='hv')
  plot[[paste0('rela_change_',tresholds)]] <- p.out
}

library(Cairo)

combined_plot <- ggarrange(plot.list=plot,ncol=2,nrow=2,
                           labels='AUTO')


CairoPNG('subgroup_sensitivity.png',
         bg='white',width=600*20,height=600*15,res=72*5)
print(combined_plot)

dev.off()
