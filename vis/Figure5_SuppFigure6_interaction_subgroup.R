# generate HR forest + ARR plots
rm(list=ls())
library(data.table)
library(ggplot2)
library(stringr)
library(Hmisc)
library(forcats)
#library(forestploter)
library(patchwork)

## input arr data
setwd('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/interactions/')

######################### 
# HR Forest
data <- readRDS('inter_tranas_select_update_vis.RDS')
name <- names(data)
spl <- unlist(strsplit(name,'*',fixed = T))
traits <- spl[seq(1,length(spl),2)]
trait <- c('Cholesterol','LDL direct','Triglyceride','Sleep duration',
           'FVC','FEV1','PEF','DBP','SBP',
           'Current smoking (Most or all days)','Current smoking (Most or all days)',
           'Current smoking','Current smoking')
prs <- spl[seq(2,length(spl),2)]
prs <- str_replace_all(prs,"Respiratory_syste",'Respir')

# Forest plot one by one
hr_forest <- list()
for(i in 1:length(name)){
  temp <- as.data.frame(data[[i]])
  temp$n.case <- as.numeric(temp[,2])[5:8]
  temp$n.ctrl <- as.numeric(temp[,2])[1:4]
  temp$Group <- c(rep(paste0(traits[i],' = 0'),2),
                  rep(paste0(traits[i],' = 1'),2))  
  temp <- temp[,c(1,3:9)]
  temp$HR <- round(temp$`exp(coef)`,2)
  temp$lower <- round(exp(temp$coef-1.96*temp$`se(coef)`),2)
  temp$upper <- round(exp(temp$coef+1.96*temp$`se(coef)`),2) 
  temp$se <- (temp$upper-temp$HR)/1.96
  
  temp$`HR (95% CI)` <- paste0(temp$HR,' (',temp$lower,'-',temp$upper,')')
  temp$`HR (95% CI)`[1] <- '1 (Ref)'
  
  # change the sequence of groups
  # HR00, HR-phe, HR-prs, HR11
  # Add group description to the 'Group' column
  temp2 <- temp[c(1,3,2,4),]
  temp2$Group1 <- c('Reference','Remaining High Risk_Exposed',
                   paste0(prs[i],' Subgroup_Unexposed'),
                   paste0(prs[i],' Subgroup_Exposed'))
  temp2$Group <- factor(temp2$Group1,levels=temp2$Group1)
  
  ###  Point estimate
  p1 <- ggplot(data=temp2,aes(y=fct_rev(Group)))+
    theme_bw() +
    geom_point(aes(x=HR), shape=15, size=3) +
    geom_linerange(aes(xmin=lower, xmax=upper)) +
    geom_vline(xintercept = 1, linetype="dashed") +
    scale_y_discrete(labels=rev(c('Reference','Remain Exposed',
                              paste0(prs[i],' Subgroup\n Unexposed'),
                              paste0(prs[i],' Subgroup\n Exposed'))))+
    theme(plot.title = element_text(hjust = 0,size=7,face='bold'),
          axis.title.y= element_blank(),
          text = element_text(size=7,face='bold')) +
    ggtitle(paste0(trait[i],'*',prs[i],' PD-PRS'))+
    coord_cartesian(xlim=c(min(temp2$lower), max(temp2$upper)))  # zoom in
  # ### Place Text
  # text <- ggplot(data=temp2)+
  #   coord_flip()+
  #   geom_text(aes(x=fct_rev(Group),y=0,label=`HR (95% CI)`),
  #             hjust=2,size=4.5,fontface='bold')+
  #   ggtitle('HR (95% CI)')+
  #   theme_bw()+
  #   theme(text = element_text(size=10,face='bold'),
  #         axis.text = element_blank(),
  #         axis.ticks = element_blank(),
  #         panel.border = element_blank(),
  #         panel.grid.major=element_blank(),
  #         panel.grid.minor=element_blank(),
  #         plot.margin = margin(5.5,1,5.5,5.5))+
  #   labs(x=NULL,y=NULL)
  
  hr_forest[[paste0(traits[i],':',prs[i])]] <- p1
  
}

######################### 
# ARR
arr <- fread('sub_arr_trans_selected_update.tsv')
arr <- as.data.frame(arr)
arr <- arr[is.na(arr$arr)==F,]
rownames(arr) <- 1:nrow(arr)

## subset data on specific phenotypes and psPRSs
temp <- arr
temp$subgroup <- factor(c('High','High','Remaining','Remaining'),
                        levels=c('High','Remaining'))
temp$prs <- arr$subgroup

## functions generating plot based on subgroup and traits
library(ggthemes)
library(ggsci)

getPlot <- function(data,name,path){
  data <- data[data$prs %in% c(path,paste0(path,'_remain')),]
  data$phenotype <- factor(c('High','Low'),
                           levels=c('High','Low'))
  
  reri <- (data$arr[1]-data$arr[3])/data$ar[4]
  p <- ggplot(data=data,aes(x=subgroup,y=ar,fill=phenotype))+
    geom_bar(stat ="identity",position = "dodge")+
    scale_fill_jama(name=name, labels=c('Exposed','Nonexposed'),)+
    guides(fill=guide_legend(title=""))+
    xlab(paste0(path,'-PD-PRS'))+
#    ylab('Incident rate of CAD')+
    theme_bw()+
    theme(legend.position='right',axis.title = element_text(size=5,face = 'plain'),
          title = element_text(size=10,face='plain'),
          legend.title = element_text(size=5, face='plain'),
          legend.background = element_rect(size=0.2, linetype='solid',colour='black'),
          axis.title.y = element_blank(),
          axis.text.x = element_text(size = 10,hjust = 0.2,face = 'plain',vjust=0),
          axis.text.y = element_text(size = 10,face = 'plain',angle=0, hjust=0.2, vjust=0))
  tinyincrease=0.01
  angle=0
  arr_label <- paste0('ARR = ', sprintf('%.1f',unique(data$arr)*100),'%')
  p_label = p +
    annotate('text',
             x=c(0.83, 1.83),
             y=c(data$ar[c(1,3)]+tinyincrease*3),
             label=arr_label,
             size= 10/3.2, angle = angle, hjust = -0.005)+
    annotate("segment", # vertical
             x=c(0.77, 1.78),
             xend=c(0.77, 1.78),
             y= data$ar[c(1,3)],
             yend=data$ar[c(1,3)] + tinyincrease * 2) +
    annotate("segment", # vertical 
             x=c(1.22, 2.22),
             xend=c(1.22, 2.22),
             y= data$ar[c(2,4)],
             yend=data$ar[c(1,3)] + tinyincrease * 2) +
    annotate("segment", # horizontal
             x=c(0.77, 1.78),
             xend=c(1.22, 2.22),
             y= data$ar[c(1,3)] + tinyincrease * 2,
             yend=data$ar[c(1,3)] + tinyincrease * 2) +
    annotate('text',
             x=1.33,
             y=max(data$ar[c(1,3)])+tinyincrease*6,
             label=paste0('RERI = ',sprintf('%.1f',reri*100),'%'),
             size= 10/3.2, angle= angle, hjust=-0.005)+
    annotate('segment', # vertical
             x=c(1, 2),
             xend=c(1, 2),
             y=c(data$ar[c(1,3)]+tinyincrease*3.5),
             yend=rep(max(data$ar[c(1,3)])+tinyincrease*5,2))+
    annotate("segment", # horizontal
             x=1,
             xend=2,
             y=max(data$ar[c(1,3)]) + tinyincrease * 5,
             yend=max(data$ar[c(1,3)]) + tinyincrease * 5) 
  return(p_label)
}

# plot figures
arr_p <- list()
trait <- temp$trait[seq(1,nrow(temp),4)]
prs <- arr$subgroup[seq(1,nrow(temp),4)]

for(i in 1:length(trait)){
  arr_p[[i]] <- getPlot(temp[temp$trait==trait[i],],trait[i],prs[i])

}

######################### 
# Combine into TWO figures
library(ggpubr)
library(gridExtra)
### Main figures: LDL, DBP, SBP, smoking current, frequent smoking current
f1 <- ggarrange(hr_forest[[8]],hr_forest[[9]],arr_p[[8]],arr_p[[9]],
                hr_forest[[10]],hr_forest[[12]],arr_p[[10]],arr_p[[12]],
                align = 'hv',
                heights = c(0.5,1,0.5,1,0.5,1,0.5,1),
                labels = c('A','B','C','D','E','F','G','H'),
                ncol=2,nrow=4,common.legend = T) + 
  theme(plot.margin = margin(0.1,2,0.2,0.1, "cm"))

print(f1)


pdf('../visualization/main_figure_update.pdf',height=10,width = 12)
print(f1)

dev.off()

## Supp
f1 <- ggarrange(hr_forest[[1]],hr_forest[[2]],hr_forest[[3]],arr_p[[1]],arr_p[[2]],arr_p[[3]],
                hr_forest[[5]],hr_forest[[6]],hr_forest[[7]],arr_p[[5]],arr_p[[6]],arr_p[[7]],
                hr_forest[[11]],hr_forest[[13]],hr_forest[[4]],arr_p[[11]],arr_p[[13]],arr_p[[4]],
                align = 'hv',
                heights = c(rep(c(0.5,1),3),rep(c(0.5,1),3),0.5,0.5,0.5,1,1),
                labels = c('A','B','C','D','E','F',
                           'G','H','I','J','K','L',
                           'M','N','O','P','Q','R'),
                ncol=3,nrow=6,common.legend = T) + 
  theme(plot.margin = margin(0.1,2,0.2,0.1, "cm"))
print(f1)

pdf('../visualization/sup_figure_update.pdf',height=14,width=12)
print(f1)

dev.off()



