# generate plot for relative change comparing 
# subgroup and the rest high-risk CAD
rm(list=ls())
setwd('~/Desktop/Research/ZHAO_HONGYU/pathwayPRS/update_V2_absolute/subgroup/')
library(data.table)
library(ggplot2)
library(Hmisc)
library(stringr)
library(readxl)

## color for each psPRS
col <- fread('../visualization//path.col')
color <- col$col[c(1:5,7,6,8:9)]

# plots
library('ggplot2')
getPlot <- function(data,group){ # input data for specific phenotype group and group name
  p <-   ggplot(data,aes(x=`pheno`,,y=`change`,fill=`psPRS`))+
    geom_bar(stat='identity',position='dodge',width = 0.8,alpha=0.7)+
    scale_fill_manual(values=color)+
    #ggtitle(paste0(group,'-related traits'))+
    theme_bw()+
    xlab('')+
    ylab(paste0('Relative change in ',group,'-related traits'))+
    guides(fill=guide_legend(title="PD-Subgroup"))+
    theme(axis.title = element_text(size=12, face='bold'),
          axis.text.x = element_text(size=12, face = 'bold'),
          axis.text.y = element_text(size = 12,face = 'bold'))+
    coord_flip()
  return(p)
}

plot <- list()
for(group in c('BP','Lipids','Obesity-obe',
               'Respiratory','T2D')){
  rela <- read_xlsx('relative_change_all.xlsx',sheet = group)
  if(group=='Obesity-obe'){group='Obesity'}
  
  ## formate name for triats
  rela$pheno <- str_replace_all(rela$pheno,'_',' ')
  rela$pheno[rela$pheno=='hypertation'] <- 'hypertension'
  rela$pheno[rela$pheno=='bmi'] <- 'BMI'
  rela$pheno <- capitalize(rela$pheno)
  path <- sort(unique(rela$psPRS)) # names for psPRS in alphabetical order
  names(color) <- path 
#  path <- setdiff(path,'Others') # remove others
  
  rela <- as.data.frame(rela)
  rela$psPRS <- factor(rela$psPRS,levels=path)
  rela <- rela[is.na(rela$psPRS)==F,]
  
  add <- getPlot(rela,group)
  plot[[group]] <- add
}

# output plots
library(ggpubr)
pdf('../visualization//rela_change.pdf',
	width=18,height=9)
ggarrange(plot[['BP']]+ggtitle('(A) BP-related traits'),
  plot[['Lipids']] + ggtitle('(B) Lipids-related traits'),
  plot[['Obesity']] + ggtitle('(C) Obesity-related traits'),
  plot[['Respiratory']] + ggtitle('(D) Smoking-related traits'),
  plot[['T2D']] + ggtitle('(E) T2D-related traits'),
  ncol=3,nrow=2,common.legend=T,
  font.label = list(size = 14, color = "black", face = "bold", family = NULL),
  legend='right',align='hv')

dev.off()


pdf('../visualization/rela_change_Lipids.pdf',
    width=8,height = 6)
print(plot[['Lipids']]+ggtitle('Relative changes in Lipids-related traits'))

dev.off()


