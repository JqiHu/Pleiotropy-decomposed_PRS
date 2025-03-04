# Generate table and figures for 
# subgroup patterns
library(data.table)
library(ggplot2)
library(stringr)

data <- fread('../subgroup/pattern.csv')
data$index <- str_replace_all(data$index,'Basic_Condition','Basic condition')
data$index <- str_replace_all(data$index,'Immune_System','Immune system')
data$index <- str_replace_all(data$index,'Respiratory_system','Respiratory system')
data$index <- str_replace_all(data$index,'_',' & ')

# Table
table <- table(data$index)
out <- cbind.data.frame(names(table),table)
colnames(out)[1] <- 'Subgroup Pattern'
library(writexl)
write_xlsx(out,'../subgroup/subgroup_pattern.xlsx')

# Figure
for_figure <- as.data.frame(out[out$Freq>=100,]) # N=38
for_figure <- for_figure[order(for_figure$Freq,decreasing = T),]
for_figure$Freq <- as.numeric(for_figure$Freq)
for_figure$`Subgroup Pattern` <- factor(for_figure$`Subgroup Pattern`,
                                        levels=for_figure$`Subgroup Pattern`)
for_figure$num_group <- '1'
for_figure$num_group[grep('&',for_figure$`Subgroup Pattern`)] <- '2'

## Histogram
library(ggsci)
p <- ggplot(for_figure,aes(x=`Subgroup Pattern`,y=Freq,fill=num_group))+
  geom_bar(stat="identity")+
  scale_fill_jama()+
  ylab('Number of Subjects')+
  theme_bw()+
  theme(legend.position="none",
        axis.title = element_text(size=12,face = 'bold'),
        axis.text.x = element_text(size = 12,angle = -45,hjust = 0,face = 'bold'),
        axis.text.y = element_text(size = 10,face = 'plain'),
        axis.line.y = element_line(colour=NULL),
        plot.margin = margin(10, 90, 10, 10))
print(p)

pdf('subgroup_pattern_le_100.pdf',width = 12,height = 7)
print(p)
dev.off()
