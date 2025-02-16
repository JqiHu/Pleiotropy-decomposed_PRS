# Combine results of clustering sensitivity analysis
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/clustering_sensitivity/hierachical_clustering/data')
library(stringr)
library(gridExtra)
library(ggplot2)

### read data
plot <- list()
clu <- c(5,7,9,15)
for(i in 1:4){
  add1 <- readRDS(paste0('./cluster_',clu[i],'/corr43_heatmap.RDS'))
  add2 <- readRDS(paste0('./cluster_',clu[i],'/prs_analysis/plot_rela_change.RDS'))
  
  gt <- add1$gtable
  
  # Adjust the search term based on your printed layout names
  row_idx <- which(gt$layout$name == "row_names")
  col_idx <- which(gt$layout$name == "col_names")

  # Replace with your new labels
  old.names <- gt$grobs[[row_idx]]$label
  new.names <- str_replace_all(old.names,'longetivity','longevity')
  new.names <- str_replace_all(new.names,'Vulvular','Valvular')
  gt$grobs[[row_idx]]$label <- new.names 
  gt$grobs[[col_idx]]$label <- new.names

  plot[[i*2-1]] <- gt
  plot[[i*2]] <- ggplotGrob(add2)
}

# Combine the two:
library(ggpubr)
# If heatmap_grob is from pheatmap and p_grob is from ggplotGrob:
pdf.options(reset = TRUE, onefile = FALSE)
pdf('All_res.pdf',width=30,height=40)
ggarrange(plot[[1]],plot[[2]],plot[[3]],plot[[4]],plot[[5]],plot[[6]],plot[[7]],plot[[8]],
	  ncol = 2, nrow=4,
	  widths = c(1,2),
	  labels = 'AUTO')

#grid.arrange(plot[[1]],plot[[2]],plot[[3]],plot[[4]],plot[[5]],plot[[6]],plot[[7]],plot[[8]],
#	     ncol=2, widths=c(1,1.8))

dev.off()
