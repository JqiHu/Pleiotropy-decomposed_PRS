# Use hierachical clustering to define clusters
# using 43 traits
setwd('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/review/PLOS_CB/clustering_sensitivity/hierachical_clustering/data')
library(factoextra) # for elbow figure
library(pheatmap)

# read correlations
dat <- readRDS('/gpfs/gibbs/pi/zhao/jh2875/ps_PRS/clean/data/gnova_43/data_for_supF1.RDS')
## correlation matrix
cor.coef <- dat[[1]]
cor_matrix <- as.matrix(cor.coef)
p_matrix <- as.matrix(dat[[2]])

#### Hierachical clustering
dist.matrix <- as.dist(1-abs(cor.coef))
hc <- hclust(dist.matrix)
# Compute the optimal number of clusters using WSS
pdf('hc_elbow.pdf')
fviz_nbclust(as.matrix(dist.matrix), FUN = hcut, method = "wss", k.max=30)

dev.off()

### Manually set the number of clusters 
#### Function to get cluster information
getCluster <- function(k,outprefix){
  clusters <- cutree(hc,k=k)
  clustered_features <- data.frame(Feature = colnames(cor.coef),
				   Cluster = clusters)
  clustered_features <- clustered_features[order(clustered_features$Cluster),]

  # Also generate heatmaps
  ## Threshold for significance (e.g., Bonferroni correction)
  threshold <- 0.05 / ((ncol(cor_matrix) * (ncol(cor_matrix) - 1)) / 2)
  ## Create a matrix with stars for significant correlations
  signif_matrix <- ifelse(p_matrix < threshold, "*", "")

  ### re-order correlation matrix
  order <- match(clustered_features$Feature,colnames(cor_matrix))
  cor_matrix <- cor_matrix[order,order]
  signif_matrix <- signif_matrix[order,order]

  # color of annotations
  ## Convert cluster labels to a factor
  cluster_labels <- as.factor(clusters[order])
  # Create annotation data frames
  annotation_df <- data.frame(Cluster = cluster_labels)
  # Ensure colors are correctly assigned
  cluster_levels <- levels(cluster_labels)
  annotation_colors <- list(Cluster = setNames(rainbow(length(cluster_levels)), cluster_levels))

  # Identify where cluster boundaries occur
  cluster_breaks <- which(duplicated(cluster_labels)[-1]==F)

  pdf(paste0(outprefix,'.pdf'),width=20,height=20)
  p <- pheatmap(cor_matrix,
           cluster_rows = F,
           cluster_cols = F,
           display_numbers = signif_matrix,
           border_color = "black",
           fontsize_row = 8,
           fontsize_col = 8,
           angle_col = 315,
           color = colorRampPalette(c("blue", "white", "red"))(100),
           breaks = seq(-1,1,length.out=100),
           annotation_row = data.frame(Cluster = as.factor(clusters[order])),
           annotation_col = data.frame(Cluster = as.factor(clusters[order])),
           annotation_colors = annotation_colors,
           main = paste0("Genetic Correlation Heatmap"),
           gaps_row = cluster_breaks,  # Add horizontal lines between clusters
           gaps_col = cluster_breaks)  # Add vertical lines between clusters
  print(p)

  dev.off()

  saveRDS(p,paste0(outprefix,'.RDS'))

  return(clustered_features)
}


###### Try 1. k=9 [same as main]
k=9
out <- getCluster(k,'./cluster_9/corr43_heatmap')
write.csv(out,'./cluster_9/cluster_info.csv',row.names=F)

###### Try 2. k=7 [seem to be the elbow point]
k=7
out <- getCluster(k,'./cluster_7/corr43_heatmap')
write.csv(out,'./cluster_7/cluster_info.csv',row.names=F)

###### Try 3. k=15 [suggested by reviewer]
k=15
out <- getCluster(k,'./cluster_15/corr43_heatmap')
write.csv(out,'./cluster_15/cluster_info.csv',row.names=F)

###### Try 4. k=5 [suggested by reviewer]
k=5
out <- getCluster(k,'./cluster_5/corr43_heatmap')
write.csv(out,'./cluster_5/cluster_info.csv',row.names=F)






