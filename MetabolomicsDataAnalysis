library(SNFtool)
library(ggplot2)
library(tidyverse)
library(readxl)
library(openxlsx)
library(WGCNA)
library(readxl)
library(reshape2)
library(janitor)
library(patchwork)
library(MetaboAnalystR)
library(preprocessCore)



setwd("/Users/yangzhang/Downloads/MP_Bio_Package/")
getwd()

c18_data <- read.xlsx("C18_ChemicalFormula.xlsx", sheet = "Compounds", colNames = TRUE)
samples_info <- as.data.frame(read_excel("Samples.xlsx", col_names =TRUE))
c18_data <- as.data.frame(c18_data) %>% 
  clean_names()
c18_x <- rownames(c18_data)

c18_x <- as.data.frame(c18_x)
c18_data_filtered <- cbind(c18_x, c18_data)

c18_data_filtered <- c18_data_filtered %>%
  filter(!is.na(checked))


test_c18_data_filtered1 <- c18_data_filtered[!is.na(c18_data_filtered[,'annot_delta_mass_ppm']),]
test_c18_data_filtered2 <- test_c18_data_filtered1[!is.na(test_c18_data_filtered1[,'rt_min']),]


# for (i in 25:35) {
#   print(paste("Hi",i, collapse =""))
# }
# colnames(test_c18_ratio_sb)
# str_sub(colnames(test_c18_ratio_sb)[1],1,4)



hist(abs(as.numeric(test_c18_data_filtered2[,"annot_delta_mass_ppm"])), xlab="Intensity", col="lightgreen", main="annot_delta_mass_ppm_c18",breaks = 40)
x <- quantile(abs(as.numeric(test_c18_data_filtered2[,"annot_delta_mass_ppm"])), probs = seq(0, 1, 0.05))

xidx <- which(abs(as.numeric(test_c18_data_filtered2[,"annot_delta_mass_ppm"])) >= x[20], arr.ind = T)
view(test_c18_data_filtered2[xidx,])
test_c18_data_filtered3 <- test_c18_data_filtered2[-xidx,]
hist(as.numeric(test_c18_data_filtered3[,"annot_delta_mass_ppm"]), xlab="Intensity", col="lightgreen", main="ppm_outliersRemoved_c18",breaks = 40)


# select ratio_omics_system blank (14 columns)
c18_data_flitered_sb <- mapply(test_c18_data_filtered3[,45:58],FUN=as.numeric)


#check nan values in the target data
nanidx <- which(c18_data_flitered_sb==NA, arr.ind = T)
length(nanidx)

# print(colnames(c18_data_flitered_sb))
# boxplot(t(c18_data_flitered_sb))
#c18_data_flitered_sb1 <- rbind(seq(1, dim(c18_data_flitered_sb)[1], 1), t(c18_data_flitered_sb))

medata <- melt(c18_data_flitered_sb) %>% 
  filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

#hist(Expresdata2, xlab="Intensity", col="lightgreen", main="ExpressionValue before log & normalization")

NLExp.plot <- stack(as.data.frame(c18_data_flitered_sb))
# ggplot(NLExp.plot, aes(x=values)) +
#   geom_histogram(binwidth=20, colour="black", fill="white")
# 
# 
# den <- density(NLExp.plot$values)
# p1 <- plot(den, frame = FALSE, col = "blue",main = "Density plot")

p1 <- NLExp.plot %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'c18_original_ratio_SysBlank', theme = theme(plot.title = element_text(hjust = 0.5)))

########### Three steps normalization for metabolomics data: 1. Sample normalization -> 2. log transformation -> 3. Feature normalization
# Sample normalization: sum normalization
c18_colSUM <- colSums(c18_data_flitered_sb, dims = 1)
c18_sample_norm <- data.frame()
numeratorM <- matrix(rep(c18_colSUM, dim(c18_data_flitered_sb)[1]), nrow=dim(c18_data_flitered_sb)[1], byrow=T)
c18_sample_norm <- c18_data_flitered_sb/numeratorM
# Sample normalization: quantile normalization
c18_quanorm_sb <- normalize.quantiles(c18_data_flitered_sb, copy = F)

# log transfermation
#log_c18_sample_norm <- log2(c18_sample_norm)
log_c18_sample_norm <- log(c18_quanorm_sb)
den <- density(stack(as.data.frame(log_c18_sample_norm))$values)
p1 <- plot(den, frame = FALSE, col = "blue",main = "Density plot")

# feature normalization
c18_rowmean <- apply(log_c18_sample_norm, 1, mean)
rowmeanM <- t(matrix(rep(c18_rowmean,dim(log_c18_sample_norm)[2]), nrow = dim(log_c18_sample_norm)[2], byrow=T))
c18_sd <- apply(log_c18_sample_norm, 1, sd)
sdM <- t(matrix(rep(c18_sd,dim(log_c18_sample_norm)[2]), nrow = dim(log_c18_sample_norm)[2], byrow=T))
c18_normalized_sb <- (log_c18_sample_norm-rowmeanM)/sdM

NLExp.plotsb <- stack(as.data.frame(c18_normalized_sb))
p1 <- NLExp.plotsb %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

medata_sb <- melt(c18_normalized_sb) %>% 
  dplyr::filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata_sb, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'c18_normalized_ratio_SysBlank', theme = theme(plot.title = element_text(hjust = 0.5)))




#######################################################     ratio: system blank resuspend
# select ratio_omics_system blank resuspend (12 columns)
c18_data_flitered_sbr <- mapply(test_c18_data_filtered3[,59:70],FUN=as.numeric)

#check nan values in the target data
nanidx <- which(c18_data_flitered_sbr==NA, arr.ind = T)
length(nanidx)

# print(colnames(c18_data_flitered_sb))
# boxplot(t(c18_data_flitered_sb))

medata <- melt(c18_data_flitered_sbr) %>% 
  filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

#hist(Expresdata2, xlab="Intensity", col="lightgreen", main="ExpressionValue before log & normalization")

NLExp.plot <- stack(as.data.frame(c18_data_flitered_sbr))
# ggplot(NLExp.plot, aes(x=values)) +
#   geom_histogram(binwidth=20, colour="black", fill="white")
# 
# 
# den <- density(NLExp.plot$values)
# p1 <- plot(den, frame = FALSE, col = "blue",main = "Density plot")

p1 <- NLExp.plot %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'c18_original_ratio_SysBlankResuspend', theme = theme(plot.title = element_text(hjust = 0.5)))

########### Three steps normalization for metabolomics data: 1. Sample normalization -> 2. log transformation -> 3. Feature normalization
c18_colSUMr <- colSums(c18_data_flitered_sbr, dims = 1)
c18_sample_normr <- data.frame()
numeratorMr <- matrix(rep(c18_colSUMr, dim(c18_data_flitered_sbr)[1]), nrow=dim(c18_data_flitered_sbr)[1], byrow=T)
c18_sample_normr <- c18_data_flitered_sbr/numeratorMr

log_c18_sample_normr <- log2(c18_sample_normr)
#den <- density(stack(as.data.frame(c18_data_flitered_sbr))$values)
#p1 <- plot(den, frame = FALSE, col = "blue",main = "Density plot")

c18_rowmeanr <- apply(log_c18_sample_normr, 1, mean)
rowmeanMr <- t(matrix(rep(c18_rowmeanr,dim(log_c18_sample_normr)[2]), nrow = dim(log_c18_sample_normr)[2], byrow=T))
c18_sdr <- apply(log_c18_sample_normr, 1, sd)
sdMr <- t(matrix(rep(c18_sdr,dim(log_c18_sample_normr)[2]), nrow = dim(log_c18_sample_normr)[2], byrow=T))
c18_normalized_sbr <- (log_c18_sample_normr-rowmeanMr)/sdMr

NLExp.plot <- stack(as.data.frame(c18_normalized_sbr))
p1 <- NLExp.plot %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

medata <- melt(c18_normalized_sbr) %>% 
  dplyr::filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'c18_normalized_ratio_SysBlankResuspend', theme = theme(plot.title = element_text(hjust = 0.5)))



#######################################################      Stats test on the merged normalized data
# checking the mean difference between two normalized categories (water-based culture solution, particle-based culture solution)
t.test(NLExp.plot$values, NLExp.plotsb$values)
# Based on the stats test, it's good to merge the data to have better power
c18_normalized <- cbind(c18_normalized_sb,c18_normalized_sbr)
# The index of rows of c18 data can be found at 'test_c18_data_filtered3'

# check the merged normalized data distribution
NLExp.plot_me <- stack(as.data.frame(c18_normalized))
p1 <- NLExp.plot_me %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

medata_me <- melt(c18_normalized) %>% 
  dplyr::filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata_me, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'c18_normalized_ratio', theme = theme(plot.title = element_text(hjust = 0.5)))


### ANALYSES
library(umap)
library(M3C)
library(labdsv)
out_pca <- prcomp(c18_normalized)

library(ggfortify)
plot(out_pca, type="l")
autoplot(out_pca, main="C18_PCA_sampleclustering") 

# M3C::tsne(as.data.frame(c18_normalized))
dist_hi <- dist(c18_normalized, method="manhattan")
library(Rtsne)
HItsne <- labdsv::tsne(dist_hi, k=3)
plot(HItsne)
title('C18_t-sne_featureclustering')

Hiumap <- umap::umap(c18_normalized, n_components=4)
plot(Hiumap$layout) 
title('C18_umap_featureclustering')



sampleNames <- str_sub(colnames(c18_normalized),7,14)
timecate <- samples_info$`Time point`
timecate <- c(timecate[1:4],timecate[25:26],timecate[9:12],timecate[17:20],timecate[5:8],timecate[13:16],timecate[21:24])
culture_solution <- samples_info$Content
culture_solution <- c(culture_solution[1:4],culture_solution[25:26],culture_solution[9:12],culture_solution[17:20],culture_solution[5:8],culture_solution[13:16],culture_solution[21:24])

#text(x=out_pca$x[,1], y=out_pca$x[,2],labels=sampleNames)
out_pca$x[,1]
ggplot(out_pca$x, aes(x=out_pca$x[,1], y=out_pca$x[,2], group = culture_solution, colour = culture_solution)) + geom_point() +
  geom_text(label=sampleNames, hjust=0, vjust=0, check_overlap = F) + stat_ellipse(aes(fill = culture_solution))


# 3D plot for out_pca$x, HItsne$points, Hiumap$layout
library(plotly)
temp <- rnorm(length(Hiumap$layout[,1]), mean=30, sd=5)
plot_ly(x=HItsne$points[,1], y=HItsne$points[,2], z=HItsne$points[,3], type="scatter3d", mode="markers", color = temp)

# based on the 3d PCA plot, there is no outlier and no obvious pattern of feature.

colnames(c18_normalized) <- sampleNames

## check the heatmap
library(gplots)
heatmap.2(x=c18_normalized)


# Hierarchical dengegram
dist_mat <- dist(c18_normalized, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 4)
plot(cut_avg)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 4, border = 2:6)
abline(h = 6.73, col = 'red')


library(dendextend)
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, h = 3)
plot(avg_col_dend)

library(dplyr)
seeds_df_cl <- mutate(as.data.frame(c18_normalized), cluster = cut_avg)
count(seeds_df_cl,cluster)

tabulate(cut_avg)

#### Try Kmeans on the data
# Decide how many clusters to look at
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)
set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(c18_normalized, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot

km.out <- kmeans(c18_normalized, centers = 4, nstart = 20)
km.out



#### WGCNA
c18_normalized_t <- t(c18_normalized) # genes should be in the columns
# Choose a set of soft-thresholding powers
powers = c(c(1:12), seq(from = 14, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(c18_normalized_t, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
#cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 12;
adjacency = adjacency(c18_normalized_t, power = softPower,type = "signed");

TOM = TOMsimilarity(adjacency,TOMType = "signed");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "metobolites clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 4, pamRespectsDendro = FALSE,
                            minClusterSize = 15);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "metobolites dendrogram and module colors")

dynamicMods1=as.data.frame(dynamicMods)
clusteridx=list()
for (i in 1:7) {
  clusteridx[i] <- list(which(dynamicMods1$dynamicMods==i, arr.ind = T))
  print(length(clusteridx[[i]]))
}

length(which(dynamicColors=="turquoise"))
# back to check processed dataset 
dim(test_c18_data_filtered3)
cluster1 <- test_c18_data_filtered3[clusteridx[[1]],]
cluster1_compounds <- as.data.frame(na.omit(unique(cluster1$name)))
colnames(cluster1_compounds) <- c("compound")
0.05/41
getwd()
write.table(cluster1_compounds, file = "/Users/yangzhang/Downloads/MP_Bio_Package/results/cluster1_compounds.csv",sep = "\t",row.names = F, col.names = F)

print(paste("cluster",1,"_compounds.csv",  sep = ""))

for (i in 1:length(clusteridx)) {
  CLusters <- test_c18_data_filtered3[clusteridx[[i]],]
  CLusters_compounds <- as.data.frame(na.omit(unique(CLusters$name)))
  colnames(CLusters_compounds) <- c("compound")
  print(length(CLusters_compounds$compound))
  pathname <- paste("/Users/yangzhang/Downloads/MP_Bio_Package/results/","cluster", i,"_compounds.csv",  sep = "")
  print(pathname)
  write.table(CLusters_compounds, file = pathname, sep = "\t",row.names = F, col.names = F)
}


HILIC_normalized1 <- cbind(c18_normalized, as.character(dynamicMods1$dynamicMods))
colnames(HILIC_normalized1)[27] <- c("cluster")
#autoplot(out_pca, data = HILIC_normalized1, aes(fill= Cluster_wgcna, color = Cluster_wgcna))
Cluster_wgcna <- as.character(dynamicMods1$dynamicMods)

Cluster_kmeans <- as.character(km.out$cluster)
ggplot(HItsne$points, aes(x=HItsne$points[,1], y=HItsne$points[,2], group = Cluster_kmeans, colour = Cluster_kmeans)) + 
  geom_point() + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","red","black","green","blue","pink"))

# 3D plot for out_pca$x, HItsne$points, Hiumap$layout
ggplot(HItsne$points, aes(x=HItsne$points[,1], y=HItsne$points[,2], group = Cluster_kmeans, colour = Cluster_kmeans)) + 
  geom_point() + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","red","black","green","blue","pink"))


ggplot(Hiumap$layout, aes(x=Hiumap$layout[,1], y=Hiumap$layout[,2], group = Cluster_kmeans, colour = Cluster_wgcna)) + 
  geom_point() + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","red","black","green","blue","pink"))


temp <- rnorm(length(Hiumap$layout[,1]), mean=30, sd=5)
plot_ly(x=out_pca$x[,1], y=out_pca$x[,2], z=out_pca$x[,3], type="scatter3d", mode="markers", color = Cluster_wgcna)




###########################################################################################################################################################
############################################################## HILIC data           ##############################################################
############################################################################################################################


HILIC_data <- read_excel("HILIC_ChemicalFormula.xlsx", sheet = "Compounds", col_names =TRUE)

HILIC_data <- as.data.frame(HILIC_data) %>% 
  clean_names()
HILIC_x <- rownames(HILIC_data)

HILIC_x <- as.data.frame(HILIC_x)
HILIC_data_filtered <- cbind(HILIC_x, HILIC_data)

HILIC_data_filtered <- HILIC_data_filtered %>%
  filter(!is.na(checked))

colnames(HILIC_data_filtered)
test_HILIC_data_filtered1 <- HILIC_data_filtered[!is.na(HILIC_data_filtered[,'annot_delta_mass_ppm']),]
test_HILIC_data_filtered2 <- test_HILIC_data_filtered1[!is.na(test_HILIC_data_filtered1[,'rt_min']),]

hist(abs(as.numeric(test_HILIC_data_filtered2[,"annot_delta_mass_ppm"])), xlab="Intensity", col="lightgreen", main="annot_delta_mass_ppm_HILIC",breaks = 40)
x <- quantile(abs(as.numeric(test_HILIC_data_filtered2[,"annot_delta_mass_ppm"])), probs = seq(0, 1, 0.05))
x
xidx <- which(abs(as.numeric(test_HILIC_data_filtered2[,"annot_delta_mass_ppm"])) >= x[20], arr.ind = T)
view(test_HILIC_data_filtered2[xidx,])
test_HILIC_data_filtered3 <- test_HILIC_data_filtered2[-xidx,]
hist(as.numeric(test_HILIC_data_filtered3[,"annot_delta_mass_ppm"]), xlab="Intensity", col="lightgreen", main="ppm_outliersRemoved_HILIC",breaks = 40)

colnames(test_HILIC_data_filtered3)
# select ratio_omics_system blank (14 columns)
HILIC_data_flitered_sb <- mapply(test_HILIC_data_filtered3[,45:58],FUN=as.numeric)

#check nan values in the target data
nanidx <- which(HILIC_data_flitered_sb==NA, arr.ind = T)
length(nanidx)


medata <- melt(HILIC_data_flitered_sb) %>% 
  filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

#hist(Expresdata2, xlab="Intensity", col="lightgreen", main="ExpressionValue before log & normalization")

NLExp.plot <- stack(as.data.frame(HILIC_data_flitered_sb))
p1 <- NLExp.plot %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'HILIC_original_ratio_SysBlank', theme = theme(plot.title = element_text(hjust = 0.5)))

########### Three steps normalization for metabolomics data: 1. Sample normalization -> 2. log transformation -> 3. Feature normalization
# Sample normalization: sum normalization
HILIC_colSUM <- colSums(HILIC_data_flitered_sb, dims = 1)
HILIC_sample_norm <- data.frame()
HInumeratorM <- matrix(rep(HILIC_colSUM, dim(HILIC_data_flitered_sb)[1]), nrow=dim(HILIC_data_flitered_sb)[1], byrow=T)
HILIC_sample_norm <- HILIC_data_flitered_sb/HInumeratorM
# Sample normalization: quantile normalization
HILIC_quanorm_sb <- normalize.quantiles(HILIC_data_flitered_sb, copy = F)

# log transfermation
#log_c18_sample_norm <- log2(c18_sample_norm)
log_HILIC_sample_norm <- log(HILIC_quanorm_sb)
den <- density(stack(as.data.frame(log_HILIC_sample_norm))$values)
p1 <- plot(den, frame = FALSE, col = "blue",main = "Density plot")

# feature normalization
HILIC_rowmean <- apply(log_HILIC_sample_norm, 1, mean)
rowmeanM <- t(matrix(rep(HILIC_rowmean,dim(log_HILIC_sample_norm)[2]), nrow = dim(log_HILIC_sample_norm)[2], byrow=T))
HILIC_sd <- apply(log_HILIC_sample_norm, 1, sd)
sdM <- t(matrix(rep(HILIC_sd,dim(log_HILIC_sample_norm)[2]), nrow = dim(log_HILIC_sample_norm)[2], byrow=T))
HILIC_normalized_sb <- (log_HILIC_sample_norm-rowmeanM)/sdM

NLExp.plotsb <- stack(as.data.frame(HILIC_normalized_sb))
p1 <- NLExp.plotsb %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

medata_sb <- melt(HILIC_normalized_sb) %>% 
  dplyr::filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata_sb, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'HILIC_normalized_ratio_SysBlank', theme = theme(plot.title = element_text(hjust = 0.5)))



#######################################################     ratio: system blank resuspend
# select ratio_omics_system blank resuspend (12 columns)
HILIC_data_flitered_sbr <- mapply(test_HILIC_data_filtered3[,59:70],FUN=as.numeric)
HILIC_data_flitered_sbr <- HILIC_data_flitered_sbr[-352,]

#check nan values in the target data
which(is.na(HILIC_data_flitered_sbr), arr.ind = T)
length(nanidx)

# print(colnames(c18_data_flitered_sb))
# boxplot(t(c18_data_flitered_sb))

medata <- melt(HILIC_data_flitered_sbr) %>% 
  filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

#hist(Expresdata2, xlab="Intensity", col="lightgreen", main="ExpressionValue before log & normalization")

NLExp.plot <- stack(as.data.frame(HILIC_data_flitered_sbr))
p1 <- NLExp.plot %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'HILIC_original_ratio_SysBlankResuspend', theme = theme(plot.title = element_text(hjust = 0.5)))

########### Three steps normalization for metabolomics data: 1. Sample normalization -> 2. log transformation -> 3. Feature normalization
# Sample normalization: sum normalization
HILIC_colSUM <- colSums(HILIC_data_flitered_sbr, dims = 1)
HILIC_sample_norm <- data.frame()
HInumeratorM <- matrix(rep(HILIC_colSUM, dim(HILIC_data_flitered_sbr)[1]), nrow=dim(HILIC_data_flitered_sbr)[1], byrow=T)
HILIC_sample_norm <- HILIC_data_flitered_sbr/HInumeratorM
# Sample normalization: quantile normalization
HILIC_quanorm_sbr <- normalize.quantiles(HILIC_data_flitered_sbr, copy = F)

# log transfermation
#log_c18_sample_norm <- log2(c18_sample_norm)
log_HILIC_sample_norm <- log(HILIC_quanorm_sbr)
den <- density(stack(as.data.frame(log_HILIC_sample_norm))$values)
p1 <- plot(den, frame = FALSE, col = "blue",main = "Density plot")

# feature normalization
HILIC_rowmean <- apply(log_HILIC_sample_norm, 1, mean)
rowmeanM <- t(matrix(rep(HILIC_rowmean,dim(log_HILIC_sample_norm)[2]), nrow = dim(log_HILIC_sample_norm)[2], byrow=T))
HILIC_sd <- apply(log_HILIC_sample_norm, 1, sd)
sdM <- t(matrix(rep(HILIC_sd,dim(log_HILIC_sample_norm)[2]), nrow = dim(log_HILIC_sample_norm)[2], byrow=T))
HILIC_normalized_sbr <- (log_HILIC_sample_norm-rowmeanM)/sdM

HILIC_normalized_sbr <- HILIC_normalized_sbr[-352,]
NLExp.plotsbr <- stack(as.data.frame(HILIC_normalized_sbr))
p1 <- NLExp.plotsbr %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

medata_sbr <- melt(HILIC_normalized_sbr) %>% 
  dplyr::filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata_sbr, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'HILIC_normalized_ratio_SysBlankResuspend', theme = theme(plot.title = element_text(hjust = 0.5)))



#######################################################      Stats test on the merged normalized data
# checking the mean difference between two normalized categories (water-based culture solution, particle-based culture solution)
t.test(NLExp.plotsb$values, NLExp.plotsbr$values)
# Based on the stats test, it's good to merge the data to have better power
HILIC_normalized <- cbind(HILIC_normalized_sb, HILIC_normalized_sbr)
naidx <- which(is.na(HILIC_normalized), arr.ind = T)[1][1]
test_HILIC_data_filtered3 <- test_HILIC_data_filtered3[-naidx,]
dim(test_HILIC_data_filtered3)
HILIC_normalized <- na.omit(HILIC_normalized)
dim(HILIC_normalized)
# The index of rows of c18 data can be found at 'test_c18_data_filtered3'

# check the merged normalized data distribution
NLExp.plot_me <- stack(as.data.frame(HILIC_normalized))
p1 <- NLExp.plot_me %>%
  ggplot(aes(x=values)) +
  geom_density() + 
  theme_test()

medata_me <- melt(HILIC_normalized) %>% 
  dplyr::filter(0<Var1 & Var1<=40 )
p2 <- ggplot(data = medata_me, aes(x=Var1, y=value)) + geom_boxplot(aes(group=Var1)) +
  coord_flip() + theme_test()

p1+p2 + plot_layout(ncol=1, heights = c(1,3)) + plot_annotation(title = 'HILIC_normalized_ratio', theme = theme(plot.title = element_text(hjust = 0.5)))

########## checking nan values in numeric matrices of metabolomics
dim(c18_normalized)
dim(HILIC_normalized)
which(is.na(c18_normalized), arr.ind = T)
which(is.na(HILIC_normalized), arr.ind = T)


### ANALYSES
library(umap)
library(M3C)
library(labdsv)
out_pca <- prcomp(HILIC_normalized)
# M3C::tsne(as.data.frame(HILIC_normalized))
dist_hi <- dist(HILIC_normalized, method="manhattan")
library(Rtsne)
HItsne <- labdsv::tsne(dist_hi, k=3)
plot(HItsne)

Hiumap <- umap::umap(HILIC_normalized, n_components=4)
plot(Hiumap$layout)

plot(out_pca, type="l")
sampleNames <- str_sub(colnames(HILIC_normalized),7,14)
timecate <- samples_info$`Time point`
timecate <- c(timecate[1:4],timecate[25:26],timecate[9:12],timecate[17:20],timecate[5:8],timecate[13:16],timecate[21:24])
culture_solution <- samples_info$Content
culture_solution <- c(culture_solution[1:4],culture_solution[25:26],culture_solution[9:12],culture_solution[17:20],culture_solution[5:8],culture_solution[13:16],culture_solution[21:24])
library(ggfortify)
autoplot(out_pca, data = HILIC_normalized) 
#text(x=out_pca$x[,1], y=out_pca$x[,2],labels=sampleNames)
out_pca$x[,1]
ggplot(out_pca$x, aes(x=out_pca$x[,1], y=out_pca$x[,2], group = culture_solution, colour = culture_solution)) + geom_point() +
  geom_text(label=sampleNames, hjust=0, vjust=0, check_overlap = F) + stat_ellipse(aes(fill = culture_solution))


# 3D plot for out_pca$x, HItsne$points, Hiumap$layout
library(plotly)
temp <- rnorm(length(Hiumap$layout[,1]), mean=30, sd=5)
plot_ly(x=Hiumap$layout[,1], y=Hiumap$layout[,2], z=Hiumap$layout[,3], type="scatter3d", mode="markers", color = temp)
# based on the 3d PCA plot, there is no outlier and no obvious pattern of feature.

## check the heatmap
library(gplots)
colnames(HILIC_normalized) <- sampleNames
hit <- heatmap.2(x=HILIC_normalized)

# Hierarchical dengegram
dist_mat <- dist(HILIC_normalized, method = 'euclidean')
hclust_avg <- hclust(dist_mat, method = 'average')
plot(hclust_avg)
cut_avg <- cutree(hclust_avg, k = 3)
plot(cut_avg)
plot(hclust_avg)
rect.hclust(hclust_avg , k = 3, border = 2:6)
abline(h = 6.8, col = 'red')


library(dendextend)
avg_dend_obj <- as.dendrogram(hclust_avg)
avg_col_dend <- color_branches(avg_dend_obj, h = 3)
plot(avg_col_dend)

library(dplyr)
seeds_df_cl <- mutate(as.data.frame(HILIC_normalized), cluster = cut_avg)
count(seeds_df_cl,cluster)

tabulate(cut_avg)

#### Try Kmeans on the data
# Decide how many clusters to look at
n_clusters <- 10

# Initialize total within sum of squares error: wss
wss <- numeric(n_clusters)
set.seed(123)

# Look over 1 to n possible clusters
for (i in 1:n_clusters) {
  # Fit the model: km.out
  km.out <- kmeans(HILIC_normalized, centers = i, nstart = 20)
  # Save the within cluster sum of squares
  wss[i] <- km.out$tot.withinss
}

# Produce a scree plot
wss_df <- tibble(clusters = 1:n_clusters, wss = wss)

scree_plot <- ggplot(wss_df, aes(x = clusters, y = wss, group = 1)) +
  geom_point(size = 4)+
  geom_line() +
  scale_x_continuous(breaks = c(2, 4, 6, 8, 10)) +
  xlab('Number of clusters')
scree_plot

km.out <- kmeans(HILIC_normalized, centers = 6, nstart = 20)
km.out

#### WGCNA
HILIC_normalized_t <- t(HILIC_normalized) # genes should be in the columns
# Choose a set of soft-thresholding powers
powers = c(c(1:12), seq(from = 14, to=30, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(HILIC_normalized_t, powerVector = powers, verbose = 5, networkType = "signed")
# Plot the results:
#sizeGrWindow(9, 5)
#par(mfrow = c(1,2));
#cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

softPower = 20;
adjacency = adjacency(HILIC_normalized_t, power = softPower,type = "signed");

TOM = TOMsimilarity(adjacency,TOMType = "signed");
dissTOM = 1-TOM

# Call the hierarchical clustering function
geneTree = hclust(as.dist(dissTOM), method = "average");
# Plot the resulting clustering tree (dendrogram)
#sizeGrWindow(12,9)
plot(geneTree, xlab="", sub="", main = "metobolites clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04);

# Module identification using dynamic tree cut:
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = 25);
table(dynamicMods)

# Convert numeric lables into colors
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
# Plot the dendrogram and colors underneath
#sizeGrWindow(8,6)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")

dynamicMods1=as.data.frame(dynamicMods)
clusteridx_HI=list()
for (i in 0:7) {
  clusteridx_HI[i+1] <- list(which(dynamicMods1$dynamicMods==i, arr.ind = T))
  print(length(clusteridx_HI[[i+1]]))
}

length(which(dynamicColors=="turquoise"))
# back to check processed dataset 
dim(test_HILIC_data_filtered3)
cluster1 <- test_HILIC_data_filtered3[clusteridx[[1]],]
cluster1_compounds <- as.data.frame(na.omit(unique(cluster1$name)))
colnames(cluster1_compounds) <- c("compound")
0.05/41
getwd()
write.table(cluster1_compounds, file = "/Users/yangzhang/Downloads/MP_Bio_Package/results/cluster1_compounds.csv",sep = "\t",row.names = F, col.names = F)

print(paste("cluster",1,"_compounds.csv",  sep = ""))

for (i in 1:length(clusteridx_HI)) {
  CLusters <- test_HILIC_data_filtered3[clusteridx_HI[[i]],]
  CLusters_compounds <- as.data.frame(na.omit(unique(CLusters$name)))
  colnames(CLusters_compounds) <- c("compound")
  print(length(CLusters_compounds$compound))
  pathname <- paste("/Users/yangzhang/Downloads/MP_Bio_Package/results/","cluster", i,"_compounds_HILIC.csv",  sep = "")
  print(pathname)
  write.table(CLusters_compounds, file = pathname, sep = "\t",row.names = F, col.names = F)
}

HILIC_normalized1 <- cbind(HILIC_normalized, as.character(dynamicMods1$dynamicMods))
colnames(HILIC_normalized1)[27] <- c("cluster")
#autoplot(out_pca, data = HILIC_normalized1, aes(fill= Cluster_wgcna, color = Cluster_wgcna))
Cluster_wgcna <- as.character(dynamicMods1$dynamicMods)

Cluster_kmeans <- as.character(km.out$cluster)
ggplot(out_pca$x, aes(x=out_pca$x[,1], y=out_pca$x[,2], group = Cluster_kmeans, colour = Cluster_kmeans)) + 
  geom_point() + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","red","black","green","blue","pink"))

# 3D plot for out_pca$x, HItsne$points, Hiumap$layout
ggplot(HItsne$points, aes(x=HItsne$points[,1], y=HItsne$points[,2], group = Cluster_kmeans, colour = Cluster_kmeans)) + 
  geom_point() + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","red","black","green","blue","pink"))


ggplot(Hiumap$layout, aes(x=Hiumap$layout[,1], y=Hiumap$layout[,2], group = Cluster_kmeans, colour = Cluster_wgcna)) + 
  geom_point() + scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9","red","black","green","blue","pink"))


temp <- rnorm(length(Hiumap$layout[,1]), mean=30, sd=5)
plot_ly(x=HItsne$points[,1], y=HItsne$points[,2], z=HItsne$points[,3], type="scatter3d", mode="markers", color = Cluster_wgcna)

#+
#scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9","red","black","green","blue","pink")) +
# geom_text(label=sampleNames, hjust=0, vjust=0, check_overlap = F) + stat_ellipse(aes(fill = culture_solution))



#################################################################################### SNFtool: sample integration

K = 20; ##number of neighbors, must be greater than 1. usually (10~30)
alpha = 0.5; ##hyperparameter, usually (0.3~0.8)
T = 20; ###Number of Iterations, usually (10~50)

#create data list
data_list <- list(c18_normalized, HILIC_normalized)
names(data_list) <- c("C18", "HILIC")

## Create distance matrices for each data type
data_dist_list <- lapply(X = data_list,
                         function(x){dist2(t(x),t(x))})

## Create affinity matrices for each data type
data_aff_list <-  lapply(X = data_dist_list,
                         function(x){affinityMatrix(x,K,alpha)})
W = SNF(data_aff_list,K,T)
estimateNumberOfClustersGivenGraph(W = W,2:10)
clustering2 = spectralClustering(W,2)
table(clustering2)
barplot(table(clustering2),col = c('darkorchid4','dodgerblue4'))
displayClustersWithHeatmap(W = W, group = clustering2)



library(cutpointr)
rocdata <- as.data.frame(cbind(clustering2, clustering2))
colnames(rocdata) <- c("media", "cluster")
roc_curve <- roc(data = rocdata, x = media, class = clustering2, pos_class = "15% FBS SMSC Media", neg_class = "Serum Free Media")
plot_roc(roc_curve)
mcp <- multi_cutpointr(data = rocdata, x = clustering2, class = media, pos_class = "15% FBS SMSC Media", neg_class = "Serum Free Media")
library(starter)

