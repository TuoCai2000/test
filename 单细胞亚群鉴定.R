install.packages("devtools")
library(devtools)
install_github("ggjlab/scHCL")
install_github("ggjlab/scMCA")
install.packages("pheatmap")
install.packages("shinythemes")
library(pheatmap)
library(shinythemes)
library(Seurat)
library(scHCL)
library(tidyverse)
#直接使用UMAP或t-SNE降维聚类后保存的rds文件
rds <- readRDS("/Users/caita/Downloads/filtered_gene_bc_matrices/pbmc_tutorial.rds")
counts <- as.matrix(rds@assays$RNA@counts)
dim(counts)
mca_result <- scHCL(scdata = counts,numbers_plot = 3)
str(mca_result)
#展示结果
#shinny的展示形式
scHCL_vis(mca_result)
#换一种方法
#定义gettissue函数
gettissue <- function(x,Num=3){
  top_markers <- order(x,decreasing = TRUE)[1:Num]
  return(top_markers)
}
cors <- mca_result$cors_matrix #结果的相关系数矩阵
view(cors)
cors_index <- apply(cors,2,gettissue) #对列使用gettissue函数，
#返回最大的三个值的索引(这个碱基序列对应的最有可能的三种细胞类型)
cors_index <-sort(unique(as.integer(cors_index)))#将索引转化为整数型，去除重复，将所有索引值从小到大排序
data=cors[cors_index,]#去除无用的索引，生成新矩阵
pbmc@meta.data$stim <- factor("normal")
annotation_col <- data.frame(Celltype=rds@meta.data$seurat_clusters)
rownames(annotation_col) <- rownames(rds@meta.data)
order_cells <- annotation_col %>% dplyr::mutate(.,barcod=rownames(.)) %>% dplyr::arrange(.,Celltype)
gaps_col <- c()
m <- 0
for (i in 1:max(as.numeric(annotation_col[,c("Celltype")]))) {
  gaps_col <- c(gaps_col,table(annotation_col[,c("Celltype")])[i]+m)
  m <- table(annotation_col[,c("Celltype")])[i]+m
}
library(pheatmap)
library(RColorBrewer)
pheatmap(data[1:76,as.vector(order_cells$barcod)],show_rownames = T,
         show_colnames = F,cluster_cols = F,cluster_rows = F,
         annotation_col = annotation_col,treeheight_col = annotation_col,
         treeheight_row = 0,cellheight = 8,fontsize = 8,
         cellwidth = 420/length(rownames(annotation_col)),
         border=FALSE,gaps_col = gaps_col,
         color = colorRampPalette(c("gray","white","red"))(100),main = "")


#SingleR
if(!requireNamespace("BiocManager",quietly = TRUE))
install.packages("BiocManager")  
BiocManager::install("SingleR")
library(SingleR)
install.packages("celldex")
library(celldex)
install.packages("remotes")
library(remotes)
remotes::install_github("LTLA/celldex")
ref_data <- HumanPrimaryCellAtlasData()#参考数据集
ref_data
install.packages("scRNAseq")
library(celldex)
library(scRNAseq)
install.packages("SingleCellExperiment")
querry <- as.matrix(rds@assays$RNA@data)
pred.hesc <- SingleR(test = querry,ref = ref_data,labels = ref_data$label.fine)
pred.hesc
annotation_col <- data.frame(Celltype=rds@meta.data$seurat_clusters)
rownames(annotation_col) <- rownames(rds@meta.data)
order_cells <- annotation_col %>% dplyr::mutate(.,barcod=rownames(.)) %>% dplyr::arrange(.,Celltype)
pred.hesc$clusters <- annotation_col[as.vector(rownames(pred.hesc)),]$Celltype
pred.hesc$clusters <- order_cells$Celltype
pred.hesc[as.vector(order_cells$barcod),]
pdf('/Users/caita/Downloads/filtered_gene_bc_matrices/stim_celltype.pdf',w=12,h=8)
p <- plotScoreHeatmap(pred.hesc,annotation_col=annotation_col,cells.order=pred.hesc$clusters)
print(p)
dev.off()