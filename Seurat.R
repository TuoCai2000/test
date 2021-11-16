install.packages("Seurat")
library(Seurat)
library(dplyr)
library(tidyverse)
#载入数据
pbmc.data <- Read10X(data.dir = "/Users/caita/Downloads/filtered_gene_bc_matrices/hg19")
#初始化Seurat object，设置过滤指标
pbmc <- CreateSeuratObject(counts = pbmc.data,project = "pbmc3k",min.cells = 3,
                           min.features = 200)
#counts:未标准化的数据，如原始计数；project：设置Seurat对象的项目名称；
#min.cells：包含至少在这些细胞检测到的features；min.features：包含至少检测到这些features的细胞
#过滤检测少于200个基因的细胞和少于3个细胞检测出的基因
pbmc
#1个数据集，包含2700个细胞，以及13714个基因
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc,pattern = "MT-")
#查看三个基因的前三十行矩阵
pbmc.data[c("CD3D","TCL1A","MS4A1"),1:30]
#大多数表达量为0，节省内存和速度
#展示前五个细胞的QC指标
head(pbmc@meta.data,5)
#用小提琴图可视化QC指标
VlnPlot(pbmc,features = c("nFeature_RNA","nCount_RNA","percent.mt"),ncol = 3)
#基因数、count数和线粒体比例相关性
plot1 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "percent.mt")
plot2 <- FeatureScatter(pbmc,feature1 = "nCount_RNA",feature2 = "nFeature_RNA")
#合并两图
CombinePlots(plots = list(plot1,plot2))
#数据再过滤和标准化
pbmc <- subset(pbmc,subset=nFeature_RNA>200 & nFeature_RNA<2500 &percent.mt<5)
#筛选出200<nFeature_RNA<2500和percent.mt<5的数据
#数据标准化
pbmc <- NormalizeData(pbmc,normalization.method = "LogNormalize",scale.factor = 10000)
#找出高可变基因
pbmc <- FindVariableFeatures(pbmc,selection.method = "vst",nfeatures = 2000)
#查看最高变的10个基因
top10 <- head(VariableFeatures(pbmc),10)
#画出不带标签或带标签基因点图
plot1 <- VariableFeaturePlot(pbmc)
plot2 <- LabelPoints(plot = plot1,points = top10,repel = TRUE)
plot1+plot2
#数据标准化（缩放）
all.genes <- rownames(pbmc)
pbmc <- ScaleData(pbmc,features = all.genes)
pbmc <- ScaleData(pbmc,vars.to.regress = "percent.mt")
#PCA线性降维
pbmc <- RunPCA(pbmc,features = VariableFeatures(object = pbmc))
#查看PC1～PC5的前五个基因
print(pbmc[["pca"]],dims=1:5,nfeatures = 5)
#查看特定PC中的基因和基因对应的贡献度
VizDimLoadings(pbmc,dims = 1:2,reduction = "pca")
#绘制PCA散点图
DimPlot(pbmc,reduction = "pca")
# 绘制PCA热图
DimHeatmap(pbmc,dims = 1:15,cells = 500,balanced = TRUE)
#确定数据的维度
pbmc <- JackStraw(pbmc,num.replicate = 100)
pbmc <- ScoreJackStraw(pbmc,dims = 1:20)
JackStrawPlot(pbmc,dims = 1:15)
ElbowPlot(pbmc)
#确定选择10个主成分用于后续分析
#细胞聚类
pbmc <- FindNeighbors(pbmc,dims = 1:10)
pbmc <- FindClusters(pbmc,resolution = 0.5)
#查看前五个细胞的分类ID
head(Idents(pbmc),5)
#非线性降维
#UMAP
pbmc <- RunUMAP(pbmc,dims = 1:10)
DimPlot(pbmc,reduction = "umap")
#显示聚类标签
DimPlot(pbmc,reduction = "umap",label = TRUE)
#t-SNE
pbmc <- RunTSNE(pbmc,dims = 1:10)
DimPlot(pbmc,reduction = "tsne",label = TRUE)
#保存结果
saveRDS(pbmc,file = "/Users/caita/Downloads/filtered_gene_bc_matrices/pbmc_tutorial.rds")
#找差异表达基因
#cluster1的标记基因
cluster1.markers <- FindMarkers(pbmc,ident.1 = 1,min.pct = 0.25)
head(cluster1.markers,5)
cluster0.markers <- FindMarkers(pbmc,ident.1 = 0,min.pct = 0.25)
head(cluster0.markers,5)
#找出区分cluster5与cluster0和cluster3的所有标记基因
cluster5.markers <- FindMarkers(pbmc,ident.1 = 1,ident.2 = c(0,3),
                                min.pct = 0.25)
head(cluster5.markers,5)
#找出全部聚类的marker基因 
pbmc.markers <- FindAllMarkers(pbmc,only.pos = TRUE,min.pct = 0.25,
                               logfc.threshold = 0.25)
pbmc.markers %>% group_by(cluster) %>% top_n(n=2,wt=avg_log2FC)
#小提琴图展示marker基因在各个聚类中的表达量
VlnPlot(pbmc,features = c("MS4A1","CD79A"))
#将marker基因的表达量映射在聚类降维结果中
FeaturePlot(pbmc,features = c("MS4A1","GNLY","CD3E","CD14","FCER1A","FCGR3A",
                              "LYZ","PPBP","CD8A"))
#热图
top10 <- pbmc.markers %>% group_by(cluster) %>% top_n(n=10,wt=avg_log2FC)
DoHeatmap(pbmc,features = top10$gene)+NoLegend()
#细胞类型替换和结果重新绘制
new.cluster.ids <- c("Naive CD4 T","CD14+ Mono","Memory CD4 T","B","CD8 T",
                     "FCGR3A+ Mono","NK","DC","Platelet")
names(new.cluster.ids) <- levels(pbmc)
pbmc <- RenameIdents(pbmc,new.cluster.ids)
DimPlot(pbmc,reduction = "umap",label = TRUE,pt.size = 0.5)+NoLegend()
saveRDS(pbmc,file="/Users/caita/Downloads/filtered_gene_bc_matrices/pbmc3k_final.rds")
