#批次效应
BiocManager::install("TENxPBMCData")
library(TENxPBMCData) #由三个pbmc组成的数据集
all.sce <- list(
  pbmc3k=TENxPBMCData('pbmc3k'),
  pbmc4k=TENxPBMCData('pbmc4k')
)
all.sce
unfiltered <- all.sce #备份未过滤的数据
BiocManager::install('scater')
library(scater)
stats <- high.mito <- list() #新建2个空列表stats和high.mito，放质控结果
#names(all.sce)是pbmc.3k和pbmc.4k
all.sce[[2]]
for (n in names(all.sce)){
  current <- all.sce[[n]] #pbmc.3k和pbmc.4k的信息
  is.mito <- grep('MT',rowData(current)$Symbol_TENx) #Symbol_TENx是所有基因列表 筛选出所有线粒体基因
  stats[[n]] <- perCellQCMetrics(current,subsets=list(Mito=is.mito)) #细胞水平指标 计算该数据矩阵的一些属性,来识别和删除可能存在问题的单元格
  high.mito[[n]] <- isOutlier(stats[[n]]$subsets_Mito_percent,type = 'higher') #查找上面数据矩阵的离群值放入high.mito TRUE/FALSE
  all.sce[[n]] <- current[,!high.mito[[n]]] #过滤线粒体基因比例较高的细胞
}
#查看线粒体过滤的结果
lapply(high.mito,summary)
#质控图
qcplots <- list()
for (n in names(all.sce)){
  current <- unfiltered[[n]] #未过滤数据
  colData(current) <- cbind(colData(current),stats[[n]]) 
  current$discard <- high.mito[[n]]
  qcplots[[n]] <- plotColData(current,x='sum',y='subsets_Mito_percent',
                              colour_by = 'discard')+scale_x_log10()
}
dev.new()
do.call(gridExtra::grid.arrange,c(qcplots,ncol=2))
#批量归一化
all.sce <- lapply(all.sce,logNormCounts) #常规方法可以不用计算size.factor，直接用logNormCounts
#如果使用去卷积方法，还是要计算size.factor
lapply(all.sce,function(x) summary(sizeFactors(x)))
#批量找高变异基因
BiocManager::install('scran')
library(scran)
all.dec <- lapply(all.sce,modelGeneVar)
all.hvgs <- lapply(all.dec,getTopHVGs,prop=0.1) #筛选高变异基因
#作图
par(mfrow=c(1,3))
for (n in names(all.dec)) {
  curdec <- all.dec[[n]]
  plot(curdec$mean,curdec$total,pch=16,cex=0.5,main=n,
       xlab='Mean of log-expression',ylab='Variance of log-expression')
  curfit <- metadata(curdec)
  curve(curfit$trend(x),col='dodgerblue',add = TRUE,lwd=2)
}
#每个点表示一个基因，蓝线代表技术因素导致的偏差，纵坐标表示总偏差：技术偏差+生物因素偏差
#衡量一个基因的生物因素偏差大小，看纵坐标减去对应蓝线的值
#批量降维
#选择与PCA近似的SVD算法（奇异值分解）对稀疏矩阵SVD算法更实用（节省空间）
library(BiocSingular)
set.seed(10000)
all.sce <- mapply(FUN = runPCA,x=all.sce,subset_row=all.hvgs,
                  MoreArgs = list(ncomponents=25,BSPARAM=RandomParam()),
                  SIMPLIFY = FALSE)
set.seed(100000)
all.sce <- lapply(all.sce,runTSNE,dimred='PCA')
set.seed(1000000)
all.sce <- lapply(all.sce,runUMAP,dimred='PCA')
#批量聚类（最近邻）
for (n in names(all.sce)) {
  g <- buildSNNGraph(all.sce[[n]],k=10,use.dimred='PCA')
  clust <- igraph::cluster_walktrap(g)$membership
  colLabels(all.sce[[n]]) <- factor(clust)
}
#查看各自分了多少群
lapply(all.sce,function(x) table(colLabels(x)))
###去除批次效应
#取测序基因交集
universe <- intersect(rownames(all.sce[['pbmc3k']]),rownames(all.sce[['pbmc4k']])) #取测序基因的交集
length(universe) #有多少基因在两份样本中共同表达
#将两份样本中都表达的基因提取出来
pbmc3k <- all.sce[['pbmc3k']][universe,]
pbmc4k <- all.sce[['pbmc4k']][universe,]
dec3k <- all.dec[['pbmc3k']][universe,]
dec4k <- all.dec[['pbmc4k']][universe,]
#使两个批次的测序深度相统一
BiocManager::install('batchelor')
library(batchelor)
#multiBatchNorm()校正因不同批次的测序深度不同造成的表达水平差异，消除部分技术误差带来的差异
rescaled <- multiBatchNorm(pbmc3k,pbmc4k)
pbmc3k <- rescaled[[1]]
pbmc4k <- rescaled[[2]]
#合并数据集
rowData(pbmc3k) <- rowData(pbmc4k)
pbmc3k$batch <- '3k' #添加batch属性
pbmc4k$batch <- '4k'
uncorrected <- cbind(pbmc3k,pbmc4k) #合并数据集
library(scran)
combined.dec <- combineVar(dec3k,dec4k) #combineVar在考虑批次因素前提下，筛选出尽可能全面的hvg
chosen.hvgs <- combined.dec$bio>0 #放宽筛选标准
sum(chosen.hvgs) #筛选出13431个hvg
#PCA降维
library(scater)
set.seed(0010101010)
uncorrected <- runPCA(uncorrected,subset_row=chosen.hvgs,
                      BSPARAM=BiocSingular::RandomParam())
#KNN聚类
library(scran)
snn.gr <- buildSNNGraph(uncorrected,use.dimred='PCA')
clusters <- igraph::cluster_walktrap(snn.gr)$membership
tab <- table(Cluster=clusters,Batch=uncorrected$batch)
tab #KNN结果可以看出批次效应十分明显
#t-SNE可视化
set.seed(1111001)
uncorrected <- runTSNE(uncorrected,dimred='PCA')
plotTSNE(uncorrected,colour_by='batch')
#基于线性回归消除批次效应
library(limma)
rescaled <- removeBatchEffect(assay(uncorrected),batch = uncorrected$batch,design = chosen.hvgs) 
#对校正后的结果再次进行批次效应评判
set.seed(1010101010)
rescaled <- runPCA(rescaled,subset_row=chosen.hvgs,
                   exprs_values='corrected',
                   BSPARAM=BiocSingular::RandomParam())
snn.gr <- buildSNNGraph(rescaled,use.dimred='PCA')
clusters.resc <- igraph::cluster_walktrap(snn.gr)$membership
tab.resc <- table(Cluster=clusters.resc,Batch=rescaled$batch)
tab.resc
#基于线性回归校正批次效应的聚类结果
rescaled <- runTSNE(rescaled,dimred='PCA')
str(rescaled)
rescaled$batch <- factor(rescaled$batch)
plotTSNE(rescaled,colour_by='batch') #仍然存在batch-specific cluster,提示批次效应未完全去除
#batch-specific cluster可能是由于批次造成的异质性，也可能是每个batch本身具有specific subpopulation
#基于MNN消除批次效应
set.seed(1000101001)
mnn.out <- fastMNN(pbmc3k,pbmc4k,d=50,k=20,subset.row = chosen.hvgs,
                   BSPARAM = BiocSingular::RandomParam(deferred = TRUE))
mnn.out
#对校正后的结果再次进行批次效应评判
library(scran)
snn.gr <- buildSNNGraph(mnn.out,use.dimred='corrected')
clusters.mnn <- igraph::cluster_walktrap(snn.gr)$membership
tab.mnn <- table(Cluster=clusters.mnn,Batch=mnn.out$batch)
tab.mnn
#基于线性回归校正批次效应的聚类结果
library(scater)
set.seed(0010101010)
mnn.out <- runTSNE(mnn.out,dimred='corrected')
mnn.out$batch <- factor(mnn.out$batch)
plotTSNE(mnn.out,colour_by='batch') #从图中直观地看，MNN效果更好

#同一个cluster里不同batch来源的细胞组成比例
clu.prop <- colSums(tab.mnn)/sum(tab.mnn) #batch1细胞数：batch2细胞数
clu.results <- apply(tab.mnn,1,FUN=chisq.test,p=clu.prop) #进行卡方检验
p.values <- vapply(clu.results,'[[',i='p.value',0)
p.values #cluster1、3、5、7、10、16的p值不够显著

#校正后的结果相比原始batch的异质性是否丢失
#校正合并多个批次后，由于细胞数目增加，聚类分群数也会增加
#batch1分为10个cluster，合并batch2后，batch1被分为15个cluster。
#需要比较对batch1而言，分为10个cluster和15个cluster的相似性程度。ARI指标???
#总体ARI指标（越接近1越理想）整个batch
library(bluster)
#ri3k <- pairwiseRand(clusters.mnn[rescaled$batch==1],colLabels(pbmc3k),mode = 'index')
#无法进行，先对前面进行补充：批量降维、聚类
ri3k <- pairwiseRand(clusters.mnn[rescaled$batch==1],colLabels(all.sce[[1]]),mode = 'index')
ri3k
ri4k <- pairwiseRand(clusters.mnn[rescaled$batch==2],colLabels(all.sce[[2]]),mode = 'index')
ri4k
#两两cluster之间的ARI
dev.off()
dev.new()
tab <- pairwiseRand(colLabels(all.sce[[1]]),clusters.mnn[rescaled$batch==1]) #计算两两cluster之间的相关性
heat3k <- pheatmap(tab,cluster_rows = FALSE,cluster_cols = FALSE,
                   col=rev(viridis::magma(100)),main = 'PBMC3K probabilities',silent = T)
heat3k #方块颜色越深越理想。
#对角线上的方格越深表示同一cluster的cell pair再聚类之后仍属于同一cluster，非对角线上的方格颜色越深表示原本不属于同一个cluster的cell pair重新聚类后仍不属于同一cluster
#校正后，单细胞数据具有批次间表达水平统一性，不建议用校正后的表达值进行差异分析，应该使用原始的表达水平并设置批次参数
m.out <- findMarkers(uncorrected,clusters.mnn,block=uncorrected$batch,
                     direction='up',lfc=1,row.data=rowData(uncorrected)[,3,drop=FALSE])
m.out[[1]] #取第一个cluster的markergene
demo <- m.out[['1']]
as.data.frame(demo[1:20,c("Symbol","Top","p.value","FDR")])
plotExpression(uncorrected,x=I(factor(clusters.mnn)),
               features = 'ENSG00000019582',colour_by = 'batch')+facet_wrap(~colour_by)
#细胞亚群鉴定基于校正批次效应后的表达矩阵进行，再对各个批次的原始数据单独加以分析判断
