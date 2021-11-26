#安装Bioconductor
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install()
BiocManager::install(c("GenomicFeatures", "AnnotationDbi","BiocGenerics"))
#安装monocle2
install.packages("devtools")
devtools::install_github("cole-trapnell-lab/monocle-release@develop")
install.packages("parallel")
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("monocle")
#创建CellDataSet
library(monocle)
#导入注释好的Seurat对象
pbmc <- readRDS("/Users/caita/Downloads/filtered_gene_bc_matrices/pbmc3k_final.rds")
#提取表达矩阵信息
expr_matrix <- as(as.matrix(pbmc@assays$RNA@counts),"sparseMatrix")
#提取表型信息
p_data <- pbmc@meta.data
p_data$celltype <- pbmc@active.ident#整合每个细胞的细胞鉴定信息
#提取基因信息
f_data <- data.frame(gene_short_name=row.names(pbmc),row.names = row.names(pbmc))
#构建CDS对象
pd <- new("AnnotatedDataFrame",data=p_data)
fd <- new("AnnotatedDataFrame",data=f_data)
#将p_data和f_data从data.frame转换成AnnotatedDataFrame对象
cds <- newCellDataSet(expr_matrix,
                      phenoData = pd,
                      featureData = fd,
                      lowerDetectionLimit = 0.5,
                      expressionFamily = negbinomial.size())
#估计size factor和离散度
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)
#过滤低质量的细胞
#（Seurat通过基因表达量过滤细胞，这里通过表达某基因的细胞数目对基因进行过滤）
cds <- detectGenes(cds,min_expr = 0.1) #这一操作在fData(cds)中添加一列num_cells_expressed
print(head(fData(cds)))
expressed_genes <- row.names(subset(fData(cds),num_cells_expressed>=10))
#选择定义过程的基因
#使用Seurat选择的高变基因
library(Seurat)
express_genes <- VariableFeatures(pbmc)
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)
#使用clusters差异表达基因
deg.cluster <- FindAllMarkers(pbmc)
express_genes <- subset(deg.cluster,p_val_adj<0.05)$gene
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)
#使用monocle选择高变基因
disp_table <- dispersionTable(cds)
disp.genes <- subset(disp_table,mean_expression>=0.1 & dispersion_empirical>=1*dispersion_fit)$gene_id
cds <- setOrderingFilter(cds,express_genes)
plot_ordering_genes(cds)
#使用dpFeature
#expressed_genes来自过滤低质量细胞，也可输入Seurat筛选出的高变基因
diff <- differentialGeneTest(cds[expressed_genes,],fullModelFormulaStr = "~celltype",cores = 1)
head(diff)
#筛选出差异表达基因
deg <- subset(diff,qval<0.01)
deg <- deg[order(deg$qval,decreasing = F),]
head(deg)
#差异基因的结果文件保存
write.table(deg,file = "/Users/caita/Downloads/filtered_gene_bc_matrices/train.monocle.DEG.xls",
            col.names = T,row.names = T,sep = "\t",quote = F)
#轨迹构建基因可视化
ordergene <- rownames(deg)
cds <- setOrderingFilter(cds,ordergene)
pdf("/Users/caita/Downloads/filtered_gene_bc_matrices/train.ordergenes.pdf")
plot_ordering_genes(cds)
dev.off()
#降维
cds <- reduceDimension(cds,max_components = 2,method="DDRTree")
#拟时间轴轨迹构建和在拟时间内排列细胞
cds <- orderCells(cds)
#以pseudotime值上色
pdf("/Users/caita/Downloads/filtered_gene_bc_matrices/train.monocle.pseudotime.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "Pseudotime",size=1,show_backbone = TRUE)
dev.off()
#以细胞类型上色
pdf("/Users/caita/Downloads/filtered_gene_bc_matrices/train.monocle.celltype.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "celltype",size=1,show_backbone = TRUE)
dev.off()
#以细胞状态上色
pdf("/Users/caita/Downloads/filtered_gene_bc_matrices/train.monocle.state.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "State",size=1,show_backbone = TRUE)
dev.off()
#以seurat cluster排列细胞
pdf("/Users/caita/Downloads/filtered_gene_bc_matrices/seurat.clusters.pdf",width = 7,height = 7)
plot_cell_trajectory(cds,color_by = "seurat_clusters",size=1,show_backbone = TRUE)
dev.off()
#以细胞状态上色（拆分）“分面”轨迹图
pdf("/Users/caita/Downloads/filtered_gene_bc_matrices/train.monocle.state.faceted.pdf",width = 10,height = 7)
plot_cell_trajectory(cds,color_by = "State")+facet_wrap("~State",nrow = 1)
dev.off()
#指定基因的可视化
#选择前四个top基因并将其对象取出
keygenes <- head(ordergene,4)
cds_subset <- cds[keygenes,]
#可视化：以state/celltype/pseudotime进行
p1 <- plot_genes_in_pseudotime(cds_subset,color_by = "State")
p2 <- plot_genes_in_pseudotime(cds_subset,color_by = "celltype")
p3 <- plot_genes_in_pseudotime(cds_subset,color_by = "Pseudotime")
plotc <- p1|p2|p3
ggsave("/Users/caita/Downloads/filtered_gene_bc_matrices/Genes_pseudotimeplot.pdf",plot = plotc,width=16,height = 8)
#指定基因
s.genes <- c("SELL","CCR7","IL7R","CD84","CCL5","S100A4")
p1 <- plot_genes_jitter(cds[s.genes,],grouping = "State",color_by = "State")
p2 <- plot_genes_violin(cds[s.genes,],grouping = "State",color_by = "State")
p3 <- plot_genes_in_pseudotime(cds[s.genes,],color_by = "State")
plotc <- p1|p2|p3
ggsave("/Users/caita/Downloads/filtered_gene_bc_matrices/Genes_Jitterplot.pdf",plot=plotc,width=16,height=8)
#寻找拟时相关的基因
#把排序基因ordergene提取出来做回归分析，来找它们是否跟拟时间有显著的关系
#如果不设置，就会用所有基因来做它们与拟时间的相关性
Time_diff <- differentialGeneTest(cds[ordergene,],cores = 1,
                                  fullModelFormulaStr = "~sm.ns(Pseudotime)")
Time_diff <- Time_diff[,c(5,2,3,4,1,6,7)]#更改列的顺序，把gene放前面，也可以不改。
write.csv(Time_diff,"/Users/caita/Downloads/filtered_gene_bc_matrices/Time_diff_all.csv",
          row.names = F)
#拟时差异基因热图绘制（提取了前100个）
library(tidyverse)
Time_genes <- top_n(Time_diff,n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character()
p=plot_pseudotime_heatmap(cds[Time_genes,],num_clusters = 4,show_rownames = T,return_heatmap = T)
ggsave("/Users/caita/Downloads/filtered_gene_bc_matrices/Time_heatmap.pdf",p,width = 5,height = 10)
#显著差异基因按热图结果排序并保存
hp.genes <- p$tree_row$labels[p$tree_row$order]
Time_diff_sig <- Time_diff[hp.genes,c("gene_short_name","pval","qval")]
write.csv(Time_diff_sig,"/Users/caita/Downloads/filtered_gene_bc_matrices/Time_diff_sig.csv",row.names = F)
#也可以手动选择基因来绘制热图
marker_genes <- row.names(subset(fData(cds),
                                 gene_short_name %in% c("MEF2C","MEF2D","MYF5",
                                                        "ANPEP","PDGFRA","MYOG",
                                                        "TPM1","TPM2","MYH2",
                                                        "MYH3","NCAM1","TNNT1",
                                                        "TNNT2","TNNC1","CDK1",
                                                        "CDK2","CCNB1","CCNB2",
                                                        "CCND1","CCNA1","ID1")))
diff_test_res <- differentialGeneTest(cds[marker_genes,],
                                      fullModelFormulaStr = "~sm.ns(Pseudotime)")
sig_gene_names <- row.names(subset(diff_test_res,qval<0.1))
plot_pseudotime_heatmap(cds[sig_gene_names,],
                        num_clusters = 6,
                        cores = 1,
                        show_rownames = T)
#单细胞轨迹的“分支”分析
#之前找的拟时相关的基因是全局的（拟时起点和终点相关的基因），这一步找和分叉点相关的基因
plot_cell_trajectory(cds,color_by = "State")
#BEAM进行统计分析
BEAM_res <- BEAM(cds[ordergene,],branch_point = 1,cores = 2)
BEAM_res <- BEAM_res[order(BEAM_res$qval),]
BEAM_res <- BEAM_res[,c("gene_short_name","pval","qval")]
plot_genes_branched_heatmap(cds[row.names(subset(BEAM_res,qval<1e-4)),],
                            branch_point = 1,
                            num_clusters = 4,
                            cores = 1,
                            use_gene_short_name = T,
                            show_rownames = T)
#选前100个基因进行可视化
BEAM_genes <- top_n(BEAM_res,n=100,desc(qval)) %>% pull(gene_short_name) %>% as.character()
p <- plot_genes_branched_heatmap(cds[BEAM_genes,],branch_point = 2,num_clusters = 3,
                                 show_rownames = T,return_heatmap = T)
ggsave("/Users/caita/Downloads/filtered_gene_bc_matrices/BEAM_heatmap.pdf",p$ph_res,width=6.5,height=10)
#差异显著基因按热图结果排序并保存
hp.genes <- p$ph_res$tree_row$labels[p$ph_res$tree_row$order]
BEAM_sig <- BEAM_res[hp.genes,c("gene_short_name","pval","qval")]
write.csv(BEAM_sig,"/Users/caita/Downloads/filtered_gene_bc_matrices/BEAM_sig.csv",row.names = F)
#选择热图赏/感兴趣的基因进行可视化
genes <- row.names(subset(fData(cds),gene_short_name %in% c("CCL4","KLRD1","XCL2")))
plot_genes_branched_pseudotime(cds[genes,],branch_point = 1,
                               color_by = "State",ncol = 1)