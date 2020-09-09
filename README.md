# single-cell-R-analysis-of-17-testis

# all 17 single cell data were first merged into one seurat object

#质量控制
allsample_combined_cca <- subset(x = allsample_combined_cca, subset = nFeature_RNA > 500 & nFeature_RNA < 9000 & percent.mt < 40 & nCount_RNA <80000)

#批次效应处理
allsample_combined_cca_list <- SplitObject(allsample_combined_cca, split.by = "tech")

for (i in 1:length(allsample_combined_cca_list)) {allsample_combined_cca_list[[i]] <- NormalizeData(allsample_combined_cca_list[[i]], verbose = FALSE) allsample_combined_cca_list[[i]] <- FindVariableFeatures(allsample_combined_cca_list[[i]], selection.method = "vst",nfeatures = 1000, verbose = FALSE)}

allsample_combined_cca_anchors <- FindIntegrationAnchors(object.list = allsample_combined_cca_list, dims = 1:30,anchor.features = 1000)

allsample_combined_cca <- IntegrateData(anchorset = allsample_combined_cca_anchors, dims = 1:30)

#标准化
#DefaultAssay(allsample_combined_cca) <- "RNA"

allsample_combined_cca <- ScaleData(allsample_combined_cca)

allsample_combined_cca <- RunPCA(allsFindClustersFindClustersample_combined_cca, npcs = 30)

#聚类
allsample_combined_cca <- FindNeighbors(object = allsample_combined_cca, dims = 1:30)

allsample_combined_cca <- FindClusters(object = allsample_combined_cca, resolution = 0.2)
#降维
allsample_combined_cca <- RunUMAP(object = allsample_combined_cca, dims = 1:30,n.neighbors = 35L,seed=17)

allsample_combined_cca <- RunTSNE(object = allsample_combined_cca, dims = 1:30,n.neighbors = 35L,seed=17)

# Visualization
DimPlot(object = allsample_combined_cca, reduction = "tsne", group.by = "age",label.size = 8,pt.size = 0.1)+theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())

DimPlot(object = allsample_combined_cca, reduction = "umap", group.by = "agetype",label.size = 8,pt.size = 1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())

DimPlot(object = allsample_combined_cca, reduction = "tsne",  label = F,label.size = 8,pt.size = 0.1)+
theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())

#代谢分析
#读取基因
X_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/chr_X_gene.csv",header = F,sep = "\t",,stringsAsFactors=F)

Y_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/chr_Y_gene.csv",header = F,sep = "\t",,stringsAsFactors=F)

eXi_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/X染色体失活逃逸基因.csv",header = F,sep = "\t",,stringsAsFactors=F)

X_gene<-X_gene[,1]

Y_gene<-Y_gene[,1]

eXi_gene<-eXi_gene[,1]

OXPHOS_genes<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/氧化磷酸化-无线粒体基因 go 基因列表.csv")

glycolysis_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/经典糖酵解 go 基因列表.csv")

triglyceride_metabolic_gene<-read.csv("F:/学习/数据分析/single cell sequence/注释文件/human/甘油三酯代谢.csv")

OXPHOS_genes<-as.character(OXPHOS_genes[,1])

glycolysis_gene<-as.character(glycolysis_gene[,1])

triglyceride_metabolic_gene<-as.character(triglyceride_metabolic_gene[,1])

#计算X，Y染色体基因表达，注意这里使用的genelist应该都要在object的RNA@meta.features里面存在，所以应当筛选
allgene<-row.names(allsample_combined_cca@assays$RNA@meta.features)

OXPHOS_genes<-intersect(x=OXPHOS_genes, y = allgene)

glycolysis_gene<-intersect(x=glycolysis_gene, y = allgene)

triglyceride_metabolic_gene<-intersect(x=triglyceride_metabolic_gene, y = allgene)

X_gene<-subset(X_gene,X_gene %in% allgene)

Y_gene<-subset(Y_gene,Y_gene %in% allgene)

eXi_gene<-subset(eXi_gene,eXi_gene %in% allgene)

Xi_gene<-setdiff(X_gene, eXi_gene)

#计算比例
allsample_combined_cca[["percent.glycolysis_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= glycolysis_gene)

allsample_combined_cca[["percent.oxidative_phosphorylation_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= OXPHOS_genes)

allsample_combined_cca[["percent.triglyceride_metabolic_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)

allsample_combined_cca[["percent.X_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)

allsample_combined_cca[["percent.Y_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)

allsample_combined_cca[["percent.eXi_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)

allsample_combined_cca[["percent.Xi_gene"]] <- PercentageFeatureSet(object = allsample_combined_cca, features= triglyceride_metabolic_gene)

#各组平均表达热图
MARKERS <- AverageExpression(allsample_combined_cca, return.seurat = T)

DoHeatmap(MARKERS, features =, size = 3,slot = "data", disp.min = 0, disp.max = 5)+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkorange")

DoHeatmap(MARKERS, features = c("AMH","INHBB","INHA"), size = 6)+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkorange")

