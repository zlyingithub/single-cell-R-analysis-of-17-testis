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
DimPlot(object = allsample_combined_cca, reduction = "tsne", group.by = "age",label.size = 8,pt.size = 0.1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())
DimPlot(object = allsample_combined_cca, reduction = "umap", group.by = "agetype",label.size = 8,pt.size = 1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())
DimPlot(object = allsample_combined_cca, reduction = "tsne",  label = F,label.size = 8,pt.size = 0.1)+ theme(axis.title=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.line = element_blank())

FeaturePlot(allsample_combined_cca, features = c( "FGFR3","SYCP3", "TNP1", "SOX9", "DLK1","MYH11","VWF","NOTCH2", "CD163","TPSAB1","CD3D"), cols = c("blue", "yellow"),min.cutoff = "q9",pt.size = 0.1, reduction = "tsne") #SSC
FeaturePlot(allsample_combined_cca, features = c( "DDX4","FGFR3","SYCP3", "TNP1","VIM", "SOX9", "DLK1","MYH11","VWF","NOTCH3","PTPRC", "CD163","TPSAB1","CD3D"), cols = c("blue", "yellow"),min.cutoff = "q9",pt.size = 0.1, reduction = "tsne") #SSC
VlnPlot(allsample_combined_cca_10X_KSandOA_germcell,features = c("VWF"),pt.size = 0.1,group.by = "seurat_GERMclusters",split.by = "age")+scale_fill_manual(values =mypalette)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 30),plot.title =element_text(size = 45,face="plain") )
VlnPlot(allsample_combined_cca_10X_OA,features = c("TIMP1"),pt.size = 0,group.by = "seurat_11clusters",cols = mypalette)+geom_boxplot(fill = "black",color="gray90",size=1,width=0.05,outlier.size = 0)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )
#带趋势线
VlnPlot(allsample_combined_cca,features = c("percent.glycolysis_gene"),pt.size = 0.1)+geom_boxplot(fill = "black",color="gray90",size=1,width=0.05,outlier.size = 0)+NoLegend()+theme(axis.title=element_blank(),axis.text = element_text(size = 35),plot.title =element_text(size = 45,face="plain") )+stat_summary(aes(group=1),fun.y=mean, geom="smooth", shape=1,size=2,color="gray")+ scale_y_continuous(limits = c(0,2))

allsample_combined_cca<-subset(x = allsample_combined_cca, idents = c("0","1","2","3","4","5","6","7","8","9","11","12","13","14","15","16"))
allsample_combined_cca <- RenameIdents(object = allsample_combined_cca, `0` = "MIX_cells", `1` = "MIX_cells",`2` = "spermatid",
                                       `3` = "spermatogonia",`4` = "Sertoli_cells",`5` = "Sertoli_cells",
                                       `6` = "spermatocyte",`7` = "spermatocyte",`8` = "spermatocyte",`9` = "macrophages",`11` = "VSM_cells",
                                       `12` = "endotheliocyte",`13` = "spermatid",`14` = "mast_cells",`15` = "T_cells",`16` = "MIX_cells")
allsample_combined_cca$seurat_noa9clusters<-Idents(allsample_combined_cca)
allsample_combined_cca$age<-factor(allsample_combined_cca$age,levels = c("OA","iNOA","KS","AZFa_Del"),ordered = F)
allsample_combined_cca$seurat_11clusters<-factor(allsample_combined_cca$seurat_11clusters,levels = c("SPGs","SPCs","SPTs","SCs","LCs","PTMs","ECs","VSMs","Mo&Mφ","MCs","TCs"),ordered = F)
MARKERS <- FindAllMarkers(object = allsample_combined_cca_10X_KS, only.pos = F, min.pct = 0.1, logfc.threshold = 0.25)
MARKERS <- FindMarkers(object = allsample_combined_cca_10X_OA, only.pos = F,ident.1 = "SCs", min.pct = 0.1, logfc.threshold = 0.25)
write.csv(MARKERS,"KS 内部11群DEGs.csv")
#计算基因平均值
MARKERS <- AverageExpression(allsample_combined_cca, return.seurat = T)
#计算metadata里的平均值（表格，按照xx分组list属性，平均值或其他）
MARKERS<-aggregate(allsample_combined_cca_10X_KSandOA@meta.data,list(allsample_combined_cca_10X_KSandOA@meta.data$seurat_11clusters,allsample_combined_cca_10X_KSandOA@meta.data$age),mean)

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


#基因表达相关性回归分析
MARKER<-allsample_combined_cca_10X_KS_SC@meta.data
MARKERS<-allsample_combined_cca_10X_KS_SC@assays$RNA@data
MARKERS<-data.frame(MARKERS)
MARKERS<-t(MARKERS)
MARKER<-cbind(MARKER,MARKERS)

ggplot(data = MARKER, mapping = aes(x = XIST, y = percent.chr_neXi_gene))+ 
  geom_point(aes(color = percent.chr_X_gene,size= nCount_RNA))+ #以drv为分组设置点的颜色
  geom_smooth(method = 'lm', formula = y ~ x,colour="ORANGE",size=2)+labs(caption ="y = 4.092-0.698x  p-value = 8.185e-14")
MARKERS<-lm(percent.chr_eXi_gene~XIST,MARKER)
summary(MARKERS) 




#各组平均表达热图
MARKERS <- AverageExpression(allsample_combined_cca_10X_OA, return.seurat = T)
DoHeatmap(MARKERS, features = c("USP9Y","DDX3Y","UTY","HSFY1","HSFY2","RBMY1A1","RBMY1B","RBMY1C","RBMY1D","RBMY1E","RBMY1J","DAZ1","DAZ2","DAZ3","DAZ4","BPY2","CDY1B","PRY","CSPG4P1Y","DAZL","BOLL","SRY","ZFX","ZFY"), size = 3,slot = "data", disp.min = 0, disp.max = 5)+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkorange")
DoHeatmap(MARKERS, features = c("AMH","INHBB","INHA"), size = 6)+scale_fill_gradient2(low = "steelblue", mid = "white", high = "darkorange")

