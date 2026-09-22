
rm(list = ls());gc()

suppressPackageStartupMessages({
  library(dplyr)
  library(GEOquery)
  library(tidyverse)
  library(data.table)
})
conflicted::conflict_prefer_matching("^filter$|^select$|^arrange$|^mutate$|", "dplyr",quiet = T)

####---- 项目根目录定位（脚本位于 <项目根>/Rscript/，脚本内全部使用相对路径）----####
.script_dir <- function() {
  a <- commandArgs(FALSE)
  i <- grep("^--file=", a)
  if (length(i) > 0) return(dirname(normalizePath(sub("^--file=", "", a[i[1]]), mustWork = FALSE)))
  if (requireNamespace("rstudioapi", quietly = TRUE) && rstudioapi::isAvailable()) {
    p <- tryCatch(rstudioapi::getSourceEditorContext()$path, error = function(e) "")
    if (nzchar(p)) return(dirname(normalizePath(p, mustWork = FALSE)))
  }
  for (k in seq_len(sys.nframe())) {
    f <- tryCatch(sys.frame(k)$ofile, error = function(e) NULL)
    if (!is.null(f) && nzchar(f)) return(dirname(normalizePath(f, mustWork = FALSE)))
  }
  if (basename(getwd()) == "Rscript") getwd() else file.path(getwd(), "Rscript")
}
SCRIPT_DIR <- .script_dir()
ROOT_CAND <- unique(c(SCRIPT_DIR, dirname(SCRIPT_DIR), getwd(), dirname(getwd())))
ROOT_HIT <- ROOT_CAND[file.exists(file.path(ROOT_CAND, "00.rawdata"))]
root <- if (nzchar(Sys.getenv("PROJECT_ROOT"))) {
  normalizePath(Sys.getenv("PROJECT_ROOT"), mustWork = FALSE)
} else if (length(ROOT_HIT) > 0) {
  normalizePath(ROOT_HIT[1], mustWork = FALSE)
} else if (basename(SCRIPT_DIR) == "Rscript") {
  normalizePath(dirname(SCRIPT_DIR), mustWork = FALSE)
} else {
  normalizePath(SCRIPT_DIR, mustWork = FALSE)
}
options(project.root = root)   # 供脚本中 rm(list = ls()) 之后继续使用
cat("项目根目录：", root, "\n")
print(root)
setwd(root)
wd <- file.path(root,"14.scRNA/")
dir.create(wd,recursive = T,showWarnings = F)
setwd(wd)
work.path <- wd
saveRDS(work.path, "workpath.rds")


# 单细胞数据集 GSE225948  ----------------------

dir <- "00.rawdata"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

# library(utils)
# untar("GSE225948_RAW.tar", exdir = "./GSE225948")

library(Seurat)
library(GEOquery)
library(data.table)
library(tidyverse)
library(lance)
library(Matrix)

## 00.load data--------

dir <- "./GSE225948/"
samples=list.files('./GSE225948','counts')
sceList = lapply(samples,function(pro){ 
  print(pro) 
  ct=fread(file.path( dir ,pro),data.table = F)
  ct[1:4,1:4]
  rownames(ct)=ct[,1]
  ct=ct[,-1]
  sce =CreateSeuratObject(counts =  ct ,
                          project =   pro%>% str_sub(1,10),
                          min.cells = 3,
                          min.features = 200 )
  return(sce)
})
scRNA=merge(x=sceList[[1]],
            y=sceList[ -1 ],
            add.cell.ids =samples%>% str_sub(1,10) )


Idents(scRNA) <- "orig.ident"
table(Idents(scRNA) )
# mitochondria proportion
scRNA[["percent.mt"]] <- PercentageFeatureSet(scRNA, pattern = "^mt-") ##mouse
range(scRNA[["percent.mt"]])


## add group
library(tidyverse)
smps = Cells(scRNA) %>% str_sub(1,10)


gset <- getGEO("GSE225948", destdir = "./", GSEMatrix =  T, getGPL = F)
phen.dat <- pData(gset[[1]])
colnames(phen.dat)
table(phen.dat$`treatment:ch1`)
table(phen.dat$`tissue:ch1`)

group <- data.frame(sample = phen.dat$geo_accession,
                    group= phen.dat$`treatment:ch1`
)
group <- group%>%filter(sample%in%unique(smps))
group$group <- gsub(" \\d+$","",group$group)
table(group$group)

group$group <- ifelse(group$group == "Sham","Normal","Disease")
write.csv(group, file = "./01.GSE225948_group.csv", quote = F)

grps = read.csv("./01.GSE225948_group.csv",row.names = 1)

table(smps)
group.label = grps$group[match(smps,grps$sample)]
names(group.label) = Cells(scRNA)
scRNA = AddMetaData(scRNA,group.label,col.name = 'group')
table(scRNA$group)
system.time(save(scRNA, file = "scRNA_orig.Rdata"))






####01 quality_control---------
rm(list = ls())
setwd(readRDS("../workpath.rds"))
dir <- "01.quality_control"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

load("../00.rawdata/scRNA_orig.Rdata")

# before quality control
paste0("质控前细胞数量为", length(colnames(scRNA)))
# "质控前细胞数量为35914"
paste0("质控前基因数量为", length(rownames(scRNA)))
# "质控前基因数量为16175"

# plot theme
theme.set = theme(
  axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  # axis.text.x = element_blank(),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 18,face='bold',family = "Times"),
  text = element_text(family = "Times"))

# plot features
plot.featrures = c("nFeature_RNA", "nCount_RNA", "percent.mt")

# before control violin plot 
library(patchwork)
plots <- list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA,
                       split.by="orig.ident",
                       pt.size = 0,
                       features = plot.featrures[i]) + theme.set+ NoLegend()}

violin.before <- wrap_plots(plots = plots, nrow = 1, ncol = 3, cols = palettes(category = "random", 22, show_col=F))
violin.before

# dev.off()
ggsave("01.before_quality_control.pdf", plot = violin.before, width = 12, height = 6)
ggsave("01.before_quality_control.png", plot = violin.before, width = 12, height = 6,dpi = 600)

# quality rate
minGene <- quantile(scRNA$nFeature_RNA,.02)   
maxGene <- quantile(scRNA$nFeature_RNA,.98)   
maxUMI <- quantile(scRNA$nCount_RNA,.95)      
minUMI <- quantile(scRNA$nCount_RNA,.02)  

# quality threshold
pctMT <- 15
minGene <- 200
maxGene <- 1000
maxUMI <- 2000
minUMI <- 200

# quality control
scRNA <- subset(scRNA, subset =
                  nCount_RNA > minUMI &
                  nCount_RNA < maxUMI & 
                  nFeature_RNA > minGene & 
                  nFeature_RNA < maxGene& 
                  percent.mt < pctMT)

# after quality control
paste0("质控后细胞数量为", length(colnames(scRNA)))
# "质控后细胞数量为34928"
paste0("质控后基因数量为", length(rownames(scRNA)))
# "质控后基因数量为16175"

# after quality control
plots = list()
for(i in seq_along(plot.featrures)){
  plots[[i]] = VlnPlot(scRNA,
                       split.by="orig.ident",
                       pt.size = 0,
                       features = plot.featrures[i]) + theme.set + NoLegend()}
violin <- wrap_plots(plots = plots, nrow = 1, ncol = 3, cols = palettes(category = "random", 22, show_col=F))
# violin
ggsave("02.after_quality_control.pdf", plot = violin, width = 12, height = 6)
ggsave("02.after_quality_control.png", plot = violin, width = 12, height = 6,dpi = 600)


# save system
system.time(save(scRNA, file = "scRNA.qc.rdata"))



####----02. seurat integrate data ----####
rm(list = ls())
setwd(readRDS("../workpath.rds"))
dir <- "02.intergration"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

# load packages
library(Seurat)
library(tidyverse)

# load seurat 
load("../01.quality_control/scRNA.qc.rdata")

combined <- NormalizeData(scRNA)
combined <- FindVariableFeatures(combined)
combined <- ScaleData(combined)
combined <- RunPCA(combined)

# # Dimensionality reduction for pca plot before integrated 
# pdf("01.before_intergrate_pca.pdf", width = 7, height = 5)
# DimPlot(combined, reduction = "pca",group.by = 'orig.ident')
# dev.off()
# png("01.before_intergrate_pca.png", width = 7, height = 5, units = 'in', res = 600)
# DimPlot(combined, reduction = "pca",group.by = 'orig.ident')
# dev.off()
DimPlot(combined, reduction = "pca",group.by = 'group')
DimPlot(combined, reduction = "pca",group.by = 'orig.ident')

# # save system 
# system.time(save(combined , file = "unintergrated.rdata"))
# load("unintergrated.rdata")
# 
library(harmony)
scRNA <- combined
scRNA <- RunHarmony(scRNA,'orig.ident',lambda = 0.5, theta = 1)#保守增大lambda减少theta
scRNA <- JoinLayers(scRNA)
DimPlot(scRNA, reduction = "harmony",group.by = 'orig.ident')
system.time(save(scRNA, file = "scRNA.integrated.rdata"))

# combined <- IntegrateLayers(
#   object = combined,
#   method = CCAIntegration,
#   orig.reduction = "pca",
#   new.reduction = "integrated.cca",
#   k.anchor = 20,
#   verbose = FALSE)
# scRNA <- combined
# scRNA <- JoinLayers(scRNA)
# system.time(save(scRNA, file = "scRNA.integrated.rdata"))


# find Hypervariable genes
scRNA <- FindVariableFeatures(scRNA, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(scRNA), 10)
cat(top10, sep = ", ")
# Camp, Ngp, Ltf, Pf4, Ccl5, Ppbp, Gzma, Saa3, Gng11, Chil3

# plot theme set
theme.set= theme(#axis.title.x=element_blank(),
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 10,  family = "Times"),
  axis.text.y = element_text(size = 14,  family = "Times"),
  legend.text = element_text(size = 14, family = "Times"),
  legend.title = element_text(size = 16,face='bold',family = "Times"),
  text = element_text(family = "Times"))

# plot variable features with labels and nonelabels
p.1 <- VariableFeaturePlot(scRNA) & theme.set
p.2 <- LabelPoints(plot = p.1, 
                   points = top10, 
                   repel = T)
p <- p.1 + p.2 
ggsave(filename = "01.feature_selection.pdf", p.2, width = 7, height = 6)
ggsave(filename = "01.feature_selection.png", p.2, width = 7, height = 6, dpi = 600)

# save system
scRNA.norm <- scRNA
system.time(save(scRNA.norm, file = "scRNA.norm.rdata"))




####----03 dimensionality reduction of pca ----####
rm(list = ls());gc()
setwd(readRDS("../workpath.rds"))
dir <- "03.dimensionality"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

# load seurat data
load("../02.intergration/scRNA.norm.rdata")

# dimensionality reduction of pca
## npc Usually choose 30 , use 50 to find optimal number of cluster
scRNA.norm.pca <- RunPCA(scRNA.norm, features = VariableFeatures(object = scRNA.norm), npcs = 50)

# pca genes
print(scRNA.norm.pca[["pca"]], dims = 1:5, nfeatures = 5)

# View contributions
VizDimLoadings(scRNA.norm.pca, dims = 1:2, reduction = "pca")

# The abscissa is the value of the score, and the greater the absolute value of the score, the greater the correlation
# The two groupings together represent no heterogeneity
DimPlot(scRNA.norm.pca, reduction = "pca")

# heatmap 
#DimHeatmap(scRNA.norm.pca, dims = 1:5, cells = 2000, balanced = TRUE)

# plot theme set
theme.set = theme(
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 18,face='bold',family = "Times"),
  text = element_text(family = "Times"))

# according to decision pc numbers
pdf(file="01.elbowplot.pdf",width=8,height=6)
ElbowPlot(scRNA.norm.pca,reduction="pca",ndims = 50) +theme.set 
dev.off()

png(file="01.elbowplot.png",width=8,height=6,units='in',res=600)
ElbowPlot(scRNA.norm.pca,reduction="pca",ndims = 50)+theme.set  
dev.off()


# find the optimal number of clusters
## NOTE: This process can take a long time for big datasets, comment out for expediency. 
## More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
## "num.replicate": number of replicate, "dims": numbers of pca
scRNA.norm.pca <- JackStraw(scRNA.norm.pca, num.replicate = 100, dims = 50)


## ScoreJackStraw()：It was used to quantify the significance intensity of principal components, 
## the principal components with more genes with low p.value were more statistically significant
scRNA.norm.pca <- ScoreJackStraw(scRNA.norm.pca, dims = 1:50)
system.time(save(scRNA.norm.pca, file = "scRNA.norm.pca.Jack.Rdata"))
load("scRNA.norm.pca.Jack.Rdata")

# Compare the p.value distribution and the $ uniform distribution for each principal component
plot.jackstraw <- JackStrawPlot(scRNA.norm.pca, dims = 1:50,ymax = 0.3) & theme.set
plot.jackstraw
ggsave(filename = '02.pca_jackstraw.pdf', plot.jackstraw, width = 13, height = 8)
ggsave(filename = '02.pca_jackstraw.png', plot.jackstraw, width = 13, height = 8,dpi = 600)


####--- dimensionality reduction of tsne ----###
###---- load packages ----###
library(scales)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
# if (!requireNamespace("IOBR", quietly = TRUE))
# devtools::install_github("IOBR/IOBR")
library(IOBR)

###---- tsne ----###
# use the flattening quantity according to the graph
pc.select <- 15

# calculate proximity distances
scRNA.norm <- FindNeighbors(object = scRNA.norm.pca, dims = 1:pc.select)

library(clustree)
library(patchwork)
seq <- seq(0.1, 1.2, by = 0.1)
scRNA.norm.pca.clu <- scRNA.norm 
for(res in seq){
  
  scRNA.norm.pca.clu <- FindClusters(scRNA.norm.pca.clu,
                                     resolution = res,
                                     graph.name = "RNA_snn",
                                     algorithm = 1,        # Louvain
                                     verbose = FALSE)
}
# system.time(save(scRNA.norm.pca.clu, file = "scRNA.norm.pca.clu.Rdata"))  
p1 <- clustree(scRNA.norm.pca.clu, prefix = 'RNA_snn_res.') + coord_flip()
ggsave("clustree.pdf",p1,width = 10,height = 10)

# group the cells
# resolution determines the number of cluster, ranging from 0.4 to 1.2
# the bigger the resolution, the more clusters you get
scRNA.norm <- FindClusters(object = scRNA.norm, resolution = 0.5)

# dimensionality reduction of tsne
scRNA.norm  <- RunTSNE(scRNA.norm , dims = 1:pc.select)
scRNA.norm <- RunUMAP(scRNA.norm , dims = 1:pc.select)
table(scRNA.norm$seurat_clusters)
# plot theme set
theme.set <- theme(
  axis.title = element_text(size = 20, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  legend.title = element_text(size = 20,face = 'bold',family = "Times"),
  text = element_text(family = "Times"))

# 配色：优先使用本地 assets/mycolor.R 中的 color 向量；缺失时使用内置备用配色
mycolor_file <- file.path(getOption("project.root", root), "assets", "mycolor.R")
if (file.exists(mycolor_file)) {
  source(mycolor_file)
  mycols <- color[51:70]
} else {
  warning("未找到 ", mycolor_file, "，已改用内置备用配色。", call. = FALSE)
  mycols <- rep(c(wes_palette("GrandBudapest2", 4), wes_palette("Zissou1", 3),
                  wes_palette("Darjeeling1", 5), wes_palette("Darjeeling2", 5),
                  wes_palette("GrandBudapest1", 4)), length.out = 20)
}

# all sample tsne 
pdf(file = "03.tsne_all.pdf", width = 10, height = 8, family = 'Times')
DimPlot(scRNA.norm ,
        reduction = "tsne",
        label.size = 6, label = TRUE, pt.size = 1.0, 
        cols = mycols) + theme.set
dev.off()

png(file = "03.tsne_all.png", width = 10, height = 8, family = 'Times', units = 'in', res = 600)
DimPlot(scRNA.norm ,
        reduction = "tsne",
        label.size = 6, label = TRUE, pt.size = 1.0,
        cols = mycols) + theme.set
dev.off()

# # group tsne
# pdf(file = "06.tsne_group.pdf", width = 14, height = 8, family = 'Times')
# DimPlot(scRNA.norm ,
#         reduction = "tsne",
#         label.size = 6, label = T, pt.size = 1.0,
#         split.by  = "group",
#         cols = palettes(category = "random", 24, show_col = F)) + theme.set
# dev.off()
# 
# png(file = "06.tsne_group.png", width=14, height=8, family = 'Times', units = 'in', res = 600)
# DimPlot(scRNA.norm ,
#         reduction = "tsne",
#         label.size = 6, label = T, pt.size = 1.0,
#         split.by  = "group",
#         cols = palettes(category = "random", 22, show_col = F)) + theme.set
# dev.off()

pdf(file = "03.umap_all.pdf", width = 10, height = 8, family = 'Times')
DimPlot(scRNA.norm ,
        reduction = "umap",
        label.size = 6, label = TRUE, pt.size = 1.0, 
        cols = mycols) + theme.set
dev.off()

png(file = "03.umap_all.png", width = 10, height = 8, family = 'Times', units = 'in', res = 600)
DimPlot(scRNA.norm ,
        reduction = "umap",
        label.size = 6, label = TRUE, pt.size = 1.0,
        cols = mycols) + theme.set
dev.off()

# save system
system.time(save(scRNA.norm, file = "scRNA.norm.rdata"))


####---- 04.cell annotation ----####
rm(list = ls());gc()
setwd(readRDS("../workpath.rds"))
dir <- "04.annotation"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

###---- load packages ----###
library(SingleR)
library(celldex)
# library(IOBR)
library(tidyverse)
library(colorspace)
library(Seurat)
library(reshape2)

###---- find markers ----###
# load seurat data
load("../03.dimensionality/scRNA.norm.rdata")
scRNA.norm <- JoinLayers(scRNA.norm)
# # find al markers
all.markers <- FindAllMarkers(scRNA.norm,
                              only.pos = TRUE, min.pct = 0.25,
                              logfc.threshold = 0.5, test.use = 'wilcox', return.thresh = 0.01)
# output markers
write.csv(all.markers, file = '01.all_markers.csv', row.names = T, quote = F)
# 
# # choose top 
top <- all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
# # plot theme set
# theme.set <- theme(
#   axis.title.x = element_text(size = 20, face = "bold", family = "Times"),
#   axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
#   legend.text = element_text(size = 16, face = "bold", family = "Times"),
#   legend.title = element_text(size = 18,face = 'bold',family = "Times"))
# 
# pdf(file = '02.top_heatmap.pdf', width = 8, height = 8)
# DoHeatmap(scRNA.norm, features = top$gene, label = F)+
#   labs(title = "", x = "Cells separated by clusters", y = "", size = 40) + theme.set
# dev.off()
# png(file = '02.top_heatmap.png', width = 8, height = 8, units = 'in', res = 600)
# DoHeatmap(scRNA.norm, features = top$gene, label = F)+
#   labs(title = "", x = "Cells separated by clusters", y = "", size= 40) & theme.set
# dev.off()
# 


###---- annotation by self ----###
DefaultAssay(scRNA.norm) <- "RNA"

##marker list #PMID: 38177281
list.gene <- list(
  `B cells` = c('Igkc', 'Ighd', 'Ms4a1', 'Ly6d'),
  Neutrophils = c('S100a8','S100a9','Retnlg','Cxcr2'),#,'Camp'
  `T cells` = c('Trbc1','Cd3d','Trac','Ms4a4b'),
  Monocytes = c('Ccr2', 'S100a4', 'Chil3', 'Ms4a6c','Plac8','Cybb'),
  `NK cells` = c('Gzma', 'Ccl5', 'Klra8', 'Nkg7'),
  Bas.Eos = c('Ccl3','Ms4a2','Gata2','Ccr3'),
  Dendritics  = c('Cd209a','Bst2','H2-DMb1')
  # Hem.pre = c( 'Mki67', 'Top2a', 'Pclaf')
)

# DotPlot(scRNA.norm, features = list.gene) + RotatedAxis() + theme.set+
#   scale_color_gradientn(colours = c('#0f59a4','#EFC971','#d11a2d'))

# plot theme set
theme.set <- theme(
  axis.title = element_text(size = 18, face = "bold", family = "Times"),
  axis.text.x = element_text(angle = 45, vjust = 1,size = 14, family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))


# experssion dot plot 
Idents(scRNA.norm) <- "seurat_clusters"

pdf(file = "02.dotplot_clusters.pdf", width = 10, height = 5, family = 'Times')
DotPlot(scRNA.norm, features = list.gene) + RotatedAxis() + theme.set+
  scale_color_gradientn(colours = c('#0f59a4','#EFC971','#d11a2d'))
dev.off()

png(file = "02.dotplot_clusters.png", width = 10, height = 5, family = 'Times', units = 'in', res = 600)
DotPlot(scRNA.norm, features = list.gene) + RotatedAxis() + theme.set+
  scale_color_gradientn(colours = c('#0f59a4','#EFC971','#d11a2d'))
dev.off()

# change cluster name 
cluster.ids = c(
  "0" = "B cells",  
  "1" = "Neutrophils",   
  "2" = "B cells" ,
  "3" = "Neutrophils", 
  "4" = "Neutrophils",   
  "5" = "T cells",  
  "6" = "Monocytes",
  "7" = "Monocytes", 
  "8" = "NK cells", 
  "9" = "Neutrophils",
  "10"= "Bas.Eos",
  "11"= "Dendritics",
  "12"= "Bas.Eos"
)

sce <- scRNA.norm
# cell name to cluster
scRNA.norm <- RenameIdents(scRNA.norm, cluster.ids)

scRNA.norm$celltype <- Idents(scRNA.norm)
table(scRNA.norm@meta.data$celltype)
# B cells Neutrophils     T cells   Monocytes    NK cells     Bas.Eos  Dendritics 
# 11894       13995        4111        3878         750         178         122 

# plot theme set
theme.set2 = theme(
  axis.title = element_text(size = 16, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 12, face = "bold", family = "Times"),
  axis.text.y = element_text(size = 12, face = "bold", family = "Times"),
  legend.text = element_text(size = 14, face = "bold", family = "Times"),
  legend.title = element_blank(),
  text = element_text(family = "Times"))

# mycols <- palettes(category = "random", 20,show_col = F)
library(wesanderson)
mycols <-  c(
  wes_palette("GrandBudapest2", 4),
  wes_palette("Zissou1", 3),
  wes_palette("Darjeeling1", 5),
  wes_palette("Darjeeling2", 5),
  wes_palette("GrandBudapest1", 4)
  
)
# annotation  plot
pdf(file = "04.celltype.pdf", width = 8, height = 6, family = 'Times')
DimPlot(scRNA.norm, reduction = "umap",
        label.size = 5,label.color = 'black',label = T,
        cols = mycols) + theme.set2
dev.off()

png(file = "04.celltype.png", width = 8, height = 6, family = 'Times', units = 'in', res = 600)
DimPlot(scRNA.norm, reduction = "umap",
        label.size = 5, label.color = 'black', label =T,
        cols = mycols) + theme.set2
dev.off()


# marker gene experssion dot plot
Idents(scRNA.norm) <- "celltype"
pdf(file = '03.marker.exp.pdf', width = 10,height = 5, family = 'Times')
DotPlot(scRNA.norm, features = list.gene)+
  RotatedAxis() + theme.set+
  scale_color_gradientn(colours = c('#0f59a4', 'white', '#d11a2d'))
dev.off()
png(file = '03.marker.exp.png', width = 10, height = 5, units = 'in', family = 'Times', res = 600)
DotPlot(scRNA.norm, features = list.gene)+
  RotatedAxis() + theme.set+
  scale_color_gradientn(colours = c('#0f59a4','white','#d11a2d'))
dev.off()

# save system
scRNA.celltype <- scRNA.norm
system.time(save(scRNA.celltype, file = "scRNA.celltype.Rdata"))

# marker gene experssion
library(dplyr)
library(tidyr)

df_marker <- list.gene %>%
  tibble::enframe(name = "cell_type", value = "gene") %>% 
  unnest(gene) %>% 
  group_by(cell_type) %>%
  summarise(gene = paste(gene, collapse = ","), .groups = "drop") %>%
  left_join(
    cluster.ids %>%
      tibble::enframe(name = "cluster", value = "cell_type") %>%
      mutate(first_clus = as.integer(cluster)) %>%   # 提取第一个编号
      group_by(cell_type) %>%
      summarise(cluster = paste(cluster, collapse = ","),
                first_clus = min(first_clus), .groups = "drop")
  ) %>%
  arrange(first_clus) %>%   # 按第一个 cluster 排序
  select(cell_type, marker=gene, cluster)
df_marker
write.csv(df_marker,file = '05.final markers.csv', row.names = F)




####---- 05.key cells ----####
rm(list = ls());gc()
setwd(readRDS("../workpath.rds"))
dir <- "05.key_cells"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

####load packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(ggsignif)
library(tidyverse)
library(reshape2)
library(tidyr)
library(tibble)

load("../04.annotation/scRNA.celltype.Rdata")
scRNA.norm <- scRNA.celltype
cell_type <- data.frame(sample=scRNA.norm$orig.ident,Celltype=scRNA.norm$celltype,group=scRNA.norm$group)
cell_type$sample <- substr(rownames(cell_type),1,10)
colnames(cell_type)
cell_type$Celltype <- as.character(cell_type$Celltype)

celltype <- as.data.frame(table(cell_type$Celltype))
celltype$Var1 <- as.character(celltype$Var1)


# 定义函数来确保包含两个组
ensure_two_groups <- function(data, groups = c("Disease", "Normal")) {
  # 检查数据中是否包含所有组
  missing_groups <- setdiff(groups, unique(data$group))
  # 如果有缺失的组，添加虚拟的零计数行
  if (length(missing_groups) > 0) {
    for (group in missing_groups) {
      data <- rbind(data, data.frame(Celltype = unique(data$Celltype), group = group, value = 0))
    }
  }
  
  return(data)
}


library(rstatix)
result <- list()
for (i in c(1:nrow(celltype))) {
  mydata <- cell_type[which(cell_type$Celltype==celltype$Var1[i]),]
  ka.data <- xtabs(~mydata$Celltype+mydata$group,data = mydata)
  # 检查列联表是否包含两个组
  if (ncol(ka.data) == 1) {
    # 如果只有一个组，添加一个虚拟的零计数列
    ka.data <- cbind(ka.data, 0)
    colnames(ka.data) <- c("Disease", "Normal")
  }
  
  stat_res <- chisq_test(ka.data) %>%
    adjust_pvalue(method = "BH") %>%
    add_significance("p.adj")
  result[[i]] <- stat_res
}

names(result) <- celltype$Var1
result <- do.call(rbind,lapply(result,data.frame))
result$celltype <- celltype$Var1
write.csv(result,file = '01.res.chiseq.csv')

sample_num<- length(unique(cell_type$sample))
plot.celltype <- table(cell_type$Celltype, cell_type$sample) %>% as.data.frame() %>% recast(Var1 ~ Var2)
colnames(plot.celltype)[1] = "Cell"
plot.celltype <- cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))
percent.celltype<- plot.celltype[-c(2:(sample_num+1))]
percent.celltype[2:ncol(percent.celltype)]<- percent.celltype[2:ncol(percent.celltype)]*100

dat.celltype <- reshape2::melt(percent.celltype)
colnames(dat.celltype) <- c('Cell','Sample','Percentage')
dat.celltype$Group <- cell_type$group[match(dat.celltype$Sample,cell_type$sample)]
dat.celltype$Group<- factor(dat.celltype$Group,levels = c("Normal","Disease"))


barplot <- ggplot(dat.celltype,aes(x =Cell, y = Percentage, fill = Group)) +
  geom_boxplot(width=0.5,
               alpha=0.8,
               position = position_dodge(0.9),
               outlier.shape = NA)+ 
  annotate(geom = "text", x = result$celltype, y = 100, size = 3, family = "Times",
           label =as.character(result$p.signif)) +
  # stat_compare_means(aes(group = Group),
  #                    method = "wilcox.test",
  #                    symnum.args = list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "ns")),
  #                    label = "p.signif")+
  theme_classic()+
  theme(legend.position = "top")+
  theme(axis.title.x =element_text(size=15,family = "Times", face = "bold"),
        axis.text.x =element_text(angle=45,size=10,hjust = 1,family = "Times", face = "bold"
        ),
        axis.title.y =element_text(size=15,family = "Times", face = "bold"),
        axis.text.y=element_text(size=15,family = "Times", face = "bold"))+
  theme(legend.title=element_text(size=15, family = "Times", face = "bold") ,
        legend.text=element_text(size=14, family = "Times", face = "bold"))+
  ggsci::scale_fill_d3()


ggsave(barplot,filename = '01.DEcell.pdf',w=6,h=4)
ggsave(barplot,filename = '01.DEcell.png',w=6,h=4)




library(celldex)
library(reshape2)
library(tidyverse)
library(colorspace)
library(Seurat)

load('../04.annotation/scRNA.celltype.Rdata')
scRNA.norm <- scRNA.celltype

# cell rate graph data
all.count <- data.frame(celltype = scRNA.norm$celltype, sample=scRNA.norm$group)
table(all.count$celltype)

plot.celltype <- table(all.count$celltype, all.count$sample) %>% as.data.frame() %>% recast(Var1 ~ Var2)
colnames(plot.celltype)[1] = "Cell"
plot.celltype <- cbind(plot.celltype, apply(plot.celltype[-1], 2, proportions))

colnames(plot.celltype)[2:5] = c( "Count_Disease","Count_Normal", "Proportion_Disease", "Proportion_Normal")###
plot.celltype <- plot.celltype[order(plot.celltype$Count_Normal, decreasing = T), ]
write.csv(plot.celltype, "02.cellcountrate.csv", row.names = F, quote = F)

dat.celltype <- plot.celltype[c(1,4:5)]
dat.celltype[2:3] <- 100 * dat.celltype[2:3]

dat.plot <- reshape2::melt(dat.celltype, id.vars = "Cell", variable.name = "Sample", value.name = "value")
dat.plot$Cell = factor(dat.plot$Cell, levels = plot.celltype$Cell)
dat.plot$label = round(dat.plot$value, 2) %>% paste0(.,"%")

library(wesanderson)
mycols <-  c(
  wes_palette("GrandBudapest2", 4),
  wes_palette("Zissou1", 3),
  wes_palette("Darjeeling1", 5),
  wes_palette("Darjeeling2", 5),
  wes_palette("GrandBudapest1", 4)
  
)
# cell rate graph 
ggplot(data = dat.plot, mapping = aes(x = value, y = Cell)) + 
  geom_text(aes(x = value + 1, label = label), hjust = 0) +
  geom_bar(aes(fill = Cell), stat = "identity") + 
  facet_wrap(~Sample) +
  scale_fill_manual(values = mycols) + 
  xlim(c(0,60)) +
  xlab("Cell Fraction") + ylab("Cell Type") +
  theme_minimal(base_size = 16) + 
  theme(panel.border = element_rect(fill = "transparent"),
        legend.position = "none", 
        panel.grid = element_blank(),
        text = element_text(family = "Times"),
        # legend.text = element_text(size = 18, face = "bold", family = "Times"),
        axis.text.x = element_text(size = 12, colour = 'black',family = "Times"),
        axis.text.y = element_text(size = 12, colour = 'black',family = "Times"),
        axis.title = element_text(size = 16, face = "bold", family = "Times"))

ggsave("02.cellproportion.pdf", width = 9, height = 5)
ggsave("02.cellproportion.png", width = 9, height = 5, units = "in", dpi = 600, bg = "white")


hubgene <- read.csv("../../06.expressionROC/keygenes.csv")[,1] %>% stringr::str_to_title()
setdiff(hubgene,rownames(scRNA.norm))
intersect(rownames(scRNA.norm), c("CP1B", "CYPIB1", "P4501b1"))
hubgene <- hubgene[-1]

table(scRNA.norm$group)
list_gene <- list(Key_genes = hubgene)

theme.set = theme(
  axis.title = element_text(size = 18, face = "bold", family = "Times"),
  axis.text.x = element_text(angle = 45, vjust = 1,size = 14, family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 18, face = "bold", family = "Times"),
  legend.title = element_blank(),
  legend.position = "left",
  text = element_text(family = "Times"))

pdf(file = '03.keygene.dotexp.pdf',width=6,height=6,family='Times')
DotPlot(scRNA.norm,features = list_gene)+
  RotatedAxis()+theme.set+
  scale_color_gradientn(colours = c('darkblue','orange','#FF2000'))
dev.off()
png(file = '03.keygene.dotexp.png',width=6,height=6, units = 'in',family='Times', res = 600)
DotPlot(scRNA.norm,features = list_gene)+
  RotatedAxis()+theme.set+
  scale_color_gradientn(colours = c('darkblue','orange','#FF2000'))
dev.off()

Idents(scRNA.norm) <- scRNA.norm$celltype

hub <- as.data.frame(hubgene)
colnames(hub) <- 'Gene'

p1 = FeaturePlot(scRNA.norm,reduction = "umap",
                 features = hub$Gene,repel = T, label = T,label.size = 5,
                 slot = "data", cols = c("lightgrey","red"), ncol = 3, order = T)

pdf(file = '04.keygene_cellmap.pdf',w=14,h=6,family='Times')
p1
dev.off()
png(file = '04.keygene_cellmap.png',w=14,h=6,units = "in",res=600,family='Times')
p1
dev.off()
scRNA.norm <- JoinLayers(scRNA.norm, assay = "RNA")

#####diff wilcox
for (i in hub$Gene) {
  expression_matrix <- data.frame(expr=GetAssayData(scRNA.norm, slot = "data")[i,])
  dat.test <- data.frame(id=scRNA.norm$orig.ident,Cell=scRNA.norm$celltype,expr=expression_matrix$expr,group=scRNA.norm$group)
  table(dat.test$Cell,dat.test$group)
  dat.test_summary <- dat.test %>%
    group_by(group, Cell) %>%
    summarise(expr_sum = sum(expr))
  filtercell <- dat.test_summary$Cell[dat.test_summary$expr_sum>0]%>%as.character()
  dat.test <- dat.test%>%filter(Cell%in%filtercell)#
  library(rstatix)
  stat_res <- dat.test %>%
    group_by(Cell) %>%
    wilcox_test(expr ~ group) %>%
    adjust_pvalue(method = "BH") %>%  # method BH == fdr
    add_significance("p")
  dat.test2 <- dat.test
  dif_cell2 <- stat_res
  library(ggpubr)
  p <- ggviolin(dat.test2 , x = "Cell", y = "expr",
                fill= "group",
                palette = c( "#6495ed","#dc143c"),
                #outlier.shape = 0,
                outlier.size = 0.8,
                bxp.errorbar = T,
                #add = "mean_se"
  ) +
    ylim(0,max(dat.test$expr)+max(dat.test$expr)/5)+
    stat_pvalue_manual(dif_cell2 ,
                       x = "Cell",
                       y.position = max(dat.test$expr)+max(dat.test$expr)/10,
                       size =3.5,
                       color = "black",
                       family = "Times",
                       # label = paste0("~italic(p)=={p.adj2}"),
                       label = "p.signif"
                       #parse = T
    ) +
    theme_bw()+
    labs(x = "", y = "expr", color = "") +
    labs(title = paste0(i))+
    theme(axis.title.x = element_text(size = 18, face = "bold", family = "Times"),
          axis.title.y = element_text(size = 18, face = "bold", family = "Times"),
          axis.text.x = element_text(size = 13, color = 'black',face = "bold",
                                     family = "Times", angle = 45,
                                     vjust = 1, hjust = 1),
          axis.text.y = element_text(size = 13, color = 'black',face = "bold",
                                     family = "Times"),
          legend.text = element_text(size = 15, family = "Times",face = "bold"),
          legend.title = element_text(size = 18, family = "Times",hjust = 0.5,face = "bold"),
          #plot.margin = ggplot2::margin(t=.3,b=0,l=2,r=.5, unit = "cm"),
          text = element_text(family = "Times"),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())+
    theme(legend.position = "right")
  p
  ggsave(filename = paste0("05.",i,"_violinplot.pdf"), width = 10, height = 7)
  ggsave(filename = paste0("05.",i,"_violinplot.png"), width = 10, height = 7,units = "in", dpi = 600)
}


###AUCell
rm(list = ls())

library(Seurat) 
library(tidyverse)
library(AUCell)
pdf.options(family="Times")


load("../04.annotation/scRNA.celltype.Rdata")

sce2 <- scRNA.celltype
DefaultAssay(sce2) <- "RNA"

hubgene <- read.csv("../../06.expressionROC/keygenes.csv")[,1] %>% stringr::str_to_title()
setdiff(hubgene,rownames(scRNA.norm))
intersect(rownames(scRNA.norm), c("CP1B", "CYPIB1", "P4501b1"))
hubgene <- hubgene[-1]

WNT_features <- list(hubgene)
names(WNT_features) <- "keygene"

exprMatrix <- GetAssayData(sce2, layer = "data")
rownames(exprMatrix) = Features(sce2)
colnames(exprMatrix) = Cells(sce2)


# cells_rankings <- AUCell_buildRankings(exprMatrix,splitByBlocks=TRUE) 
# cells_AUC <- AUCell_calcAUC(WNT_features, cells_rankings, 
#                             aucMaxRank=nrow(cells_rankings)*0.1)#只使用排名中前10%的基因

#也可以直接使用 AUCell_run 函数得到上述2个步骤相同的结果。
cells_AUC2 <- AUCell_run(exprMatrix, WNT_features)
saveRDS(cells_AUC2,"AUCcells.rds")

cells_assignment <- AUCell_exploreThresholds(cells_AUC2, plotHist=TRUE, assign=TRUE) 
auc_thr = cells_assignment$keygene$aucThr$selected 
auc_thr
sce <- sce2
Idents(sce) <- sce$celltype
sce$auc_score = as.numeric(getAUC(cells_AUC2))
sce$auc_group = ifelse(sce$auc_score>auc_thr,"high","low") 

dat<- data.frame(sce@meta.data, 
                 sce@reductions$umap@cell.embeddings,
                 seurat_annotation = sce@active.ident)
class_avg <- dat %>%
  group_by(seurat_annotation) %>% #按照seurat_annotation列(即细胞的分类)对数据进行分组。
  summarise(
    umap_1 = median(umap_1),
    umap_2 = median(umap_2) #对每个分组计算UMAP坐标的中位数 画label
  )

library(ggpubr)
library(ggrepel)
pdf.options(family="Times")
ggplot(dat, aes(umap_1, umap_2))  +
  geom_point(aes(colour  = auc_score)) +
  viridis::scale_color_viridis(option="D") +
  ggrepel::geom_label_repel(aes(label = seurat_annotation),
                            data = class_avg,
                            label.size = 0,
                            segment.color = NA)+
  theme_bw()

ggplot(dat, aes(umap_1, umap_2)) +
  geom_point(aes(colour = auc_score)) +
  scale_color_gradient2(low = "grey90", mid = "orange", high = "darkred", midpoint = 0.5) +
  ggrepel::geom_label_repel(aes(label = seurat_annotation),
                            data = class_avg,
                            label.size = 0,
                            segment.color = NA) +
  theme_bw()

ggsave(filename="03.Aucell_celltype.pdf",width = 9,height = 7)
ggsave(filename="03.Aucell_celltype.png",width = 9,height = 7,units = 'in',dpi=600)

# 高低分组
# ggplot(dat,aes(x = umap_1,y = umap_2))+
#   geom_point(aes(color = auc_group),size = 0.5)+
#   theme_classic()
ggplot(dat, aes(x = umap_1, y = umap_2)) +
  geom_point(aes(color = auc_group), size = 0.5, alpha = 0.6) +
  geom_text_repel(data = class_avg,
                  aes(x = umap_1, y = umap_2, label = seurat_annotation),
                  size = 4.5, fontface = "bold",
                  box.padding = 0.5,  # 标签间距
                  point.padding = 0.3,
                  max.overlaps = Inf) +
  scale_color_manual(values = c("low" ="lightgrey",  "high" = "#7F7FFF")) +
  theme_classic()

ggsave(filename="04.Aucell_aucgroup.pdf",width = 9,height = 7)
ggsave(filename="04.Aucell_aucgroup.png",width = 9,height = 7,units = 'in',dpi=600)

head(sce@meta.data)

ggviolin(sce@meta.data, x="celltype", y="auc_score", width = 0.6, 
         color = "black",#轮廓颜色
         fill="celltype",#填充
         palette = "npg",
         add = 'mean_sd',
         xlab = F, #不显示x轴的标签
         ylab = "AUCell Score", #不显示x轴的标签
         bxp.errorbar=T,#显示误差条
         bxp.errorbar.width=0.5, #误差条大小
         size=1, #箱型图边线的粗细
         outlier.shape=NA, #不显示outlier
         legend = "right")
ggsave(filename="05.Aucell_violin.pdf",width = 10,height = 5)
ggsave(filename="05.Aucell_violin.png",width = 10,height = 5,units = 'in',dpi=600)

# 按 celltype 分组，计算每个细胞类型在 disease vs Normal 间的差异
library(dplyr)
# 获取所有细胞类型
cell_types <- unique(sce$celltype)

# 存储结果
results <- data.frame()

for(ct in cell_types) {
  # 提取该细胞类型的数据
  ct_data <- subset(sce@meta.data, celltype == ct)
  
  # 检查是否有两组数据
  if(length(unique(ct_data$group)) < 2) next
  
  # 提取两组的 auc_score
  disease_auc <- ct_data$auc_score[ct_data$group == "Disease"]
  normal_auc <- ct_data$auc_score[ct_data$group == "Normal"]
  
  # Wilcoxon 秩和检验（非参数，更适合 auc_score 这类有偏分布）
  test_result <- wilcox.test(disease_auc, normal_auc)
  
  # 计算中位数和差异
  disease_median <- median(disease_auc)
  normal_median <- median(normal_auc)
  
  # 保存结果
  results <- rbind(results, data.frame(
    celltype = ct,
    disease_median = disease_median,
    normal_median = normal_median,
    median_diff = disease_median - normal_median,
    p_value = test_result$p.value,
    statistic = test_result$statistic
  ))
}

# 校正 p 值（FDR）
results$padj <- p.adjust(results$p_value, method = "BH")

# 按 p 值排序
results <- results[order(results$p_value), ]
print(results)

library(ggplot2)
library(ggpubr)

# 准备数据
plot_data <- as.data.frame(sce@meta.data) %>%
  dplyr::select(celltype, group, auc_score) %>%
  filter(!is.na(auc_score), group %in% c("Disease", "Normal")) %>%
  mutate(
    group = factor(group, levels = c("Normal", "Disease")),
    celltype = factor(celltype, levels = results$celltype[order(results$median_diff)])
  )


results <- results %>%
  arrange(median_diff) %>%
  mutate(celltype = factor(celltype, levels = celltype))

# 再设置 plot_data 的因子水平
plot_data <- as.data.frame(sce@meta.data) %>%
  dplyr::select(celltype, group, auc_score) %>%
  filter(!is.na(auc_score), group %in% c("Disease", "Normal")) %>%
  mutate(
    group = factor(group, levels = c("Normal", "Disease")),
    celltype = factor(celltype, levels = results$celltype)  # 使用 results 的排序
  )

#绘制（annotate 的 x 也要对应排序后的位置） ==
p_ggpubr <- ggplot(plot_data, aes(x = celltype, y = auc_score, fill = group)) +
  geom_violin(trim = FALSE, scale = "width", alpha = 0.9, position = position_dodge(width = .92)) +
  # geom_boxplot(width = 0.15, position = position_dodge(width = 0.7), 
  #              outlier.size = 0.5, alpha = 0.8) +
  scale_fill_manual(values = c("Normal" = "#38bdf8", "Disease" = "#f35e5a")) +
  labs(
    title = "AUC Score Distribution by Cell Type",
    x = NULL,
    y = "AUC Score",
    fill = "Group"
  ) +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12, face = "bold"),
    legend.position = "top",
    panel.grid.major.x = element_blank()
  ) +
  # annotate x 用 seq_along 对应排序后的 results
  annotate("text", 
           x = seq_along(results$celltype), 
           y = max(plot_data$auc_score, na.rm = TRUE) * 1.12,
           label = sapply(results$p_value, function(p) {
             if (p < 1e-100) return("p<1e-100\n****")
             else if (p < 0.0001) return(sprintf("p=%.1e\n****", p))
             else if (p < 0.001) return(sprintf("p=%.1e\n***", p))
             else if (p < 0.01) return(sprintf("p=%.1e\n**", p))
             else if (p < 0.05) return(sprintf("p=%.1e\n*", p))
             else return(sprintf("p=%.1e\nns", p))
           }),
           size = 3, vjust = 0.5)

print(p_ggpubr)

ggsave("06.plot_diff.pdf", p_ggpubr, width = 12, height = 6)
ggsave("06.plot_diff.png", p_ggpubr, width = 12, height = 6)




####---- 06.cell chat ----####
rm(list = ls())
setwd(readRDS("../workpath.rds"))
dir <- "06.cell_chat"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)
###---- load packages ----###
library(CellChat)
library(Seurat)
library(ggplot2)
library(tidyverse)

# load seurat data
load("../04.annotation/scRNA.celltype.Rdata")

# different cell seurat
scRNA.diff.cell <- scRNA.celltype
# scRNA.diff.cell <- JoinLayers(scRNA.celltype, assay = "RNA")
# normalized data matrix
data.input <- scRNA.diff.cell@assays$RNA$data 
meta <- scRNA.diff.cell@meta.data # a dataframe with rownames containing cell mata data
table(meta$group)

###---- Normal cell ----###
con.cell.use <- rownames(meta)[meta$group == "Normal"]
con.dat <- data.input[, con.cell.use]
con.meta <- meta[con.cell.use, ]

table(con.meta$celltype)
unique(con.meta$celltype)

# CellChat need matrix data and meta 
con.cellchat <- createCellChat(object = con.dat, meta = con.meta, group.by = "celltype")

# add meta data
con.cellchat <- addMeta(con.cellchat, meta = con.meta)

# change labels as your needed
con.cellchat <- setIdent(con.cellchat, ident.use = "celltype") 
levels(con.cellchat@idents)
unique(con.cellchat@idents)

con.cellchat@idents <- droplevels(con.cellchat@idents,
                                  exclude = setdiff(levels(con.cellchat@idents),
                                                    unique(con.cellchat@idents)))

# cell numbers
group.size <- as.numeric(table(con.cellchat@idents)) 

# human or mouse
# CellChatDB <- CellChatDB.human
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.interaction <- CellChatDB$interaction

# use “Secreted Signaling” anylsis 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
con.cellchat@DB <- CellChatDB.use

# The expression data is further preprocessed to save computing power
con.cellchat <- subsetData(con.cellchat)

# overexpression genes (ligand-receptor) are identified first
con.cellchat <- identifyOverExpressedGenes(con.cellchat)

# the overexpressed interactions between the overexpressed ligands and receptors are then identified
con.cellchat <- identifyOverExpressedInteractions(con.cellchat)

# # project gene expression data onto PPI network (optional) 
# # tip: human or mouse 
# # cellchat <- projectData(cellchat, PPI.human)
# con.cellchat <- projectData(con.cellchat, PPI.human)

# Calculate the probability of intercellular communication and predict the communication network
con.cellchat <- computeCommunProb(con.cellchat)

# output net 
con.net <- subsetCommunication(con.cellchat) 
write.csv(con.net,file = '01.net.Normal.csv', row.names = F, quote = F)

# the level of signaling pathways further infers intercellular communication and calculates aggregation networks
con.cellchat <- computeCommunProbPathway(con.cellchat)
con.cellchat <- aggregateNet(con.cellchat)
group.size <- as.numeric(table(con.cellchat@idents))

# # cell chat heatmap
# netVisual_heatmap(con.cellchat)
# pdf('01.number_of_interactions_Normal.pdf', width = 6, height = 5, family = "Times")
# netVisual_heatmap(con.cellchat)
# dev.off()
# png('01.number_of_interactions_Normal.png', width = 500, height = 400, family = "Times")
# netVisual_heatmap(con.cellchat)
# dev.off()

# save rds 
saveRDS(con.cellchat, file = 'Normalcellchat.rds')

# cell chat net of count 
pdf(file = '02.net_number_Normal.pdf', width = 7, height = 9, family = "Times")
netVisual_circle(con.cellchat@net$count,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()

png(filename = '02.net_number_Normal.png', width = 7, height = 9, units = 'in', res = 600, family = "Times")
netVisual_circle(con.cellchat@net$count,
                 vertex.weight = group.size, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
dev.off()

# cell chat net of weight
pdf(file = '03.net_weight_Normal.pdf', width = 7, height = 9, family = "Times")
netVisual_circle(con.cellchat@net$weight,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weight/strength")
dev.off()

png(filename = '03.net_weight_Normal.png',w=7,h=9,units='in',res=600, family="Times")
netVisual_circle(con.cellchat@net$weight,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weight/strength")
dev.off()

levels(con.cellchat@idents)
# show all the significant interactions (L-R pairs)
# dot plot
con.cellchat@data.signaling

pdf(file = '04.buble_Normal.pdf', width = 7, height = 5)
netVisual_bubble(con.cellchat, remove.isolate = FALSE, angle.x = 45)
dev.off()
png('04.buble_Normal.png', width = 7, height = 5, units = 'in', res = 600)
netVisual_bubble(con.cellchat, remove.isolate = FALSE,angle.x = 45#,sources.use = 1, targets.use = c(1:9)
                 )
dev.off()


###---- Disease cell ----###
table(meta$group)
di.cell.use <- rownames(meta)[meta$group == "Disease"]
di.dat <- data.input[, di.cell.use]
di.meta <- meta[di.cell.use, ]

table(di.meta$celltype)
unique(di.meta$celltype)

# CellChat need matrix data and meta 
di.cellchat <- createCellChat(object = di.dat, meta = di.meta, group.by = "celltype")

# add meta data
di.cellchat <- addMeta(di.cellchat, meta = di.meta)

# change labels as your needed
di.cellchat <- setIdent(di.cellchat, ident.use = "celltype") 
levels(di.cellchat@idents)
unique(di.cellchat@idents)

di.cellchat@idents <- droplevels(di.cellchat@idents,
                                 exclude = setdiff(levels(di.cellchat@idents),
                                                   unique(di.cellchat@idents)))

# cell numbers
group.size <- as.numeric(table(di.cellchat@idents)) 

# human or mouse
# CellChatDB <- CellChatDB.human
CellChatDB <- CellChatDB.mouse
showDatabaseCategory(CellChatDB)

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.interaction <- CellChatDB$interaction

# use “Secreted Signaling” anylsis 
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling")
di.cellchat@DB <- CellChatDB.use

# The expression data is further preprocessed to save computing power
di.cellchat <- subsetData(di.cellchat)

# overexpression genes (ligand-receptor) are identified first
di.cellchat <- identifyOverExpressedGenes(di.cellchat)

# the overexpressed interactions between the overexpressed ligands and receptors are then identified
di.cellchat <- identifyOverExpressedInteractions(di.cellchat)

# # project gene expression data onto PPI network (optional) 
# # tip: human or mouse 
# # cellchat <- projectData(cellchat, PPI.human)
# di.cellchat <- projectData(di.cellchat, PPI.human)

# Calculate the probability of intercellular communication and predict the communication network
di.cellchat <- computeCommunProb(di.cellchat)

# output net 
di.net <- subsetCommunication(di.cellchat) 
write.csv(di.net,file = '06.net.Disease.csv', row.names = F, quote = F)

# the level of signaling pathways further infers intercellular communication and calculates aggregation networks
di.cellchat <- computeCommunProbPathway(di.cellchat)
di.cellchat <- aggregateNet(di.cellchat)
group.size <- as.numeric(table(di.cellchat@idents))

# # cell chat heatmap
# netVisual_heatmap(di.cellchat)
# pdf('05.number_of_interactions_Disease.pdf', width = 6, height = 5, family = "Times")
# netVisual_heatmap(di.cellchat)
# dev.off()
# png('05.number_of_interactions_Disease.png', width = 500, height = 400, family = "Times")
# netVisual_heatmap(di.cellchat)
# dev.off()

# save rds 
saveRDS(di.cellchat, file = 'Diseasecellchat.rds')

# cell chat net of count 
pdf(file = '06.net_number_Disease.pdf', width = 7, height = 9, family = "Times")
netVisual_circle(di.cellchat@net$count,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Number of interactions")
dev.off()

png(filename = '06.net_number_Disease.png', width = 7, height = 9, units = 'in', res = 600, family = "Times")
netVisual_circle(di.cellchat@net$count,
                 vertex.weight = group.size, 
                 weight.scale = T, label.edge= F, 
                 title.name = "Number of interactions")
dev.off()

# cell chat net of weight
pdf(file = '07.net_weight_Disease.pdf', width = 7, height = 9, family = "Times")
netVisual_circle(di.cellchat@net$weight,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weight/strength")
dev.off()

png(filename = '07.net_weight_Disease.png',w=7,h=9,units='in',res=600, family="Times")
netVisual_circle(di.cellchat@net$weight,
                 vertex.weight = group.size,
                 weight.scale = T, label.edge= F,
                 title.name = "Interaction weight/strength")
dev.off()

levels(di.cellchat@idents)
# show all the significant interactions (L-R pairs)
# dot plot
di.cellchat@data.signaling
pdf(file = '08.buble_Disease.pdf', width = 7, height = 5)
netVisual_bubble(di.cellchat, remove.isolate = FALSE,angle.x = 45#,sources.use = 1, targets.use = c(1:9)
                 )
dev.off()
png('08.buble_Disease.png', width = 7, height = 5, units = 'in', res = 600)
netVisual_bubble(di.cellchat, remove.isolate = FALSE,angle.x = 45)
dev.off()


####---- 07.scMetabolism ----####
rm(list = ls());gc()
setwd(readRDS("../workpath.rds"))
dir <- "07.scMetabolism"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)



library(scMetabolism)
library(ggplot2)
library(rsvd)
library(dplyr)
library(Seurat)


load("../04.annotation/scRNA.celltype.Rdata")
table(scRNA.celltype$celltype)
scRNA <- scRNA.celltype

library(biomaRt)
library(Seurat)
mouse_genes <- rownames(scRNA)
human <- useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = "https://dec2021.archive.ensembl.org/") 
mouse <- useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://dec2021.archive.ensembl.org/")

MtoH <- getLDS(attributes = "mgi_symbol", # 要转换符号的属性，这里基因名（第3步是基因名）
               filters = "mgi_symbol", #参数过滤
               mart = mouse, #需要转换的基因名的种属来源，也就是第2步的mouse
               values = mouse_genes, #要转换的基因集
               attributesL = "hgnc_symbol", #要同源转换的目标属性，这里还是转为基因名，也可加其他
               martL = human, #要同源转换的目标种属，也就是第2步的human
               uniqueRows = TRUE)
head(MtoH)
MtoH_unique <- MtoH[!duplicated(MtoH$MGI.symbol), ]
write.csv(MtoH,"MtoH.csv")


MtoH <- read.csv("MtoH.csv",row.names = 1)
# 建立映射
MtoH_unique <- MtoH[!duplicated(MtoH$MGI.symbol), ]
gene_map <- setNames(MtoH_unique$HGNC.symbol, MtoH_unique$MGI.symbol)

# 获取可匹配的小鼠基因
genes_to_convert <- intersect(rownames(scRNA), names(gene_map))

# 提取子集
scRNA <- subset(scRNA, features = genes_to_convert)

# 获取新行名
new_rownames <- gene_map[rownames(scRNA)]
# 检查哪些 new_rownames 重复
dup_genes <- new_rownames[duplicated(new_rownames)]

# 方法：如果有重复，只保留第一个
if (length(dup_genes) > 0) {
  # 标记不重复的位置
  keep <- !duplicated(new_rownames)
  
  # 重新 subset scRNA，只保留不重复的
  scRNA <- subset(scRNA, features = rownames(scRNA)[keep])
  
  # 重新生成 new_rownames
  new_rownames <- gene_map[rownames(scRNA)]
  
  message("Removed ", sum(!keep), " duplicated genes, kept ", sum(keep))
}
rownames(scRNA) <- new_rownames
# rownames(scRNA) <- toupper(rownames(scRNA))

scRNA[['RNA']] <- as(scRNA[['RNA']],'Assay')
scRNA <- SeuratObject::UpdateSeuratObject(scRNA)

# 可选的本地基因集文件（仅用于构建 geneSetsList，后续代谢分析由 scMetabolism 自行读取基因集）
cellreports_file <- file.path(getOption("project.root", root), "assets", "CellReports.txt")
if (file.exists(cellreports_file)) {
  geneSets <- readLines(cellreports_file)
  geneSetsSplit <- strsplit(geneSets, split = "\t")
  geneSetsList <- list()
  for (i in seq_along(geneSetsSplit)) {
    setName <- geneSetsSplit[[i]][1]
    genes <- geneSetsSplit[[i]][3:length(geneSetsSplit[[i]])]
    geneSetsList[[setName]] <- genes
  }
  geneSets <- geneSetsList
} else {
  message("未找到可选基因集文件 ", cellreports_file, "，跳过 geneSetsList 构建（不影响后续分析结果）。")
  geneSetsList <- list()
}

# # trace('sc.metabolism.Seurat',edit = T,where = asNamespace('scMetabolism'))
if (file.exists("AUCell_results.rds")) {
  
  res <- readRDS("AUCell_results.rds")
} else {
  res <- sc.metabolism.Seurat(obj = scRNA, method = "AUCell",  imputation = F, ncores = 2, metabolism.type = "KEGG")
  
  saveRDS(res,"AUCell_results.rds")
}


score <- res@assays$METABOLISM$score
rownames(score)
score$sum <- rowSums(score)
score_sorted <- score[order(-abs(score$sum)), , drop = FALSE]%>%dplyr::select(sum)

input.pathways <- rownames(score_sorted)[1:10]

kycell <- colnames(scRNA)[scRNA$celltype=="Neutrophils"]
kycell <- gsub("-",".",kycell)
setdiff(kycell,names(score))
sub.score <- score[,kycell]
sub.score$sum <- rowSums(sub.score)
score_sorted <- sub.score[order(-abs(sub.score$sum)), , drop = FALSE]%>%dplyr::select(sum)
input.pathways <- rownames(score_sorted)[1:10]

input.pathways <- c("Citrate cycle (TCA cycle)", "Oxidative phosphorylation", "Fatty acid degradation",
                    "Valine, leucine and isoleucine degradation","Metabolism of xenobiotics by cytochrome P450",
                    "Drug metabolism - cytochrome P450",
                    "Glycolysis / Gluconeogenesis",  "Pyruvate metabolism",
                    "Fructose and mannose metabolism", "Galactose metabolism",
                    "Pentose phosphate pathway" #"Butanoate metabolism",
)

pdf("01.group_metaDotPlot.pdf", width = 7, height = 5, family = "Times")
print(DotPlot.metabolism(obj = res, pathway = input.pathways, phenotype = "group", norm = "y"))
dev.off()
png("01.group_metaDotPlot.png", width = 7, height = 5, res = 300,units = "in", family = "Times")
print(DotPlot.metabolism(obj = res, pathway = input.pathways, phenotype = "group", norm = "y"))
dev.off()

pdf("02.celltype_metaDotPlot.pdf", width = 9, height = 5, family = "Times")
print(DotPlot.metabolism(obj = res, pathway = input.pathways, phenotype = "celltype", norm = "y"))
dev.off()

png("02.celltype_metaDotPlot.png", width = 9, height = 5,units = "in", res = 300, family = "Times")
print(DotPlot.metabolism(obj = res, pathway = input.pathways, phenotype = "celltype", norm = "y"))
dev.off()


write.csv(score_sorted,"results.csv")
write.csv(input.pathways,"pathways.csv")







####---- 08.pseudotime anlysis ----####
rm(list = ls())
setwd(readRDS("../workpath.rds"))
dir <- "08.pseudotime"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)
###---- load packages ----###
library(tidyverse)
library(scales)
library(RColorBrewer)
library(ggsci)
library(wesanderson)
library(IOBR)
library(Seurat)
library(monocle)
library(harmony)

# load seurat data 
load("../04.annotation/scRNA.celltype.Rdata")# 

table(scRNA.celltype$celltype)
key.cell <- c('Neutrophils')

i=1

# hubgene
hubgene <-  read.csv("../../06.expressionROC/keygenes.csv")[,1] %>% stringr::str_to_title()
hubgene <- hubgene[-1]
scRNA.single.cell <- scRNA.celltype[ ,scRNA.celltype$celltype %in% c(key.cell[i])]

# dimensionality reduction
scRNA.single.cell <- NormalizeData(scRNA.single.cell,
                                   normalization.method = "LogNormalize", scale.factor = 10000)
scRNA.single.cell <- FindVariableFeatures(scRNA.single.cell,
                                          selection.method = 'vst', nfeatures = 2000)
scRNA.single.cell <- ScaleData(scRNA.single.cell)
scRNA.single.cell <- RunPCA(scRNA.single.cell, features = VariableFeatures(object = scRNA.single.cell))
scRNA.single.cell <- FindNeighbors(scRNA.single.cell, dims = 1:20)
scRNA.single.cell <- FindClusters(scRNA.single.cell, resolution = 0.2)
scRNA.single.cell <- RunUMAP(scRNA.single.cell, dims = 1:20)
table(scRNA.single.cell$seurat_clusters)

theme.set <- theme(
  plot.title = element_text(size = 20,face = "bold", family = "Times"),
  axis.title = element_text(size = 18, face = "bold", family = "Times"),
  axis.text.x = element_text(size = 14,  face = "bold", family = "Times"),
  axis.text.y = element_text(size = 14,  face = "bold", family = "Times"),
  legend.text = element_text(size = 16, face = "bold", family = "Times"),
  text = element_text(family = "Times"))
p <- DimPlot(scRNA.single.cell ,
             reduction = "umap",
             group.by = "seurat_clusters",
             label.size = 6, label = TRUE, pt.size=1.0,
             cols = palettes(category = "random", 14, show_col = F)) + theme.set+
  ggtitle(paste0('Subtype of ', key.cell[i]))
ggsave(filename = paste0('01.subtype_', key.cell[i], '.pdf'), p, width = 6, height = 5)
ggsave(filename = paste0('01.subtype_', key.cell[i], '.png'), p, width = 6, height = 5, dpi = 600, bg = "white")
system.time(save(scRNA.single.cell, file = paste0("scRNA.", key.cell[i], ".Rdata")))
# load(paste0("scRNA.", key.cell[i], ".Rdata"))



# pseudotime anlysis
single.matrix <- GetAssayData(scRNA.single.cell, layer = "count", assay = "RNA")
feature.ann <- data.frame(gene_id = rownames(single.matrix), gene_short_name = rownames(single.matrix))

# gene feature
rownames(feature.ann) <- rownames(single.matrix)

# AnnotatedDataFrame
single.fdat <- new("AnnotatedDataFrame", data = feature.ann)

sample.ann <- scRNA.single.cell@meta.data
cell.type <- Idents(scRNA.single.cell) %>% data.frame(.)
sample.ann <- cbind(sample.ann, cell.type)

#add subcluster cell information
colnames(sample.ann)[ncol(sample.ann)] = "Subtype"

single.pdat <- new("AnnotatedDataFrame", data = sample.ann)

single.cds <- newCellDataSet(single.matrix, phenoData = single.pdat, featureData = single.fdat,
                             expressionFamily = negbinomial.size())
single.cds <- estimateSizeFactors(single.cds)
single.cds <- estimateDispersions(single.cds)

# Free resources that aren't needed 
WGCNA::collectGarbage()

# add experssion
express.genes <- VariableFeatures(scRNA.single.cell)
single.exp.cds <- setOrderingFilter(single.cds, express.genes)

single.exp.cds <- reduceDimension(single.exp.cds, max_components = 2,
                                  verbose = T, norm_method = "log")   #, residualModelFormulaStr = "~nFeature_RNA+group"
save(single.exp.cds,file="single.exp.cds.rda")

# trace('project2MST',edit = T,where = asNamespace('monocle'))
single.exp.cds <- orderCells(single.exp.cds)

save(single.exp.cds, file = paste0(key.cell[i],"_cds_order.Rdata"))
load(paste0(key.cell[i],"_cds_order.Rdata"))

p1 <- plot_cell_trajectory(single.exp.cds, color_by = "Pseudotime", cell_size = 1, theta = 180,
                           size = 1, show_backbone=TRUE, show_branch_points = F) +
  theme(text = element_text(size = 10))
ggsave(paste0("02.", key.cell[i],"_Trajectory.png"), p1, width = 8, height = 5, units = "in", dpi = 600,path="./singleplot")
ggsave(paste0("02.", key.cell[i],"_Trajectory.pdf"), p1, width = 8, height = 5, units = "in",path="./singleplot")


p2 <- plot_cell_trajectory(single.exp.cds, color_by ="seurat_clusters", cell_size = 1, theta = 180,
                           size = 1, show_backbone=TRUE, show_branch_points = F) +
  scale_color_simpsons()+
  theme(text = element_text(size = 14))
ggsave(paste0("03.", key.cell[i], "_type_Trajectory.png"), p2, width = 12, height = 6, units = "in", dpi = 600,path="./singleplot")
ggsave(paste0("03.", key.cell[i], "_type_Trajectory.pdf"), p2, width = 12, height = 6, units = "in",path="./singleplot")

p3 <- plot_cell_trajectory(single.exp.cds, color_by = "State", cell_size = 1, theta = 180,
                           size = 1, show_backbone = TRUE, show_branch_points = F) +
  scale_color_npg() +
  theme(text = element_text(size = 14))
ggsave(paste0("04.", key.cell[i], "_State_Trajectory.png"), p3, width = 6, height = 6, units = "in", dpi = 600,path="./singleplot")
ggsave(paste0("04.", key.cell[i], "_State_Trajectory.pdf"), p3, width = 6, height = 6, units = "in",path="./singleplot")

p4 <- plot_cell_trajectory(single.exp.cds, color_by = "group", cell_size = 1, theta = 180,
                           size = 1, show_backbone = TRUE, show_branch_points = F) +
  scale_color_simpsons()+
  theme(text = element_text(size = 14))#+
  # facet_wrap(~ group, nrow = 1)
ggsave(filename = '05.celltype_Trajectory(group).pdf', p4, width = 10, height = 6, units = "in", dpi = 600,path="./singleplot")
ggsave(filename = '05.celltype_Trajectory(group).png', p4, width = 10, height = 6,path="./singleplot")

plot <- (p1+p2)/(p3+p4)
ggsave(filename = paste0('02.celltype_Trajectory_', key.cell[i], '.pdf'), plot, width = 12, height = 8, units = "in", dpi = 600)
ggsave(filename = paste0('02.celltype_Trajectory_', key.cell[i], '.png'), plot, width = 12, height = 8)

# cell experssion for pseudotime plot
cds.subset <- single.exp.cds[hubgene,]
theme.set <-   theme(legend.position = "right")+
  theme(axis.title.x = element_text(size = 22, color = 'black', face = "bold", family = 'Times'),
        axis.text.x = element_text(size = 18, family = 'Times'),
        axis.title.y = element_text(size = 22,color='black', face = "bold", family = 'Times'),
        axis.text.y = element_text(size = 18,  family = 'Times'),
        legend.title=  element_text(size = 20, color = 'black', face = "bold", family = 'Times'),
        legend.text = element_text(size = 18,  face = "bold", family = 'Times'),
        title = element_text(size = 20, color = 'black', face = "bold", family = 'Times'))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())



j1 <- monocle::plot_genes_in_pseudotime(cds.subset , color_by = "Pseudotime") + theme.set

j2 <- monocle::plot_genes_in_pseudotime(cds.subset , color_by = "State") + theme.set

j3 <-monocle::plot_genes_in_pseudotime(cds.subset , color_by = "group") + theme.set

# j4 <-monocle::plot_genes_in_pseudotime(cds.subset , color_by = "seurat_clusters") + theme.set
plot <- j1|j2|j3

ggsave(paste0('03.trajectory_gene_', key.cell[i], '.pdf'), plot, height = 8, width = 18, family='Times')
ggsave(paste0('03.trajectory_gene_', key.cell[i], '.png'), plot, height = 8, width = 18, dpi = 600)


# trajectory_heatmap
p <- plot_pseudotime_heatmap(
  cds.subset,cluster_rows = F,
  num_clusters = 4,
  cores = 1,
  show_rownames = TRUE ,return_heatmap = T
)

pdf(paste0('03.gene_pseud_heatmap_', key.cell[i], '.pdf'), width = 4.5, height = 3.5, onefile = FALSE, family = 'Times')
print(p)
dev.off()

png(paste0('03.gene_pseud_heatmap_', key.cell[i], '.png'), width = 4.5, height = 3.5, units = 'in',res = 600, family = 'Times')
print(p)
dev.off()

### 拟时序monocle3 ----
rm(list = ls())
setwd(readRDS("../workpath.rds"))
dir <- "08.pseudotime"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)

library(Seurat)
library(harmony)
library(tidyverse)
library(patchwork)
library(ggplot2)
library(viridis)
library(dplyr) 
library(tidyr)
library(tinyarray)
library(tidyverse)
library(GEOquery)
library(Biobase)
library(AnnoProbe) 
library(clusterProfiler)
library(readr)


# 重启R加载环境包
library(monocle3)

rm(list = ls());gc()

key.cell <- "Neutrophils"
i=1

load("scRNA.Neutrophils.Rdata")
UMAP <- scRNA.single.cell

table(UMAP@meta.data$seurat_clusters)
Idents(UMAP) <- "seurat_clusters"
unique(Idents(UMAP))

data <- GetAssayData(UMAP, assay = 'RNA', slot = 'counts')
cell_metadata <- UMAP@meta.data
gene_annotation <- data.frame(gene_short_name = rownames(data))
rownames(gene_annotation) <- rownames(data)
cds <- new_cell_data_set(data,
                         cell_metadata = cell_metadata,
                         gene_metadata = gene_annotation)
cds <- preprocess_cds(cds)
cds <- reduce_dimension(cds,preprocess_method = "PCA")
cds <- cluster_cells(cds)
cds <- learn_graph(cds,use_partition = F)


table(UMAP@meta.data$seurat_clusters)
Idents(UMAP) <- "seurat_clusters"
unique(Idents(UMAP))

cds.embed <- cds@int_colData$reducedDims$UMAP
int.embed <- Embeddings(UMAP, reduction = "umap")
int.embed <- int.embed[rownames(cds.embed),]
cds@int_colData$reducedDims$UMAP <- int.embed
plot_cells(cds, reduction_method="UMAP", color_cells_by="seurat_clusters") + ggtitle('int.umap')
table(cds@colData$seurat_clusters)


root_cells <- colnames(cds)[cds@colData$seurat_clusters == "0"]
cds <- order_cells(cds, root_cells  = root_cells)

p1 <- plot_cells(cds, color_cells_by = "pseudotime", label_roots=F,
                 label_leaves=F,
                 label_branch_points=F,
                 cell_size=0.8, 
                 min_expr = 5,
                 trajectory_graph_color='grey40',
                 group_label_size=3,
                 rasterize=F)+scale_color_viridis_c(option = "plasma",
                                                    name = "Pseudotime")
p1
p2 <- plot_cells(cds, color_cells_by = "group", label_roots=F,
                 label_leaves=F,
                 label_branch_points=F,
                 label_cell_groups=F,
                 cell_size=0.8,
                 trajectory_graph_color='grey40',
                 group_label_size=3,
                 rasterize=F)+ 
  # scale_color_manual(values = c('Control' ='#74a893' ,'DFU' ="red" ))
  scale_color_brewer(palette = "Set2", name = "group")
p2

p3 <- plot_cells(cds, color_cells_by = "seurat_clusters", label_roots=F,
                 label_leaves=F,
                 label_branch_points=F,
                 label_cell_groups=F,
                 cell_size=0.8,
                 trajectory_graph_color='grey40',
                 group_label_size= 3,
                 rasterize=F)+
  guides(color = guide_legend(title = "subtype"))
p3


p_all<- p1|p2|p3
p_all
ggsave(filename = paste0("04.monocle3_pseudotime_",key.cell[i],".pdf"),plot = p_all,device = cairo_pdf,width = 14,height = 5)
ggsave(filename = paste0("04.monocle3_pseudotime_",key.cell[i],".png"),plot = p_all,width = 14,height = 5,units = 'in',dpi = 600)


##hub基因表达--
library("ggiraph")
library("ggraph")
library("patchwork")
colors <-c("#dc8e97","#e3d1db","#74a893","#ac9141","#5ac6e9",
                    "#ebce8e","#f9766e","#4e79a6","#7587b1","#c7deef",
                    "#e1a4c6","#916ba6","#cb8f82","#7db3af","#d2e0ac",
                    "#d5231d","#3777ac","#4ea64a","#8e4c99","#e88f18",
                    "#e47faf","#b698c5","#a05528","#58a6d6","#1f2d6f",
                    "#279772","#add387","#d9b71a","#fbbab6","#e97371",
                    "#e1c548","#f0e2a3","#aedd2f","#d7ee96","#a199be",
                    "#5fa664","#abd0a7","#ca6a6b","#e5b5b5","#e5c06e",
                    "#bac4d0","#45337f")
                    
gene <- read.csv("../../06.expressionROC/keygenes.csv")[,1] %>% stringr::str_to_title()
gene <- gene[-1]

p6 <- plot_genes_in_pseudotime(cds[c(gene), ],
                               color_cells_by = "pseudotime",
                               ncol = 1) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"),
        legend.text = element_text(size = 8))+
  scale_color_viridis(option = "D")+guides()
p6

p7 <- plot_genes_in_pseudotime(cds[c(gene), ],
                               color_cells_by = "seurat_clusters",
                               ncol = 1) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))+
  scale_color_manual(values = colors)+guides(
    color = guide_legend(
      nrow = 4,
      byrow=T,
      reverse = F))
p7


p8 <- plot_genes_in_pseudotime(cds[c(gene), ],
                               color_cells_by = "group",
                               ncol = 1) +
  theme(axis.title = element_text(face = "bold"),
        strip.text = element_text(face = "bold"))+
  scale_color_manual(values = colors)+guides(
    color = guide_legend(
      nrow = 2,
      byrow=T,
      reverse = F))
p8

all <- (p6+p7+p8) +
  plot_layout(nrow = 1) &
  theme(
    legend.position = "top",
    legend.direction = "horizontal")
all

ggsave(filename = paste0("05.monocle3_sdheat_hubgene",key.cell[i],".pdf"),plot = all,device = cairo_pdf,width = 10,height = 7)
ggsave(filename = paste0("05.monocle3_sdheat_hubgene",key.cell[i],".png"),plot = all,width = 10,height = 7,units = 'in',dpi = 600)

save.image("monocle3.Rdata")


####---- 09.ReactomeGSA ----####
rm(list = ls());gc()

setwd(readRDS("../workpath.rds"))
dir <- "09.ReactomeGSA"
if (! dir.exists(dir)){dir.create(dir)}
setwd(dir)


###---- load packages ----###
library(ReactomeGSA)
library(Seurat)

# load seurat data
load("../04.annotation/scRNA.celltype.Rdata")
DefaultAssay(scRNA.celltype) <- "RNA"
table(Idents(scRNA.celltype))

scRNA.diff.cell <- scRNA.celltype
# gsva anlysis
gsva.result <- analyse_sc_clusters(scRNA.diff.cell, verbose = TRUE)
saveRDS(gsva.result,"gsva.result.rds")

# pathway
pathway.expression <- ReactomeGSA::pathways(gsva.result)
colnames(pathway.expression) <- gsub("\\.Seurat", "", colnames(pathway.expression))

write.csv(pathway.expression, file = "01.gsva_results.csv", quote = F, row.names = T)

# find the maximum differently expressed pathway
max.difference <- do.call(rbind, apply(pathway.expression, 1, function(row) {
  values <- as.numeric(row[2:length(row)])
  return(data.frame(name = row[1], min = min(values), max = max(values)))
}))

max.difference$diff <- max.difference$max - max.difference$min
max.difference <- max.difference[order(max.difference$diff, decreasing = T), ]

write.csv(max.difference, "02.max_difference.csv", quote = F, row.names = T)
max.difference$name[1:20]

# gsva heatmap
plot.num <- 20
plot.gsva <- pathway.expression[rownames(max.difference[1:plot.num,]),]

macro.heat <- pheatmap::pheatmap(plot.gsva[,-1],      # Use a data matrix and do not transpose
                                 scale = "column",    # Normalize the columns
                                 angle_col = 45,      #
                                 cellwidth = 25,      #
                                 cellheight = 15,     #
                                 labels_row = plot.gsva[,1],  # labels
                                 color = colorRampPalette(c("#4f87ff", "white", "#F06A21"))(50)  #
)
dev.off()

ggsave(filename = "03.gsva_heatmap.pdf", macro.heat, width = 10, height = 6)
ggsave(filename = "03.gsva_heatmap.png", macro.heat, width = 10, height = 6, dpi = 600, bg = "white")




