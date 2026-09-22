
print(Sys.time())
rm(list = ls());gc()
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

wd <- file.path(root, '09.GSEA')
dir.create(wd, showWarnings = F, recursive = T)
setwd(wd)

####---- load packages ----####
library(lance)
library(readr)
library(msigdbr)
library(biomaRt)
library(tidyverse)
library(patchwork)
library(enrichplot)
library(data.table)
library(org.Hs.eg.db)
library(clusterProfiler)
library(dplyr)
library(corrplot)
library(GseaVis)
library(RColorBrewer)
library(psych)

# 如需指定本机 R 库路径，可在此处设置：.libPaths(c("你的R库路径", .libPaths()))


# KEGG 基因集：优先读取本地 assets/geneset 下的 gmt 文件（v7.5.1 版本）；
# 若本地没有，则用 msigdbr 获取 C2:CP:KEGG 基因集并缓存到该路径，供下次直接使用。
kegg_file <- file.path(root, "assets", "geneset", "c2.cp.kegg.v7.5.1.symbols.gmt")
if (file.exists(kegg_file)) {
  kegg <- read.gmt(kegg_file)
} else {
  message("未找到本地 KEGG 基因集文件：", kegg_file, "，改用 msigdbr 获取 C2:CP:KEGG 基因集。")
  kegg <- as.data.frame(msigdbr(species = "Homo sapiens", category = "C2", subcategory = "CP:KEGG")[, c("gs_name", "gene_symbol")])
  colnames(kegg) <- c("term", "gene")
  dir.create(dirname(kegg_file), showWarnings = FALSE, recursive = TRUE)
  write.table(kegg, kegg_file, sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
}
# kegg <- kegg[c(3,4,5)]

data <-  read.csv('../00.rawdata/GSE16561/01.data_GSE16561.csv',row.names = 1) 

gene <- read.csv('../06.expressionROC/keygenes.csv')[,1]


rt<-data

for (ii in gene){  
  tryCatch({
    set.seed(1)
    tar.exp <- rt[ii,]
    y <- as.numeric(tar.exp)
    
    rdata_file <- paste0(ii, '.cor.Rdata')
    if (file.exists(rdata_file)) {
      
      load(rdata_file)
    } else {
      
      data1 <- data.frame() 
      for (i in rownames(rt)) {
        dd <- corr.test(as.numeric(rt[i,]), y, method = "spearman", adjust = "fdr")
        data1 <- rbind(data1, data.frame(gene = i, cor = dd$r, p.value = dd$p))
      }
      data1 <- data1[order(data1$cor, decreasing = TRUE), ]
      gene1 <- data1$cor
      names(gene1) <- mapIds(org.Hs.eg.db, keys = data1$gene, column = 'ENTREZID',
                             keytype = 'SYMBOL', multiVals = 'filter')
      names(gene1) <- data1$gene
      gene1 <- na.omit(gene1)
      save(gene1, file = rdata_file)
    }
    
    kk=GSEA(gene1, TERM2GENE=kegg[,1:2], pvalueCutoff = 1)
    kkTab=as.data.frame(kk)
    kkTab=kkTab[kkTab$p.adjust<0.05 & abs(kkTab$NES)>1,]
    kkTab <- kkTab%>% arrange(pvalue)
    kkTab$Description <- gsub('KEGG_','',kkTab$Description )
    kkTab$Description <- gsub('_',' ',kkTab$Description )
    kkTab$Description <- stringr::str_to_title(kkTab$Description)##TitleCase
    write.csv(kkTab,paste0(ii,'_gsea.csv'),row.names = F)
    
    target_ids <- kkTab$ID[1:5]
    # p1=gseaNb(object = kk,geneSetID = target_ids,curveCol = (rainbow(10)),subPlot = 3)+
    #   labs(title = ii)
    p <- gseaplot2( kk,  kkTab$ID[1:5],
                    base_size =10,
                    color = c('#7B68EE', "#ef99a6", "#34B5B6", "#F0AD4E", "#38bdf8"),
                    rel_heights = c(1.5, 0.3, 0.5),
                    title = paste0("Gene Set Enrichment Analysis of ",ii))
    
    # revise legend position
    p[[1]] <- p[[1]] + theme(legend.position = "right",
                             legend.text = element_text(size = 12, colour = "black"),
                             axis.text.y = element_text(size = 12, colour = "black"),
                             axis.title.y = element_text(size = 14, face = "bold", colour = "black"),
                             plot.title = element_text(size = 18, face = "bold", color = "black", hjust = .5))
    p[[3]] <- p[[3]] + theme(axis.text.x = element_text(size = 12, colour = "black"),
                             axis.text.y = element_text(size = 12, colour = "black"),
                             axis.title.x = element_text(size = 14, face = "bold", colour = "black"),
                             axis.title.y = element_text(size = 14, face = "bold", colour = "black"))
    
    pdf(paste0(ii,'.GSEA.kegg.pdf'),width = 12,height = 6,family = 'Times')
    print(p)
    dev.off()
    png(paste0(ii,'.GSEA.kegg.png'),width =12,height =6,units = 'in',res = 600,family = 'Times')
    print(p)
    dev.off()
    
  },error=function(e){})
}



library(dplyr)
library(KEGGREST)
library(ggtext)
library(GSEABase)
library(aplot)  # 
library(ComplexHeatmap)
pdf.options(family="Times")


# com pathway------------------

gene <- read.csv('../06.expressionROC/keygenes.csv')[,1]

files <- paste0(gene,"_gsea.csv")
files

features <- gene
features


df_list <- lapply(files,function(file){
  read.csv(file)
})
df_list <- stats::setNames(df_list,features)
sapply(df_list, nrow)
# CYP1B1   DEGS1   HSDL2 OSBPL1A    PTEN 
# 27      25      30      78      76 
description_list <- lapply(df_list, function(df) df$Description)
description_list <- stats::setNames(description_list, features)

# Venn
library(VennDiagram)
library(venn)
library(ggvenn)


common_description <- Reduce(intersect, description_list)
length(common_description)  # 

write.csv(common_description, 'common_term.csv', row.names = FALSE)

venn_list <- description_list
library(ggVennDiagram)

conflicted::conflicts_prefer(ggplot2::alpha)

# ColorBrewer Dark2（5色）
my_colors <- c("#1B9E77", "#D95F02", "#7570B3", "#E7298A", "#66A61E")

p2 <- ggVennDiagram(venn_list, 
                    set_color = my_colors,
                    label = "count") +           
  scale_fill_gradient(low = "white",       # 连续填充：白到指定色
                      high = "white") +
  guides(fill = "none", color = "none")




# p <- ggvenn(venn_list,show_percentage = FALSE,digits = 0,stroke_size = 0.2,text_size = 6)#digits
# 
# p2 <- p + scale_fill_manual(values = c("#98CFE6", "#ADE7A8","#F39F4E","#EEB7D3")) 

pdf('Venn.pdf',width = 8,height = 8)
p2
dev.off()
png("Venn.png", width =8, height = 8, units = 'in',  res = 600) #   
p2
dev.off() 


common_description
common_description %>% length() 
if (length(common_description) > 20) {
  common_description <- common_description[1:20]
}

nes_list <- lapply(names(df_list), function(feature) {
  df <- df_list[[feature]]
  df_fiter <- df %>% dplyr::filter(Description %in% common_description) %>% dplyr::select(Description,NES) %>% dplyr::rename( !!feature := NES)
})
nes_list <- stats::setNames(nes_list,names(df_list))

nes_wide <- Reduce(function(x, y) merge(x, y, by = "Description", all = TRUE), nes_list) %>% tibble::column_to_rownames("Description")
write.csv(nes_wide,file = "NES_Heatmap.csv")

nes_list <- lapply(names(df_list), function(feature) {
  df <- df_list[[feature]]
  df_fiter <- df %>% dplyr::filter(Description %in% common_description) %>% dplyr::select(Description,p.adjust) %>% dplyr::rename( !!feature := p.adjust)
})
nes_list <- stats::setNames(nes_list,names(df_list))

nes_p <- Reduce(function(x, y) merge(x, y, by = "Description", all = TRUE), nes_list) %>% tibble::column_to_rownames("Description")
write.csv(nes_p,file = "NES_p_Heatmap.csv")

mat <- as.matrix(nes_wide)
pval <- as.matrix(nes_p)

lim <- max(abs(mat), na.rm = TRUE)
library(circlize)
library(viridisLite)


# # # 7级颜色：深蓝 → 蓝 → 浅蓝 → 白 → 浅红 → 红 → 深红
# col_fun <- colorRamp2(
#   breaks = c(-3, -1.5, -1, -0.5, 0, 0.5, 1, 1.5, 3),
#   colors = c("#053061", "#2166AC", "#4393C3", "#92C5DE",
#              "#F7F7F7",
#              "#FDBB84", "#F67E4B", "#D6604D", "#B2182B")
# )
# # # 11级 RdBu 配色（红-蓝 diverging）
# col_fun <- colorRamp2(
#   breaks = seq(-3, 3, length.out = 11),
#   colors = rev(brewer.pal(11, "RdBu"))  # rev() 让蓝色在负值，红色在正值
# )
mat[mat>3] <- 3
mat[mat<(-3)] <- (-3)
col_fun <- colorRamp2(
  breaks = c(-3, -1,  0,  1,  3),
  colors = c( "#4393C3", "#92C5DE",
              "#F7F7F7",
              "#FDBB84", "#F67E4B")
)

ht <- Heatmap(
  mat,
  name = "NES",
  col = col_fun,
  na_col = "#F2F2F2",
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  rect_gp = gpar(col = NA),
  row_names_gp = gpar(fontfamily = "Times", fontsize = 10),
  show_column_names = FALSE,  # 隐藏默认列名
  row_names_max_width = max_text_width(rownames(mat), gp = gpar(fontsize = 10)),
  heatmap_legend_param = list(
    title = "NES",
    at = c(-lim, 0, lim),
    labels = c(sprintf("%.2f", -lim), "0", sprintf("%.2f", lim))
  ),
  cell_fun = function(j, i, x, y, width, height, fill) {
    # 只在第一行绘制列名（居中）
    if (i == 1) {
      grid.text(
        colnames(mat)[j], 
        x = x, 
        y = unit(1, "npc") + unit(2, "mm"),  # 放在热图上方
        gp = gpar(fontfamily = "Times", fontsize = 10, fontface = "bold"),
        just = "center"  # 关键：文本居中
      )
    }
    
    # 原有的显著性标记
    pval_ij <- pval[i, j]
    stars <- ""
    if (!is.na(pval_ij)) {
      if (pval_ij < 0.001) stars <- "***"
      else if (pval_ij < 0.01) stars <- "**"
      else if (pval_ij < 0.05) stars <- "*"
    }
    if (stars != "") {
      grid.text(stars, x, y, gp = gpar(fontsize = 13, col = "black"))
    }
  }
)
draw(ht, padding = unit(c(1, 1, 10, 1), "mm"))  # 下、左、上、右边距


pdf("NES_Heatmap.pdf", width = 8, height = 4, family = "Times")
draw(ht, padding = unit(c(1, 1, 10, 1), "mm"))
dev.off()

png("NES_Heatmap.png", width = 8, height = 4,family = 'Times',units = "in",res = 600)
draw(ht, padding = unit(c(1, 1, 10, 1), "mm"))
dev.off()





