

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
setwd(root)
dir <- "01.different_genes"
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)



###---- load packages ----#####
library(limma); library(magrittr); library(stringr); library(tidyverse); library(readr);
library(psych); library(ggplot2); library(ggthemes); library(RColorBrewer); library(Ipaper);
library(scales); library(ComplexHeatmap); library(ggrepel)

###---- data set ----####
# 数据集与分组标签
gse <- "GSE16561"
disease <- "Disease"
control <- "Control"
###---- load data ----####
data <- read.csv(file = "../00.rawdata/GSE16561/01.data_GSE16561.csv", row.names = 1)
range(data)

group <- read.csv(file = "../00.rawdata/GSE16561/02.group_GSE16561.csv", row.names = 1)

# factor group, design.mat
table(group$group)
group$group <- factor(group$group, levels = c('Control', 'Disease'))
design.mat <- cbind(Control = ifelse(group$group == "Control", 1, 0), 
                    Disease = ifelse(group$group == "Control", 0, 1))

# limma anlysis
contrast.mat <- makeContrasts(contrasts = "Disease-Control", levels = design.mat)
fit <- lmFit(data, design.mat)
fit <- contrasts.fit(fit, contrast.mat)
fit <- eBayes(fit)
fit <- topTable(fit, coef = 1, number = Inf, adjust.method = "fdr")
DEG <- na.omit(fit)

# fc and pvalue cutoff
logfc.cutoff <- 0.5
pvalue.cutoff <- 0.05

# filter the data
DEG$change <- as.factor(ifelse(DEG$adj.P.Val < pvalue.cutoff & abs(DEG$logFC) > logfc.cutoff,
                               ifelse(DEG$logFC > logfc.cutoff ,'UP','DOWN'), 'NOT')
)
sig.diff <- subset(DEG, DEG$adj.P.Val < pvalue.cutoff & abs(DEG$logFC) > logfc.cutoff)

nrow(sig.diff)
paste0("最终获得 ",nrow(sig.diff)," 个差异表达基因")
# "最终获得 601 个差异表达基因"
table(sig.diff$change)
# DOWN  NOT   UP 
#  138    0  463 

# ouput data
write.csv(DEG, file = paste0("01.DEGs_all(", gse, ").csv"))
write.csv(sig.diff, file = paste0("02.DEGs_sig(", gse, ").csv"))



###---- volcano plot of different gene ----####
# extract top10 for up and down
dat.rep <- DEG[rownames(DEG) %in%
                 rownames(rbind(head(sig.diff[order(sig.diff$logFC,
                                                    decreasing = T),], 10),
                                head(sig.diff[order(sig.diff$logFC,
                                                    decreasing = F),], 10))), ]
# plot data 
DEG$Symbols <- rownames(DEG)
up_DEG <- DEG[which(DEG$change == "UP"),]
up_DEG <- up_DEG[order(-up_DEG$logFC, up_DEG$P.Value),]
down_DEG <- DEG[which(DEG$change == "DOWN"),]
down_DEG <- down_DEG[order(down_DEG$logFC, down_DEG$P.Value),]
DEG$change <- factor(DEG$change, levels =  c("UP", "NOT", "DOWN"))
max_lfc <- max(DEG$logFC)
min_lfc <- min(DEG$logFC)

# output volcano plot 
if(T){
  volcano_plot <- ggplot(data = DEG, 
                         aes(x = logFC,
                             y = -log10(adj.P.Val), 
                             color = change)) +
    scale_color_manual(values = c("#ed556a", "darkgray","#15559a")) +
    geom_point(size = 2, alpha = 0.6, na.rm=T) +
    theme_bw(base_size = 12) +
    geom_vline(xintercept = c(-logfc.cutoff, logfc.cutoff),
               lty = 4,col = "darkgray", lwd = 0.6)+
    geom_hline(yintercept = -log10(0.05),
               lty = 4,col = "darkgray", lwd = 0.6)+
    theme(legend.position = "right",
          panel.grid = element_blank(),
          legend.title = element_blank(),
          legend.text = element_text(face="bold",color="black",size=13),
          plot.title = element_text(hjust = 0.5, face = "bold", color = "black", size = 18),
          axis.text.x = element_text(face = "bold",color = "black",size = 15),
          axis.text.y = element_text(face = "bold",color = "black",size = 15),
          axis.title.x = element_text(face = "bold",color = "black",size = 15),
          axis.title.y = element_text(face = "bold",color = "black",size = 15),
          plot.subtitle = element_text(hjust = 0.5, size = 12, face = "italic", colour = "black")) +
    geom_label_repel(
      data = dat.rep,
      aes(label = rownames(dat.rep)),
      max.overlaps = 20,
      size = 4,
      box.padding = unit(0.5, "lines"),
      min.segment.length = 0,
      point.padding = unit(0.8, "lines"), segment.color = "black", show.legend = FALSE )+
    labs(x = "log2(Fold Change)",y = "-log10 (adj.P.Val)",  ### P.Value ， adj.P.Val
         title = paste0(disease," vs ", control),
         subtitle = paste(sprintf('adj.P.Val: %.2f;', pvalue.cutoff), 
                          sprintf('log2 (FoldChange): %.2f;', logfc.cutoff),
                          sprintf('Up: %1.0f; Down: %1.0f;', nrow(up_DEG), nrow(down_DEG)),
                          sprintf('Total: %1.0f', nrow(sig.diff)))) 
  volcano_plot
  ggsave(paste0('03.volcano(',gse,').pdf'), volcano_plot, width = 8, height = 6)
  ggsave(paste0('03.volcano(',gse,').png'), volcano_plot, width = 8, height = 6, dpi = 600, units = "in")
}

##---- different gene heatmap ----####
hmap.dat <- data
hmap.dat <- hmap.dat[rownames(dat.rep), ]
hmap.group <- group[order(group$group, decreasing = T),]
hmap.dat <- hmap.dat[ ,hmap.group$sample]
hmap.dat <- na.omit(hmap.dat)

# scale data
plot.dat <- t(scale(t(hmap.dat)))

# assignment -2 to 2
plot.dat[plot.dat < (-2)] <- (-2)
plot.dat[plot.dat > 2] <- 2


pdf(paste0('04.heatmap(',gse,').pdf'), width = 7, height = 6, family='Times')
densityHeatmap(plot.dat, title = "Distribution as heatmap", ylab = " ", height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = hmap.group$group, col = list(Group = c("Disease" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(plot.dat,
          # row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          name = "expression",
          cluster_rows = T,
          height = unit(10, "cm"),
          col = colorRampPalette(c("#104680", "white","#D80305"))(100))
dev.off()

png(paste0('04.heatmap(',gse,').png'), width = 7, height = 6, units = 'in', res = 600, family = 'Times')
densityHeatmap(plot.dat, title = "Distribution as heatmap", ylab = " ", height = unit(3, "cm")) %v%
  HeatmapAnnotation(Group = hmap.group$group, col = list(Group = c("Disease" = "#B72230", "Control" = "#104680"))) %v%
  Heatmap(plot.dat,
          row_names_gp = gpar(fontsize = 9),
          show_column_names = F,
          show_row_names = T,
          name = "expression",
          cluster_rows = T,
          height = unit(10, "cm"),
          col = colorRampPalette(c("#104680", "white", "#D80305"))(100))
dev.off()



# ###cir.heatmap-------------------------
# 
# library(ggnewscale)
# library(reshape2)
# dat.rep <- dat.rep[order(dat.rep$change),]
# 
# group$group = factor(group$group, levels = c("Control", "Disease"))
# group <- group[order(group$group,decreasing = T),]
# 
# 
# dat <- data
# diff <- dat[rownames(dat.rep),group$sample] %>% data.frame()
# diff <- data.frame(t(scale(t(diff))))
# diff[diff < (-2)] <- (-2)
# diff[diff > 2] <- 2
# 
# 
# annotation_col <- data.frame(group=group$group)
# rownames(annotation_col) <- colnames(diff)
# annotation_col$group <- factor(annotation_col$group,levels = c("Control", "Disease"))
# 
# annotation_row <- data.frame(dat.rep$change)
# rownames(annotation_row) <- rownames(dat.rep)
# colnames(annotation_row) <- " "
# 
# diff$gene <- rownames(diff)
# diff <- reshape2::melt(diff)
# 
# 
# res <- diff %>% dplyr::filter(variable  == colnames(dat)[1])
# 
# 
# {
#   res$ang<-NA
#   res$ang[1]<-30
#   for (i in (2:nrow(annotation_row))){
#     res$ang[i]<-  res$ang[i-1]-c(360/24)
#   }
#   
#   res$hjust <- 0
#   res$hjust[which(res$ang < -90)] <- 1
#   res$ang[which(res$ang < -90)] <- (180+res$ang)[which(res$ang < -90)]
#   
#   diff$var <- rep(1:nrow(annotation_row),nrow(annotation_col))
#   range(diff$value)
#   median(diff$value)
#   
# }
# 
# 
# p <- ggplot() +
#   geom_bar(data = annotation_col,stat = 'identity',
#            aes(x = 0.25,y = 1,fill = group),
#            width = 0.5,
#            color = NA) +
#   scale_fill_manual(name = 'group',
#                     values = c(Disease ="#ff9200",Control="#52b9d8"))+
#   new_scale("fill") +
#   geom_tile(data = diff[which(diff$variable == group$sample[1]),],
#             aes(x = 1:nrow(annotation_row),y = 0.5,fill = diff[which(diff$variable == group$sample[1]),]$value),
#             color = 'white')
# 
# for (i in 2:nrow(annotation_col)) {
#   p <- p+geom_tile(data = diff[which(diff$variable == group$sample[i]),],
#                    aes_string(x = 1:nrow(annotation_row),y = i-0.5,
#                               fill = diff[which(diff$variable == group$sample[i]),]$value),
#                    color = 'white')
# }
# p<- p+
#   # scale_fill_gradientn(colors = c(colorRampPalette(colors = c("#3C8DAD","white"))(length(seq(0.191,0.451,by=0.001))),
#   # colorRampPalette(colors = c("white","#FF6767"))(length(seq(0.452,4.324,by=0.001))))) +
#   scale_fill_gradient2(name=" ",
#                        midpoint = median(diff$value),
#                        low = 'blue3',
#                        mid = "white",
#                        high = 'red') +
#   coord_polar(theta = 'x') +
#   theme_void() +
#   xlim(-4,nrow(annotation_row)+1) +
#   theme(legend.text = element_text(size = 12),
#         legend.title = element_text(size = 14),
#         plot.margin = margin(10, 10, 10, 10))+
#   geom_text(data = res,
#             aes(x = as.numeric(rownames(res)),
#                 y = nrow(annotation_col)+1,
#                 label = gene, angle = ang, hjust = hjust),
#             size = 3)+
#   ylim(-8,nrow(annotation_col)+1)+ 
#   theme(panel.background = element_rect(fill = "white", colour = NA),
#         plot.background = element_rect(fill = "white", colour = NA))
# 
# ggsave(paste0('04.Train_cir.heatmap(',gse,').pdf'), p, width = 9, height = 7)
# ggsave(paste0('04.Train_cir.heatmap(',gse,').png'), p, width = 8, height = 7)
# 




