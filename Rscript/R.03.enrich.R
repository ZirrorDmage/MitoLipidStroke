
print(Sys.time())
rm(list = ls())
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
wd <- file.path(root,'03.enrich')
dir.create(wd,showWarnings = F,recursive = T)
setwd(wd)

rm(list=ls());gc()
# 如需指定本机 R 库路径，可在此处设置：.libPaths(c("你的R库路径", .libPaths()))
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(enrichplot)
library(stringr)
library(tidyr)
library(ggplot2)
library(colorspace)
library(circlize)
library(RColorBrewer)


candidatGene <- read.csv("../02.venn/01.cadi.genes.csv",header = T,check.names = F)

Gene = bitr(candidatGene$symbol,fromType = "SYMBOL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")  
Gene <- na.omit(Gene) 


ego <- enrichGO(gene = Gene$ENTREZID, 
                OrgDb = org.Hs.eg.db,
                keyType = "ENTREZID", 
                ont = "ALL",
                pAdjustMethod = "BH",
                pvalueCutoff = 0.2,
                qvalueCutoff = 0.2,
                readable = TRUE)

kk <- enrichKEGG(gene = Gene$ENTREZID,
                 keyType = "kegg",
                 organism = "hsa",
                 pAdjustMethod = "BH",
                 pvalueCutoff = 0.2, # adjusted pvalue cutoff on enrichment tests to report, default is 0.05
                 qvalueCutoff = 0.2
)
save.image("enrich.rdata")
kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")

dim(kk@result)
ekegg_all <- kk
ego_all <- ego


ego_readable <- setReadable(ego_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
GO <- as.data.frame(ego_readable)
go_result <- GO
go_result <- go_result[go_result$p.adjust<0.05,]
go_result <- go_result[order(go_result$p.adjust,decreasing = F),]
BP <- go_result[go_result$ONTOLOGY=='BP',]
CC <- go_result[go_result$ONTOLOGY=='CC', ]
MF <- go_result[go_result$ONTOLOGY=='MF', ]
paste0("得到",dim(go_result)[[1]],"个结果，其中",dim(BP)[[1]],"个生物过程，",dim(CC)[[1]],"个细胞组分，",dim(MF)[[1]],"个分子功能")
# 得到144个结果，其中94个生物过程，13个细胞组分，37个分子功能

write.csv(go_result,file = "01.go_ALL.csv",  row.names = T)

ekegg_readable <- setReadable(ekegg_all, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
KEGG <- as.data.frame(ekegg_readable)
kk_result <- KEGG
kk_result <- kk_result[kk_result$p.adjust<0.05,]
kk_result <- kk_result[order(kk_result$p.adjust,decreasing = F),]
nrow(kk_result)# 8
write.csv(kk_result,"02.KEGG-enrich.csv",row.names =F)



##plotdata
use_pathway <- group_by(go_result, ONTOLOGY) %>%
  top_n(5, wt = -p.adjust) %>%
  # group_by(p.adjust) %>%
  # top_n(1, wt = Count) %>%
  rbind(
    top_n(kk_result, 5, -p.adjust) %>%
      # group_by(p.adjust) %>%
      # top_n(1, wt = Count) %>%
      mutate(ONTOLOGY = 'KEGG')
  ) %>%
  ungroup() %>%
  mutate(ONTOLOGY = factor(ONTOLOGY, 
                           levels = rev(c('BP', 'CC', 'MF', 'KEGG')))) %>%
  dplyr::arrange(ONTOLOGY, p.adjust) %>%
  mutate(Description = factor(Description, levels = Description)) %>%
  tibble::rowid_to_column('index')

use_pathway <- use_pathway %>%
  mutate(geneID = sapply(str_split(geneID, "/"), function(x) paste(head(x, 20), collapse = "/")))

#leftplot
width <- 0.5
# x axis
xaxis_max <- max(-log10(use_pathway$p.adjust)) + 1

#  ONTOLOGY num
ontology_table <- table(use_pathway$ONTOLOGY)
ontology_df <- as.data.frame(ontology_table)
colnames(ontology_df) <- c("ONTOLOGY", "n")

# ymin & ymax
ontology_df$ymax <- cumsum(ontology_df$n)
ontology_df$ymin <- c(0, head(ontology_df$ymax, -1)) + 0.6
ontology_df$ymax <- ontology_df$ymax + 0.4

# and xmin & xmax
ontology_df$xmin <- -3 * width
ontology_df$xmax <- -2 * width

#  rect.data
rect.data <- ontology_df


library(ggprism)
pal <- c("#49afca", "#ee9b01", "#90d3c1",'#807DBC')


plot_enrichment <- function() {
  use_pathway_top_genes <- use_pathway %>%
    group_by(ONTOLOGY, Description) %>%
    slice_head(n = 5) %>%  
    ungroup()
  library(stringr)
  library(dplyr)
  
  use_pathway_top_genes <- use_pathway_top_genes %>%
    mutate(geneID = sapply(str_split(geneID, "/"),
                           function(x) paste0(head(x, 10), collapse = "/")))# max 10 genes
  p <- use_pathway_top_genes %>%
    ggplot(aes(-log10(p.adjust), y = index, fill = ONTOLOGY)) +
    geom_col(aes(y = Description), width = 0.6, alpha = 0.8) +
    geom_text(
      aes(x = 0.05, label = Description),
      hjust = 0, size = 5
    ) +
    geom_text(
      aes(x = 0.1, label = geneID, colour = ONTOLOGY),
      hjust = 0, vjust = 2.6, size = 3.5, fontface = 'italic',
      show.legend = FALSE
    ) +
    # gene number
    geom_point(
      aes(x = -width, size = Count),
      shape = 21
    ) +
    geom_text(
      aes(x = -width, label = Count)
    ) +
    scale_size_continuous(name = 'Count', range = c(3, 8)) +
    #  geom_rect  instead of  geom_round_rect
    geom_rect(
      aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
          fill = ONTOLOGY),
      data = rect.data,
      inherit.aes = FALSE
    ) +
    geom_text(
      aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = ONTOLOGY),
      data = rect.data,
      inherit.aes = FALSE,
      angle = 90,   # 
      hjust = 0.5,  # 
      vjust = 0.5,  # 
      size = 4
    ) +
    geom_segment(
      aes(x = 0, y = 0, xend = xaxis_max, yend = 0),
      linewidth = 1.5,
      inherit.aes = FALSE
    ) +
    labs(y = NULL) +
    scale_fill_manual(name = 'Category', values = pal) +
    scale_colour_manual(values = pal) +
    scale_x_continuous(
      breaks = seq(0, xaxis_max, 10),
      expand = expansion(c(0, 0))
    ) +
    theme_prism() +
    theme(
      axis.text.y = element_blank(),
      axis.line = element_blank(),
      axis.ticks.y = element_blank(),
      legend.title = element_text()
    )
  return(p)
}
enrichment_plot <- plot_enrichment()

# print(enrichment_plot)

ggsave("03.GO_KEGG_plot.pdf", width = length(use_pathway$ONTOLOGY)/4*3, height = length(use_pathway$ONTOLOGY)/20*9, enrichment_plot,family = "Times")
png("03.GO_KEGG_plot.png", width = length(use_pathway$ONTOLOGY)/4*3, height = length(use_pathway$ONTOLOGY)/20*9,units = "in",res = 600,family = "Times")
print(enrichment_plot)
dev.off()


