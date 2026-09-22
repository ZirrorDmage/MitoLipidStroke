
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
wd <- file.path(root,'02.venn')
dir.create(wd,showWarnings = F,recursive = T)
setwd(wd)


####---- load packages ----####
library(ggvenn); library(ggplot2); library(tidyr); library(data.table)
library(VennDiagram)

###---- function ----###
venn.plot <- function(data){
  color.set <- c("#EA9C9D","#FEBF0E","#103778", "#d64c39", "#49afca", "#F26800")
  name.color.set <- c("#DE5A5A", "#FF8204","#0897B4", "#a6666b", "#6a84a0", "#6B4945")
  venn.p <- ggvenn(data,
                   fill_color = color.set[1:length(data)],
                   show_percentage = T,
                   stroke_alpha = 0.5,
                   digits = 2,
                   stroke_size = 1,
                   text_size = 5,
                   stroke_color = "white",
                   stroke_linetype = "dashed",
                   set_name_color = name.color.set[1:length(data)],
                   set_name_size = 6,
                   text_color = 'black')
  return(venn.p)
}

####---- anlysis ----####

library(readxl)
excel_file <- "../00.rawdata/线粒体和脂质代谢相关基因.xlsx"
tmp <- readxl::read_xlsx(excel_file, sheet = 1)
LMRGs <- unique(tmp$`lipid metabolism-related genes`) ## 847

tmp <- readxl::read_xlsx(excel_file, sheet = 2)
MMRGs <- unique(tmp$`Mitochondrial Genes (n=2030)`)## 2030

LM_MM <- intersect(LMRGs,MMRGs)## 253

write.csv(LM_MM,"MLM related genes list.csv")

DEGs <- read.csv("../01.different_genes/02.DEGs_sig(GSE16561).csv", row.names = 1) %>%
  mutate(symbol = rownames(.)) %>% 
  # filter(abs(log2FoldChange)>0.5, padj<0.05) %>% 
  pull(symbol)

candi <- intersect(LM_MM, DEGs)




####
venn.ls <- list(
  "DEGs" = DEGs, 
  "MLMRGs" = LM_MM
)

venn <- venn.plot(venn.ls)
venn


ggsave(filename = "02.venn.pdf", venn, width = 6, height = 6, bg = "white")
ggsave(filename = "02.venn.png", venn, width = 6, height = 6, bg = "white", dpi = 600, units = "in")

inter <- get.venn.partitions(venn.ls)
for (i in 1:nrow(inter)) inter[i,'values'] <- paste(inter[[i,'..values..']], collapse = ', ')
cadigenes <- data.frame(symbol = inter[1,]$..values.. $`1`)

write.csv(cadigenes, file = "02.cadi.genes.csv", row.names = F, quote = F)
dev.off()

