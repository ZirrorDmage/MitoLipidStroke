
rm(list=ls())
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
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

dir <- "13.drugs"
if(!dir.exists(dir)) { dir.create(dir)}
setwd(dir)

###---- load packages ----####
library(enrichR)
# BiocManager::install('enrichR')
library(ggalluvial); library(ggplot2); library(RColorBrewer);library(dplyr);library(stringr)
# library(networkD3)

###---- drug ----####
dbs <- listEnrichrDbs()
dbs <- c("DSigDB")#


symbols <-  read.csv('../06.expressionROC/keygenes.csv')[,1]
enrichr <- enrichr(symbols, dbs)

result <- data.frame(enrichr$DSigDB)
write.csv(result,"raw.drugs_result.csv", row.names = F, quote = T)

result <- result[which(result$Adjusted.P.value < 0.05 & result$Combined.Score>1000),] ##$Combined.Score>1000
result <- result[order(result$P.value),]
table(result$Genes)



# 
# results_list <- list()
# for (symbol in symbols) {
#   enrichr <- enrichr(symbol, dbs)
#   result <- data.frame(enrichr$DSigDB)
#   result <- result[which(result$P.value < 0.05), ]#&result$Combined.Score>1000
#   result <- result[order(result$P.value), ]
#   results_list[[symbol]] <- result
# }
# 
# result2 <- do.call(rbind, results_list)
# table(result2$Genes)


plot_result <- result
plot_result <- tidyr::separate_rows(plot_result, Genes)
table(plot_result$Genes)
# CYP1B1   DEGS1   HSDL2 OSBPL1A    PTEN 
# 64       1       1       1      13 
length(unique(plot_result$Term)) ##=74


plot_result <- plot_result %>%
  mutate(Drug_name = if_else(grepl("CTD|BOSS|MCF7|PC3|HL60|TTD", Term), 
                             str_extract(Term, "^.*?(?= (CTD|BOSS|MCF7|PC3|HL60|TTD))"), 
                             Term))
write.csv(plot_result, file = '01.drugs_result.csv', row.names = F, quote = T)




plot.data <- data.frame(genes = plot_result$Genes, drugs = plot_result$Drug_name)%>% unique()

sample_cols <- unique(c(brewer.pal(8,"Dark2"),brewer.pal(9,"Set1"),brewer.pal(8,"Set2"),brewer.pal(12,"Set3")))
colors <- c("#d0bbf8", "#9cbcf7","#f1c97d","#72c2c5","#f1adad","#d3d47b","#f5b8f5","#BC8F8F",sample_cols,
            "orange","#f1f177","#8febe6","#65a10a",sample_cols,
            "#f1c97d","#b40808", sample_cols)


p <- ggplot(plot.data, aes(axis1 = genes, axis2 = drugs)) +
  geom_alluvium(aes(fill = genes), width = 0.5, alpha = 0.8, color = NA) +
  geom_stratum(aes(fill = after_stat(stratum)), width = 0.5, color = NA) +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)),
            size = 4, color = "black") +
  scale_fill_manual(values = colors) +  # color
  scale_x_discrete(limits = c("Genes", "Drugs"), expand = c(0.05, 0.05)) +
  labs(title = "Gene-Drug Interaction",
       subtitle = "Interactions between genes and drugs") +
  theme_minimal(base_size = 14) +  # Base font size
  theme(axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.text.x = element_text(colour = "black", size = 14, face = "bold"),
        axis.text.y = element_blank(),
        legend.position = "none",
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        plot.subtitle = element_text(size = 14, hjust = 0.5),
        plot.caption = element_text(size = 10, hjust = 1),
        panel.border = element_blank(),  # Remove the border
        panel.grid.major = element_blank(),  # Remove the main grid lines
        panel.grid.minor = element_blank())
p
ggsave("02.drugs_sankey.pdf", width = 11, height = 12)
ggsave("02.drugs_sankey.png", width = 11, height = 12, units = "in", dpi = 600, bg = "white")


