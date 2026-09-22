


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
wd <- file.path(root,'08.location')
dir.create(wd,showWarnings = F,recursive = T)
setwd(wd)


####---- Subcellular localization ----####
rm(list = ls())
# https://services.healthtech.dtu.dk/services/DeepLoc-2.1/

loc <- read.csv("01.subcell.results_DeepLoc2.1.csv") 
colnames(loc)
loc <- loc%>%dplyr::select(-c("Localizations", "Signals","Membrane.types",
                              "Plastid","Peroxisome",
                              "Peripheral" , "Transmembrane", "Lipid.anchor", "Soluble"))


library(tidyverse)
library(ggplot2)


dat <- loc %>%
  pivot_longer(cols = -Protein_ID,             
               names_to = "Localization",
               values_to = "Probability") %>%
  mutate(Protein_ID = factor(Protein_ID))

p_barh <- ggplot(dat, aes(x = Probability,
                          y = Localization,
                          fill = Protein_ID)) +
  geom_col(position = "dodge", width = 0.9) +
  scale_fill_manual(values = c("#F49897", "#b5d66b", "#637bbd", "#FFB716",  "#90d3c1"),
                    name = "Gene") +
  labs(title = "Subcellular Localization",
       x = "Probability",
       y = NULL) +
  scale_x_continuous(limits = c(0, 1),
                     expand = expansion(mult = c(0, 0.02))) +
  theme_bw(base_family = "Times") +
  theme(
    plot.title    = element_text(hjust = 0.5, face = "bold", size = 16),
    axis.text.y   = element_text(face = "bold", size = 12),
    axis.text.x   = element_text(face = "bold", size = 12),
    axis.title.x  = element_text(face = "bold", size = 14),
    legend.title  = element_text(face = "bold", size = 12),
    legend.text   = element_text(face = "bold", size = 11),
    legend.position = "top",
    panel.grid    = element_blank()
  )

print(p_barh)


pdf("02.localization.pdf", width = 8, height = 6, family = 'Times')
print(p_barh)
dev.off()

png("02.localization.png", width = 8, height = 6, family = 'Times', units = 'in', res = 600)
print(p_barh)
dev.off()

