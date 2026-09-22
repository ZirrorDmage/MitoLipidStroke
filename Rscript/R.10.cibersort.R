
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

wd <- file.path(root,'10.cibersort')
dir.create(wd,showWarnings = F,recursive = T)
setwd(wd)



###---- load packages ----####
library(ggpubr); library(IOBR); library(magrittr); library(RColorBrewer); library(rstatix); library(reshape2)
library(corrplot); library(psych); library(Hmisc); library(ggplot2); library(ggcor); library(ggpubr); library(igraph)
library(vegan)
library(tidyverse)
library(lance)
###---- load data ----####
train.data <- read.csv('../00.rawdata/GSE16561/01.data_GSE16561.csv',row.names = 1) %>% lc.tableToNum()

# read group data file
group <- read.csv('../00.rawdata/GSE16561/02.group_GSE16561.csv',row.names = 1)
colnames(group) <- c('id', 'group')

# check sample
dat <- train.data[, group$id]
exprSet <- dat

###---- cibersort anlysis ----####
if (file.exists("result.rds")) {
  
  cibersort.res <- readRDS("result.rds")
} else {
  cibersort <- deconvo_tme(eset = exprSet, method = "cibersort", arrays = F, perm = 200,reference = lm22)
  
  # extract cibersot result and organize result data
  cibersort.res <- cibersort %>% as.data.frame()
  colnames(cibersort.res) %<>% gsub("_CIBERSORT","",.)
  colnames(cibersort.res) %<>% gsub("_"," ",.)
  rownames(cibersort.res) <- cibersort.res$ID
  cibersort.res <- cibersort.res[,-c(1,24:26)]
  colnames(cibersort.res) <- c("naive B cells",
                               "memory B cells",
                               "Plasma cells",
                               "CD8 T cells",
                               "naive CD4 T cells",
                               "resting CD4 memory T cells",
                               "activated CD4 memory T cells",
                               "follicular helper T cells",
                               "regulatory T cells (Tregs)",
                               "gammadelta T cells",
                               "resting NK cells",
                               "activated NK cells",
                               "Monocytes",
                               "M0 Macrophages",
                               "M1 Macrophages",
                               "M2 Macrophages",
                               "resting dendritic cells",
                               "activated dendritic cells",
                               "resting Mast cells",
                               "activated mast cells",
                               "Eosinophils",
                               "Neutrophils")
  # save result file 
  write.csv(cibersort.res,'01.cibersort_result.csv')
  saveRDS(cibersort.res,"result.rds")
}



###---- cibersort stacked diagrams ----####
col_sums <- colSums(cibersort.res)
sum <- data.frame(cell=names(col_sums),
                  total=col_sums)
sum <- sum%>%arrange(desc(sum$total))
cat(sum$cell[1:3],sep = "、")
# CD8 T cells、resting CD4 memory T cells、Neutrophils

# stack data
stack.dat <- t(cibersort.res) %>% as.data.frame()
stack.dat$cell_type <- rownames(stack.dat)

# stacked diagrams colors set
mypalette <- colorRampPalette(brewer.pal(9,"Set1"))

# cibersort stacked
stack.dat %>% gather(id, fraction, -cell_type) %>%  merge(group, by = 'id') %>%
  ggplot(aes(x = id, y = fraction, fill = cell_type)) +
  geom_bar(position = 'stack',stat = 'identity')+
  scale_y_continuous(expand = c(0,0))+
  theme_bw()+
  labs(x = '', y = 'Relative Percent', fill = '')+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = 'top') +
  scale_fill_manual(values = mypalette(22))+
  facet_grid(~ group, scales = "free", space = "free")

# extract plot
save.p <- recordPlot()

# save fig
pdf('02.cibersort_stack.pdf',width = 10,height = 8, family = 'Times')
replayPlot(save.p)
dev.off()

png('02.cibersort_stack.png',w=10,h=8,units='in',res=600, family = 'Times')
replayPlot(save.p)
dev.off()

###---- cibersort wilcox test ----####
wil.dat <- cibersort.res[, which(colSums(cibersort.res) > 0)]

# order
wil.dat <- wil.dat[group$id,order(colnames(wil.dat),decreasing = F)]
identical(group$id, rownames(wil.dat))

# merge data
wil.dat <- wil.dat %>% as.data.frame %>% tibble::rownames_to_column(var = "id")
identical(wil.dat$id, group$id)

# merge wilcox result and group
wil.dat <- merge(wil.dat, group, by = "id")

# wilcox test anlysis
dat.test <- tidyr::gather(wil.dat, Cell, `Cell composition`, -c(id, group))
colnames(dat.test)
stat.res <- dat.test %>%
  group_by(Cell) %>%
  wilcox_test(`Cell composition` ~ group) %>%
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")
stat.res

# ouput result
write.csv(stat.res, file = "03.wilcoxon.csv", row.names = F, quote = F)

# extract p < 0.05 cells
dif.cell <- stat.res[stat.res$p < 0.05, ]
paste0("共获得 ", nrow(dif.cell)," 个差异免疫细胞")
# "共获得11个差异免疫细胞"
cat(dif.cell$Cell,sep = "（）、")
# CD8 T cells（）、M0 Macrophages（）、M2 Macrophages（）、Monocytes（）、Neutrophils（）、activated mast cells（）、follicular helper T cells（）、gammadelta T cells（）、naive B cells（）、regulatory T cells (Tregs)（）、resting dendritic cells

###---- cibersort wilcox test box plot ----####
# organize data
if(nrow(dif.cell) > 4){
  box.dat <- dat.test[dat.test$Cell %in% dif.cell$Cell, ]
  box.cell <- dif.cell
  angle.set <- 45
  hjust.set <- 1
}else{
  box.dat <- dat.test
  box.cell <- stat.res
  angle.set <- 60
  hjust.set <- 1
}

# box plot
if(T){
  p <- ggboxplot(box.dat, x = "Cell", y = "Cell composition", fill = "group",
                 palette = c("#4f87ff",  "#C74546"), outlier.size = 0.8, bxp.errorbar = T) +
    ylim(0,.6)+
    stat_pvalue_manual(box.cell , x = "Cell", y.position = .6, size = 3.5, color = "black",
                       family = "Times", label = "p.signif") +
    theme_bw()+
    labs(title = "Differential immune cells", x = "", y = "Cell composition", color = "") +
    theme(plot.title = element_text(hjust = 0.5, colour = "black", face = "bold", size = 16),
          axis.text.x = element_text(angle = angle.set, hjust = hjust.set, colour ="black", face = "bold", size = 12),
          axis.text.y = element_text(hjust = 0.5, colour = "black", face = "bold", size = 12),
          axis.title.x = element_text(size = 16,face = "bold"),
          axis.title.y = element_text(size = 16,face = "bold"),
          legend.text = element_text(face = "bold", hjust = 0.5, colour = "black", size = 14),
          legend.title = element_text(face = "bold", size = 14),
          text = element_text(family = 'Times'),
          legend.position = "top",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())
  p
}

# ouput box plot
ggsave(filename = "04.differential_immune_cells.pdf", p, width = 12, height = 7)
ggsave(filename = "04.differential_immune_cells.png", p, width = 12, height = 7, dpi = 600, bg = "white")

#
#keygene&cell correlation---------
library(ggplot2)
library(ggcor)
library(vegan)
library(psych)
library(Hmisc)

sample <- group$id
group <- subset(group,sample %in% rownames(cibersort.res))

cell.dat <- wil.dat %>% column_to_rownames(var = 'id')
score2 <- cell.dat[,dif.cell$Cell]

hubgene <- read.csv('../06.expressionROC/keygenes.csv')[,1]
dat1 <- dat[,sample]
hub.dat1 <- dat1[hubgene,]%>%lc.tableToNum()
hub.dat1 <- t(hub.dat1)%>%as.data.frame()

score2 <- score2[rownames(hub.dat1),]
cor.dat <- cbind(score2,hub.dat1 )

cor.res <- cor(cor.dat, method = "spearman")
corp <- cor_pmat(cor.dat)
write.csv(cor.res, '05.gene_correlation.csv', quote = F)
write.csv(corp, '05.gene_correlation_p.csv', quote = F)

if(nrow(cor.res) > 10){
  width.set <- 14
  height.set <- 14
}else{
  width.set <- 10
  height.set <- 10
  
}


if(T){
  cell_num <- nrow(dif.cell)
  gene_num <- length(hubgene)
  mycolors <- c("#8db3dd" ,  "#827db3", "white", "#ef99a6",  "#ae2213")
  corr_plot <- quickcor(cor.dat, cor.test = TRUE,method = "spearman",
                        insig = "blank",
                        outline = "white",
                        addCoef.col ="black",
                        col = col1(100),
                        tl.col = 'black',
                        tl.offset = 0.4,
                        number.font = 2,
                        cl.cex = 0.8,
                        number.cex = 1.2,
                        show.diag = TRUE,
                        type = c("lower", "upper", "full")[1]) +
    geom_square(data = get_data(type = "lower", show.diag = T)) +
    geom_mark(data = get_data(type = "lower",abs(r) > 0.3, show.diag = F), color = "black")+
    scale_fill_gradientn(colours = mycolors, limits = c(-1,1))+  ###color
    labs(title = "Genes and Cells Spearman correlation",
         fill = "cor")+
    theme(
      plot.title = element_text(color = "black", size = 20 ,face = "bold", hjust = .5),
      axis.text.x = element_text(color = c(rep(mycolors[2], cell_num), rep(mycolors[length(mycolors)], gene_num)),
                                 #face = "bold",
                                 size = 14),
      axis.text.y = element_text(color = c(rep(mycolors[length(mycolors)], gene_num), rep(mycolors[2], cell_num)),
                                 #face = "bold",
                                 size = 14),
      legend.title = element_text(size = 13, color = "black", face = "bold"),
      legend.text = element_text(size = 12, color = "black"),
      legend.position = "right"
    )
  corr_plot
  ggsave("06.correlation.pdf", corr_plot, height = height.set, width = width.set)
  ggsave("06.correlation.png", corr_plot, height = height.set, width = width.set, dpi = 600, units = "in")
}


