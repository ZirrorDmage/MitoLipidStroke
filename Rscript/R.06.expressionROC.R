

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


wd <- file.path(root,"06.expressionROC")
dir.create(wd,recursive = T,showWarnings = F)
setwd(wd)



pdf.options(family="Times")
library("ggpubr")
library("reshape2")
library(dplyr)
library(rstatix)


###############train-GSE16561----------------------
train_id <- "GSE16561"
train.dat <- read.csv('../00.rawdata/GSE16561/01.data_GSE16561.csv',row.names = 1, check.names = F) 
range(train.dat)

group1 <- read.csv('../00.rawdata/GSE16561/02.group_GSE16561.csv',row.names = 1)
train.dat <- train.dat[,group1$sample]
geneRT <- read.csv("../05.machine_learning/Best_model_gene.csv",row.names = 1)[,1]


dat1 <- na.omit( train.dat[geneRT, group1$sample] )%>% t() %>% as.data.frame()

dat1 <-  merge(dat1, group1, by.x = 'row.names', by.y = 'sample') %>% tibble::column_to_rownames(var = 'Row.names')

table(dat1$group)
dat1$group <- factor(dat1$group, levels = c('Control', 'Disease'))

train <-  dat1

colnames(train)

train_data=reshape2::melt(train,id.vars=c("group"))

colnames(train_data)=c("Type","Gene","Expression")

p_values1 <- train_data %>%
  group_by(Gene) %>%
  summarise(
    p_value = wilcox.test(Expression ~ Type)$p.value,
    median_dis = median(Expression[Type == "Disease"],  na.rm = TRUE),
    median_con = median(Expression[Type == "Control"], na.rm = TRUE)
  ) %>%
  mutate(
    change = case_when(
      p_value >= 0.05                     ~ NA_character_,
      median_dis >  median_con            ~ "UP",
      median_dis <  median_con            ~ "DOWN",
      TRUE                                ~ "NoChange"   # 
    )
  )
sig1 <- p_values1[p_values1$p_value<0.05,]%>%as.data.frame()


stat.res <- train_data %>% 
  group_by(Gene) %>% 
  wilcox_test(`Expression` ~ Type) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

p1 <- ggboxplot(train_data, x = "Gene", y = "Expression", fill = "Type",
                palette = c("#1ba784", "#fb8b05"), outlier.shape = NA, outlier.size = 0.8, bxp.errorbar = T) +
  # ylim(0,.6)+
  labs(y =paste0("Expression of ",train_id) , x = NULL, fill = "Type") +
  stat_pvalue_manual(stat.res , x = "Gene", y.position = max(train_data$Expression)*1.02, size = 3.5, color = "black",
                     family = "Times", label = "p.signif") +
  theme(
    legend.title = element_text(color = "black"),  # 
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),  # 
    panel.background = element_rect(fill = "white"),  # background
    panel.grid.major = element_blank(),  # grid
    panel.grid.minor = element_blank(),  # grid
    axis.line = element_line(color = "black"),  # axis
    axis.ticks = element_line(color = "black") ,
    legend.position = "right"
  ) 

pdf(file=paste0("01.train_",train_id,"_expr.pdf"), width=9, height=5,family = 'Times')
print(p1)
dev.off()

png(file=paste0("01.train_",train_id,"_expr.png"), width=9, height=5,family = 'Times',units = 'in',res = 600)
print(p1)
dev.off()


###############GSE58294------------------------
test_id <- "GSE58294"
test.dat <- read.csv('../00.rawdata/GSE58294/01.data_GSE58294.csv',row.names = 1) 
group2 <- read.csv('../00.rawdata/GSE58294/02.group_GSE58294.csv',row.names = 1)
test.dat <- test.dat[,group2$sample]

setdiff(geneRT,rownames(test.dat))

dat2 <-  test.dat[geneRT, group2$sample] %>% t() %>% as.data.frame()


dat2 <-  merge(dat2, group2, by.x = 'row.names', by.y = 'sample') %>% tibble::column_to_rownames(var = 'Row.names')

table(dat2$group)
dat2$group <- factor(dat2$group, levels = c('Control', 'Disease'))

test <-  dat2

colnames(test)

test_data=reshape2::melt(test,id.vars=c("group"))

colnames(test_data)=c("Type","Gene","Expression")


stat.res <- test_data %>% 
  group_by(Gene) %>% 
  wilcox_test(`Expression` ~ Type) %>% 
  adjust_pvalue(method = "BH") %>%  # method BH == fdr
  add_significance("p")

p2 <- ggboxplot(test_data, x = "Gene", y = "Expression", fill = "Type",
                palette = c("#1ba784", "#fb8b05"), outlier.shape = NA, outlier.size = 0.8, bxp.errorbar = T) +
  # ylim(0,.6)+
  labs(y =paste0("Expression of ",test_id) , x = NULL, fill = "Type") +
  stat_pvalue_manual(stat.res , x = "Gene", y.position = max(test_data$Expression)*1.02, size = 3.5, color = "black",
                     family = "Times", label = "p.signif") +
  theme(
    legend.title = element_text(color = "black"),  # 
    axis.text.x = element_text(angle = 0, hjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),  # 
    panel.background = element_rect(fill = "white"),  # background
    panel.grid.major = element_blank(),  # grid
    panel.grid.minor = element_blank(),  # grid
    axis.line = element_line(color = "black"),  # axis
    axis.ticks = element_line(color = "black") ,
    legend.position = "right"
  ) 

pdf(file=paste0("02.test_",test_id,"_expr.pdf"), width=9, height=5,family = 'Times')
print(p2)
dev.off()

png(file=paste0("02.test_",test_id,"_expr.png"), width=9, height=5,family = 'Times',units = 'in',res = 600)
print(p2)
dev.off()


p_values2 <- test_data %>%
  group_by(Gene) %>%
  summarise(
    p_value = wilcox.test(Expression ~ Type)$p.value,
    median_dis = median(Expression[Type == "Disease"],  na.rm = TRUE),
    median_con = median(Expression[Type == "Control"], na.rm = TRUE)
  ) %>%
  mutate(
    change = case_when(
      p_value >= 0.05                     ~ NA_character_,
      median_dis >  median_con            ~ "UP",
      median_dis <  median_con            ~ "DOWN",
      TRUE                                ~ "NoChange"   # 
    )
  )
sig2 <- p_values2[p_values2$p_value<0.05,]%>%as.data.frame()


sig12 <- sig1%>%filter(Gene%in%sig2$Gene)


common_genes12 <- sig12 %>%
  inner_join(
    sig2 %>% dplyr::select(Gene, change),
    by = "Gene",
    suffix = c("_test", "_train")
  ) %>%
  filter(change_train == change_test) %>%
  dplyr::select(Gene)


keygenes <- data.frame(symbol=as.character(common_genes12$Gene))
cat(keygenes$symbol,sep="、")

# write.csv(keygenes,"keygenes.csv",row.names = F,quote = F)






##############trainROC############
library(pROC)
library(ggplot2)
library(dplyr)

dat1 <- dat1 %>%
  mutate(y_bin = ifelse(group == "Disease", 1, 0))
dat1$y_bin <- factor(dat1$y_bin)
pred_vars <-keygenes$symbol  # 

roc_list <- lapply(pred_vars, function(v) {
  roc(response = dat1$y_bin,
      predictor = dat1[[v]],
      quiet = TRUE)
})


auc_values <- sapply(roc_list, function(roc_obj) roc_obj$auc)

# 创建一个新的数据框，包含基因名和对应的AUC值
Trainauc_df <- data.frame(
  Gene = pred_vars,
  AUC = auc_values,
  stringsAsFactors = FALSE
)
Trainauc_df2 <- Trainauc_df %>% filter(AUC>0.7&AUC<1)




col_set   <- c("#00B4A0FF","#0072B5FF","#E18727FF","#DE5A5A","#9F9AC7")


pdf(sprintf("03.train_ROC_%s.pdf",train_id), width = 7, height = 6)
plot(roc_list[[1]], col = col_set[1], lwd = 2.5,
     legacy.axes = TRUE, grid = TRUE, grid.col = "gray85",
     main = "", xlab = "FPR", ylab = "TPR")
for (i in 2:nrow(Trainauc_df)) {
  plot(roc_list[[i]], add = TRUE, col = col_set[i], lwd = 2.5)
}

aucs <- sapply(roc_list, function(r) round(auc(r), 3))
legend("bottomright",
       legend = sprintf("AUC of %s = %.3f", pred_vars, aucs),
       col    = col_set,
       lty    = 1,
       lwd    = 2.5,
       bty    = "n",
       cex    = 1.1)

dev.off()

png(sprintf("03.train_ROC_%s.png",train_id), width = 7, height = 6,units = 'in',res = 600)
plot(roc_list[[1]], col = col_set[1], lwd = 2.5,
     legacy.axes = TRUE, grid = TRUE, grid.col = "gray85",
     main = "", xlab = "FPR", ylab = "TPR")
for (i in 2:nrow(Trainauc_df)) {
  plot(roc_list[[i]], add = TRUE, col = col_set[i], lwd = 2.5)
}

aucs <- sapply(roc_list, function(r) round(auc(r), 3))
legend("bottomright",
       legend = sprintf("AUC of %s = %.3f", pred_vars, aucs),
       col    = col_set,
       lty    = 1,
       lwd    = 2.5,
       bty    = "n",
       cex    = 1.1)

dev.off()



##############testROC############
library(pROC)
library(ggplot2)
library(dplyr)

dat2 <- dat2 %>%
  mutate(y_bin = ifelse(group == "Disease", 1, 0))
dat2$y_bin <- factor(dat2$y_bin)


roc_list <- lapply(pred_vars, function(v) {
  roc(response = dat2$y_bin,
      predictor = dat2[[v]],
      quiet = TRUE)
})

auc_values <- sapply(roc_list, function(roc_obj) roc_obj$auc)
Testauc_df <- data.frame(
  Gene = pred_vars,
  AUC = auc_values,
  stringsAsFactors = FALSE
)

Testauc_df2 <- Testauc_df %>% filter(AUC>0.7&AUC<1)

keygenes2 <- data.frame(symbol=sort(intersect(Testauc_df2$Gene,Trainauc_df2$Gene)))

write.csv(keygenes2,"keygenes.csv",row.names = F,quote = F)

cat(keygenes2$symbol,sep = "、")
# CYP1B1、DEGS1、HSDL2、OSBPL1A、PTEN

pdf(sprintf("04.test_ROC_%s.pdf",test_id), width = 7, height = 6)
plot(roc_list[[1]], col = col_set[1], lwd = 2.5,
     legacy.axes = TRUE, grid = TRUE, grid.col = "gray85",
     main = "", xlab = "FPR", ylab = "TPR")
for (i in 2:nrow(Trainauc_df)) {
  plot(roc_list[[i]], add = TRUE, col = col_set[i], lwd = 2.5)
}

aucs <- sapply(roc_list, function(r) round(auc(r), 3))
legend("bottomright",
       legend = sprintf("AUC of %s = %.3f", pred_vars, aucs),
       col    = col_set,
       lty    = 1,
       lwd    = 2.5,
       bty    = "n",
       cex    = 1.1)

dev.off()

png(sprintf("04.test_ROC_%s.png",test_id), width = 7, height = 6,units = 'in',res = 600)
plot(roc_list[[1]], col = col_set[1], lwd = 2.5,
     legacy.axes = TRUE, grid = TRUE, grid.col = "gray85",
     main = "", xlab = "FPR", ylab = "TPR")
for (i in 2:nrow(Trainauc_df)) {
  plot(roc_list[[i]], add = TRUE, col = col_set[i], lwd = 2.5)
}

aucs <- sapply(roc_list, function(r) round(auc(r), 3))
legend("bottomright",
       legend = sprintf("AUC of %s = %.3f", pred_vars, aucs),
       col    = col_set,
       lty    = 1,
       lwd    = 2.5,
       bty    = "n",
       cex    = 1.1)

dev.off()



