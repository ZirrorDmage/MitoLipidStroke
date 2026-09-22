

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
dir <- "00.rawdata"
if(!dir.exists(dir)){dir.create(dir)}
setwd(dir)


library(readr)
library(tidyr)
library(openxlsx)
library(dplyr)
library(readxl)
library(GEOquery)
library(tinyarray)
library(tidyverse)
library(data.table)

#########GSE58294----------------array#########
rm(list = ls())
data.id <- "GSE58294"
path <- file.path(getwd(), data.id)
if (!dir.exists(path)) dir.create(path)
###---- function ----###
# 检查是否需要 log2 转换
check_log_transform <- function(data) {
  
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  
  return((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2))
  
}

# 执行 log2 转换
log_transform <- function(data) {
  data[which(data <= 0)] <- 0
  return(log2(data + 1))
}

###---- 处理 GEO 数据 ----###
# get geo matrix
geo.set <- getGEO(data.id, destdir = path, GSEMatrix = T, getGPL = F)
expr.data <- Biobase::exprs(geo.set[[1]])

# 判断是否需要 log2 转换
if (check_log_transform(expr.data)) {
  expr.data <- log_transform(expr.data)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}
#"log2 transform not needed"


# 获取 GPL 注释数据

gpls <- geo.set[[1]]@annotation

gpl.raw <- getGEO(gpls, destdir = file.path(getwd(), "GPL"))
gpl <- Table(gpl.raw)
colnames(gpl)

gpl <- gpl %>%
  dplyr::select(ID, `Gene Symbol`) %>%
  dplyr::rename(id = ID, symbol = `Gene Symbol`) %>%
  separate(symbol, into = c("symbol", NA), sep = " /// ")
gpl <- gpl %>% filter(symbol != "")
# 匹配并处理表达矩阵数据
data <- expr.data[rownames(expr.data) %in% gpl$id,]
gpl <- gpl[match(rownames(data), gpl$id),]

data <- as.data.frame(data) %>%
  tibble::rownames_to_column(var = "rowname") %>%
  mutate(rowname = gpl[match(rowname, gpl$id), "symbol"]) %>%
  mutate(rowMean = rowMeans(dplyr::select(., starts_with("GSM")))) %>%
  arrange(desc(rowMean)) %>%
  distinct(rowname, .keep_all = TRUE) %>%
  dplyr::select(-rowMean) %>%
  {rownames(.) <- NULL; .} %>%
  tibble::column_to_rownames("rowname")

range(data)
data[data < 0] <- 0

# 分组信息
phen.data <- pData(geo.set[[1]])
# View(phen.data)
names(phen.data)
Group <- data.frame(sample = phen.data$geo_accession,
                    group = phen.data$`group:ch1`)

table(Group$group)

# library(stringr)
# Group$group <- str_extract(Group$group, "^\\w+")#leave the first word
# table(Group$group)

group <- Group %>%
  na.omit() %>%
  mutate(group = ifelse(group == "Control", "Control", "Disease"))

group <- group[order(group$group),]
head(group)

data <- data[, group$sample]

write.csv(data, sprintf("%s/01.data_%s.csv", path, data.id), row.names = T, quote = T)
write.csv(group, sprintf("%s/02.group_%s.csv", path, data.id), row.names = T, quote = T)







#########GSE16561----------------array#########
rm(list = ls())



library(readr)
library(tidyr)
library(openxlsx)
library(dplyr)
library(readxl)
library(GEOquery)
library(tinyarray)
library(tidyverse)
library(data.table)


rm(list = ls())
data.id <- "GSE16561"
path <- file.path(getwd(), data.id)
if (!dir.exists(path)) dir.create(path)
###---- function ----###
# 检查是否需要 log2 转换
check_log_transform <- function(data) {
  
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  
  return((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2))
  
}

# 执行 log2 转换
log_transform <- function(data) {
  data[which(data <= 0)] <- 0
  return(log2(data + 1))
}

###---- 处理 GEO 数据 ----###
# get geo matrix
geo.set <- getGEO(data.id, destdir = path, GSEMatrix = T, getGPL = F)
expr.data <- Biobase::exprs(geo.set[[1]])

expr_matrix <- as.matrix(expr.data)  
set.seed(12345)
library(impute)
expr <- impute.knn(expr_matrix)$data
expr.data <- data.frame(expr, check.rows = F, check.names = F)

# 判断是否需要 log2 转换
if (check_log_transform(expr.data)) {
  expr.data <- log_transform(expr.data)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}
#"log2 transform not needed"

# 获取 GPL 注释数据

gpls <- geo.set[[1]]@annotation

gpl.raw <- getGEO(gpls, destdir = file.path(getwd(), "GPL"))
gpl <- Table(gpl.raw)
colnames(gpl)

gpl <- gpl %>%
  dplyr::select(ID, Symbol) %>%
  dplyr::rename(id = ID, symbol = Symbol) 

gpl <- gpl %>% filter(symbol != "")
# 匹配并处理表达矩阵数据
data <- expr.data[rownames(expr.data) %in% gpl$id,]
gpl <- gpl[match(rownames(data), gpl$id),]

data <- as.data.frame(data) %>%
  tibble::rownames_to_column(var = "rowname") %>%
  mutate(rowname = gpl[match(rowname, gpl$id), "symbol"]) %>%
  mutate(rowMean = rowMeans(dplyr::select(., starts_with("GSM")))) %>%
  arrange(desc(rowMean)) %>%
  distinct(rowname, .keep_all = TRUE) %>%
  dplyr::select(-rowMean) %>%
  {rownames(.) <- NULL; .} %>%
  tibble::column_to_rownames("rowname")

range(data)
# data[is.na(data)] <- 0
# data[data < 0] <- 0

# 分组信息
phen.data <- pData(geo.set[[1]])
# View(phen.data)
names(phen.data)
Group <- data.frame(sample = phen.data$geo_accession,
                    group = phen.data$description)

table(Group$group)



group <- Group %>%
  na.omit() %>%
  mutate(group = ifelse(group == "Control", "Control", "Disease"))

group <- group[order(group$group),]
head(group)

data <- data[, group$sample]

write.csv(data, sprintf("%s/01.data_%s.csv", path, data.id), row.names = T, quote = T)
write.csv(group, sprintf("%s/02.group_%s.csv", path, data.id), row.names = T, quote = T)



#########GSE22255----------------array#########
# 备选数据集：早期筛查时一并下载，最终未纳入本分析
rm(list = ls())
data.id <- "GSE22255"
path <- file.path(getwd(), data.id)
if (!dir.exists(path)) dir.create(path)
###---- function ----###
# 检查是否需要 log2 转换
check_log_transform <- function(data) {
  
  qx <- as.numeric(quantile(data, c(0., 0.25, 0.5, 0.75, 0.99, 1.0), na.rm = TRUE))
  
  return((qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0) || (qx[2] > 0 && qx[2] < 1 && qx[4] > 1 && qx[4] < 2))
  
}

# 执行 log2 转换
log_transform <- function(data) {
  data[which(data <= 0)] <- 0
  return(log2(data + 1))
}

###---- 处理 GEO 数据 ----###
# get geo matrix
geo.set <- getGEO(data.id, destdir = path, GSEMatrix = T, getGPL = F)
expr.data <- Biobase::exprs(geo.set[[1]])

# 判断是否需要 log2 转换
if (check_log_transform(expr.data)) {
  expr.data <- log_transform(expr.data)
  print("log2 transform finished")
} else {
  print("log2 transform not needed")
}
#"log2 transform not needed"


# 获取 GPL 注释数据

gpls <- geo.set[[1]]@annotation

gpl.raw <- getGEO(gpls, destdir = file.path(getwd(), "GPL"))
gpl <- Table(gpl.raw)
colnames(gpl)

gpl <- gpl %>%
  dplyr::select(ID, `Gene Symbol`) %>%
  dplyr::rename(id = ID, symbol = `Gene Symbol`) %>%
  separate(symbol, into = c("symbol", NA), sep = " /// ")
gpl <- gpl %>% filter(symbol != "")
# 匹配并处理表达矩阵数据
data <- expr.data[rownames(expr.data) %in% gpl$id,]
gpl <- gpl[match(rownames(data), gpl$id),]

data <- as.data.frame(data) %>%
  tibble::rownames_to_column(var = "rowname") %>%
  mutate(rowname = gpl[match(rowname, gpl$id), "symbol"]) %>%
  mutate(rowMean = rowMeans(dplyr::select(., starts_with("GSM")))) %>%
  arrange(desc(rowMean)) %>%
  distinct(rowname, .keep_all = TRUE) %>%
  dplyr::select(-rowMean) %>%
  {rownames(.) <- NULL; .} %>%
  tibble::column_to_rownames("rowname")

range(data)
data[data < 0] <- 0

# 分组信息
phen.data <- pData(geo.set[[1]])
# View(phen.data)
names(phen.data)
Group <- data.frame(sample = phen.data$geo_accession,
                    group = phen.data$`affected status (disease state):ch1`)

table(Group$group)


group <- Group %>%
  na.omit() %>%
  mutate(group = ifelse(group == "control", "Control", "Disease"))

group <- group[order(group$group),]
head(group)

data <- data[, group$sample]

write.csv(data, sprintf("%s/01.data_%s.csv", path, data.id), row.names = T, quote = T)
write.csv(group, sprintf("%s/02.group_%s.csv", path, data.id), row.names = T, quote = T)






