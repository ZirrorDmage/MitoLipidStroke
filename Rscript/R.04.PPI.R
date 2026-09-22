

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
wd <- file.path(root,'04.PPI')
dir.create(wd,showWarnings = F,recursive = T)
setwd(wd)

suppressPackageStartupMessages(library(clusterProfiler))
suppressPackageStartupMessages(library(enrichplot))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(GOplot))
suppressPackageStartupMessages(library(cowplot))
suppressPackageStartupMessages(library(patchwork))
suppressPackageStartupMessages(library(tidyverse))
library(export)
library(circlize)
library(grid)
library(graphics)
library(ComplexHeatmap)
select=dplyr::select

##function
get.string <- function(ids, score = 400, species = 9606){
  string_api_url <- "https://string-db.org/api"
  # Mapping identifiers------
  # 9606: hunman
  # 10090: mouse
  # 10116: rat
  output_format <- "tsv"
  method <- paste0("get_string_ids?species=", species, "&identifiers=")
  identifiers <- paste0(ids, collapse = '%0d')
  url <- paste0(c(string_api_url, output_format, method), collapse = '/')
  url <- paste0(url, identifiers)
  tmp <- tempfile()
  curl::curl_download(url, tmp)
  identifiers <- read.delim2(tmp) %>% pull(preferredName) %>% paste0(collapse = '%0d')
  unlink(tmp)
  # Getting the STRING network interactions-----
  # 150：（low confidence）0.15
  # 400：（medium confidence, default）0.4
  # 700：（high confidence）0.7
  # 900：（highest confidence）0.9
  output_format <- "tsv"
  method <- paste0("network?species=", species, "&required_score=", score, "&identifiers=")
  url <- paste0(c(string_api_url, output_format, method),collapse = '/')
  url <- paste0(url, identifiers)
  curl::curl_download(url, 'string_interactions.tsv')
  # Getting STRING network image-------
  output_format <- "highres_image" #  "image","svg"
  method <- paste0("network?species=", species, "&required_score=", score, "&identifiers=") 
  url <- paste0(c(string_api_url, output_format, method), collapse = '/')
  url <- paste0(url,identifiers)
  # curl::curl_download(url,'string_hires_image.png')
  print("finished")
  
}

p_fun <- function(string_data,select_dregee=1){
  # 
  filtered_data <- string_data  # 
  has_changes <- TRUE  # 
  #
  while(has_changes) {
    # 
    all_nodes <- c(filtered_data$preferredName_A, filtered_data$preferredName_B)
    # degree
    node_degrees <- as.data.frame(table(all_nodes))
    colnames(node_degrees) <- c("node", "degree")
    # degree≥1
    nodes_to_keep <- node_degrees$node[node_degrees$degree >= select_dregee]
    # filter
    new_filtered_data <- filtered_data %>%
      filter(preferredName_A %in% nodes_to_keep & 
               preferredName_B %in% nodes_to_keep)
    
    # check
    if(nrow(new_filtered_data) == nrow(filtered_data)) {
      has_changes <- FALSE  # 
    } else {
      filtered_data <- new_filtered_data  # 
    }
  }
  # print result
  print(paste("原始数据行数:", nrow(string_data)))
  print(paste("过滤后数据行数:", nrow(filtered_data)))
  # 
  final_nodes <- c(filtered_data$preferredName_A, filtered_data$preferredName_B)
  final_degrees <- as.data.frame(table(final_nodes))
  colnames(final_degrees) <- c("node", "degree")
  # 
  low_degree_nodes <- final_degrees$node[final_degrees$degree < select_dregee]
  if(length(low_degree_nodes) > 0) {
    print(paste0("警告：仍有度数<",select_dregee,"的节点:"))
    print(low_degree_nodes)
  } else {
    print(paste0("成功：所有节点的度数都≥",select_dregee))
  }
  print(paste0('包含节点：', length(unique(c(filtered_data$preferredName_A,filtered_data$preferredName_B)))))
  print(paste0('包含边：', nrow(filtered_data)))
  library(ComplexHeatmap)
  library(circlize)
  chord_data <- data.frame(
    from = filtered_data$preferredName_A,
    to = filtered_data$preferredName_B,
    value = filtered_data$score
  )
  all_genes <- unique(c(chord_data$from, chord_data$to))
  
  gene_connectivity <- table(c(chord_data$from, chord_data$to))
  gene_connectivity <- gene_connectivity[all_genes]  # order
  
  sorted_genes <- names(sort(gene_connectivity, decreasing = TRUE))
  
  #ComplexHeatmap colorRamp2()
  connectivity_range <- range(gene_connectivity)
  connectivity_col_fun <- colorRamp2(
    breaks = c(connectivity_range[1], mean(connectivity_range), connectivity_range[2]),
    colors = c("#FCDF90","#F8984E", "#F06A21")
  )
  sector_colors <- connectivity_col_fun(as.numeric(gene_connectivity))
  names(sector_colors) <- all_genes
  
  # score
  score_range <- range(chord_data$value)
  score_col_fun <- colorRamp2(
    breaks = c(score_range[1], mean(score_range), score_range[2]),
    colors = c("#A2D39A", "#0fb227", "#11b88f")
  )
  
  par(mar = c(1, 1, 3, 1))  # 
  #clear
  circos.clear()
  #  interval 
  n_sectors <- length(all_genes)
  gap_degree <- min(2, 360/(n_sectors * 3))  # 
  # 
  circos.par(
    gap.degree = gap_degree,
    start.degree = 90,  #
    # 
    canvas.xlim = c(-1, 1),
    canvas.ylim = c(-1.2, 1.2)  # 
  )
  # 
  chordDiagram(
    chord_data,
    grid.col = sector_colors,  # 
    col = score_col_fun,  # 
    transparency = 0.5,
    link.lwd = 1.5,
    link.lty = 1,
    link.sort = TRUE,
    link.decreasing = TRUE,
    annotationTrack = "grid",  #
    order = sorted_genes  #
  )
  # 
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      sector.name <- get.cell.meta.data("sector.index")
      circos.text(
        x = mean(c(get.cell.meta.data("xlim"))),
        y = get.cell.meta.data("ylim")[1] + 1.2,  
        labels = sector.name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        # cex = 0.3
        cex = 0.7
      )
    },
    bg.border = NA
  )
  title("", line = -2)
  title(
    "Protein-Protein Interaction Network",
    line = 0.3,  # 
    cex.main = 1.2,  # 
    font.main = 2    #
  )
  
  connectivity_legend <- Legend(
    title = "Degree",
    col_fun = connectivity_col_fun,
    at = round(seq(connectivity_range[1], connectivity_range[2], length = 5)),
    labels = round(seq(connectivity_range[1], connectivity_range[2], length = 5)),
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  
  score_legend <- Legend(
    title = "Interaction Score",
    col_fun = score_col_fun,
    at = round(seq(score_range[1], score_range[2], length = 5), 3),
    labels = round(seq(score_range[1], score_range[2], length = 5), 3),
    direction = "horizontal",
    title_position = "topcenter",
    title_gp = gpar(fontsize = 10, fontface = "bold"),
    labels_gp = gpar(fontsize = 8)
  )
  # two legends
  legend_list <- packLegend(connectivity_legend, score_legend, direction = "vertical")
  
  # legend
  pushViewport(viewport(x = 0.85, y = 0.15, width = 0.3, height = 0.3))
  draw(legend_list)
  popViewport()
  
  # resetting circlize
  circos.clear()
}

########PPI.net
gene <- read.csv('../02.venn/01.cadi.genes.csv') %>% pull(1)

# ppi network
options(timeout = 9999999)
get.string(gene, score = 400, species = 9606)

string_data <- read_tsv('string_interactions.tsv') %>% dplyr::select(c("preferredName_A", "preferredName_B", "score"))


library(Cairo)
CairoPDF("04.ppi.pdf", width=8, height=7)
p_fun(string_data)
dev.off()

CairoPNG("04.ppi.png", width=8, height=7, units="in", res=600)
p_fun(string_data)
dev.off()


# [1] "原始数据行数: 12"
# [1] "过滤后数据行数: 12"
# [1] "成功：所有节点的度数都≥1"
# [1] "包含节点：9"
# [1] "包含边：12"

