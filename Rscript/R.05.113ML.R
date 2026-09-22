# .rs.restartR()

rm(list = ls());gc()

print(Sys.time())
library(dplyr)
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

wd <- file.path(root,'05.machine_learning')
dir.create(wd,showWarnings = F,recursive = T)
setwd(wd)



suppressPackageStartupMessages({ library(caret); library(data.table); library(e1071); library(ggplot2); library(glmnet); library(MASS); library(openxlsx); library(pROC) })
strip_utf8_bom <- function(x) {
  gsub("^\ufeff", "", as.character(x), perl = TRUE)
}

# 候选基因超过 100 个时，按 |log2FoldChange| 保留 top100
prepare_ml_input_top100 <- function(candidate_file, deg_file, outdir) {
  ensure_dir(outdir)
  input_genes <- read_gene_vector_local(candidate_file)
  input_count <- length(input_genes)
  if (input_count <= 100) {
    keep_file <- file.path(outdir, 'selected_genes_for_ml.tsv')
    data.table::fwrite(data.frame(symbol = input_genes, stringsAsFactors = FALSE), keep_file, sep = '\t')
    data.table::fwrite(data.frame(gene = input_genes, input_order = seq_along(input_genes), log2FoldChange = NA_real_, abs_log2FoldChange = NA_real_, rank = seq_along(input_genes), stringsAsFactors = FALSE), file.path(outdir, 'ranked_genes_by_abs_logfc.tsv'), sep = '\t')
    write_status_local(file.path(outdir, 'selection_summary.tsv'), 'not_required', 'Gene count <= 100; kept all genes.', list(input_gene_count = input_count, selected_gene_count = input_count, method = 'keep_original'))
    return(list(status = 'not_required', selected_file = keep_file, gene_count = input_count))
  }
  deg_df <- data.table::fread(deg_file, data.table = FALSE)
  rank_df <- merge(data.frame(gene = input_genes, input_order = seq_along(input_genes), stringsAsFactors = FALSE), deg_df[, c('gene', 'log2FoldChange'), drop = FALSE], by = 'gene', all.x = TRUE, sort = FALSE)
  rank_df$abs_log2FoldChange <- abs(suppressWarnings(as.numeric(rank_df$log2FoldChange)))
  rank_df$abs_log2FoldChange[is.na(rank_df$abs_log2FoldChange)] <- -Inf
  rank_df <- rank_df[order(-rank_df$abs_log2FoldChange, rank_df$input_order), , drop = FALSE]
  rank_df$rank <- seq_len(nrow(rank_df))
  selected_genes <- utils::head(rank_df$gene, 100)
  selected_file <- file.path(outdir, 'selected_genes_for_ml.tsv')
  data.table::fwrite(rank_df, file.path(outdir, 'ranked_genes_by_abs_logfc.tsv'), sep = '\t')
  data.table::fwrite(data.frame(symbol = selected_genes, stringsAsFactors = FALSE), selected_file, sep = '\t')
  write_status_local(file.path(outdir, 'selection_summary.tsv'), 'ok', 'Selected top100 genes ranked by abs(log2FoldChange).', list(input_gene_count = input_count, selected_gene_count = length(selected_genes), method = 'abs_log2FoldChange_top100'))
  list(status = 'ok', selected_file = selected_file, gene_count = length(selected_genes))
}


# ---- 本地素材与项目目录 ----
# 工具包目录（需包含 113ML_functions_noneed.R、113AUCplot_adjust.R、113ML_adjust.R、methods.xlsx）
assets_dir <- file.path(root, "assets", "113MLtoolkit")
if (!file.exists(file.path(assets_dir, "113ML_functions_noneed.R"))) {
  stop("未找到 ML 工具包目录：", assets_dir,
       "\n请将 113MLtoolkit 函数集（含 113ML_functions_noneed.R、113AUCplot_adjust.R、113ML_adjust.R、methods.xlsx）放入 assets/ 下，详见 assets/README.md。",
       call. = FALSE)
}
source(file.path(assets_dir, '113ML_functions_noneed.R'))
source(file.path(assets_dir, '113AUCplot_adjust.R')) # 载入作图工具函数

ORIGINAL_DIR <- root   # 项目根目录（输入数据所在处）


expr_path <- "../00.rawdata/GSE16561/01.data_GSE16561.csv"
expr_hint <- "unknown"
group_path <- "../00.rawdata/GSE16561/02.group_GSE16561.csv"

validation_expr_path <- '../00.rawdata/GSE58294/01.data_GSE58294.csv'
validation_expr_hint <- "unknown"
validation_group_path <- "../00.rawdata/GSE58294/02.group_GSE58294.csv"
case_label <- "Disease"
control_label <- "Control"

# 建模纳入的基因数范围
min_vars <- 2L
max_vars <- 5L

cvfold <- 5   # 交叉验证折数

output <- ensure_dir(file.path(ORIGINAL_DIR, '05.machine_learning'))
ML_ROOT <- output

candidate_file <- file.path(ORIGINAL_DIR, '02.venn', "01.cadi.genes.csv")
ml_filter_dir <- ensure_dir(file.path(ML_ROOT, 'ml_input_filter'))
deg_dir <- file.path(ORIGINAL_DIR, "01.different_genes", "02.DEGs_sig(GSE16561).csv")
ml_res <- prepare_ml_input_top100(candidate_file, deg_dir, ml_filter_dir)
genes_file <- ml_res$selected_file

outdir <- ensure_dir(file.path(ML_ROOT, '113ML'))
logs_dir <- ensure_dir(file.path(outdir, '00_logs'))
input_dir <- ensure_dir(file.path(outdir, '01_input'))
methods_dir <- ensure_dir(file.path(outdir, '02_methods'))
pretrain_dir <- ensure_dir(file.path(outdir, '03_pretrain'))
model_dir <- ensure_dir(file.path(outdir, '04_models'))
output_dir <- ensure_dir(file.path(outdir, '05_output'))
plots_dir <- ensure_dir(file.path(output_dir, 'plots'))
current_stage_file <- file.path(logs_dir, 'current_stage.tsv')
stage_history_file <- file.path(logs_dir, 'stage_history.tsv')
log_con <- file(file.path(logs_dir, 'workflow.log'), 'at')
log_msg('113 loading ML assets from: ', assets_dir)
source(file.path(assets_dir, '113ML_adjust.R'))
log_msg('113ML_adjust.R loaded successfully.')
log_msg("113 wrapper start")

on.exit(close(log_con), add = TRUE)
write_table_dual <- function(df, root_name, subdir_file) {
  data.table::fwrite(df, file.path(outdir, root_name), sep = '\t')
  data.table::fwrite(df, subdir_file, sep = '\t')
}
write_status_dual <- function(df) {
  write_status_table(df, file.path(outdir, 'status.tsv'))
  write_status_table(df, file.path(output_dir, 'status.tsv'))
}
write_stage_state('startup', 'bootstrap', 'running', '113 reproducible script bootstrap initialized.')
pretrain_trace <- list()
model_trace <- list()
evaluation_trace <- list()
has_assets <- dir.exists(assets_dir) && file.exists(file.path(assets_dir, '113ML_adjust.R')) && file.exists(file.path(assets_dir, 'methods.xlsx'))


# if (!has_assets) {
#   file.copy(genes_file, file.path(outdir, 'selected_features.tsv'), overwrite = TRUE)
#   write_status_local(file.path(outdir, 'status.tsv'), 'skipped_missing_assets', paste0('113 assets not found in ', assets_dir, '; fallback selected_features.tsv was created.'), list(selected_feature_count = length(read_gene_vector_local(genes_file))))
#   data.table::fwrite(data.frame(route = '113_input_bypass', description = '113 assets missing; ML input genes were retained.', feature_count = length(read_gene_vector_local(genes_file)), stringsAsFactors = FALSE), file.path(outdir, 'best_model.tsv'), sep = '\t')
#   quit(save = 'no')
# }

expr_info <- prepare_expression_for_analysis(read_expr_matrix(expr_path), expr_path = expr_path, expr_hint = expr_hint)
group_df <- read_group_file(group_path)
aligned <- align_expr_group(expr_info$expr, group_df, case_label, control_label)
expr <- aligned$expr
group_df <- aligned$group

val_expr_info <- prepare_expression_for_analysis(read_expr_matrix(validation_expr_path), expr_path = validation_expr_path, expr_hint = validation_expr_hint)
val_group <- read_group_file(validation_group_path)
val_aligned <- align_expr_group(val_expr_info$expr, val_group, case_label, control_label)
val_expr <- val_aligned$expr
val_group <- val_aligned$group

fwrite(data.frame(symbol = candidate_genes <- intersect(read_gene_vector(genes_file), intersect(rownames(expr), rownames(val_expr)))), file.path(input_dir, "candidate_genes_in_both_datasets.tsv"), sep = "\t")
fwrite(group_df, file.path(input_dir, "train_group_aligned.tsv"), sep = "\t")
fwrite(val_group, file.path(input_dir, "validation_group_aligned.tsv"), sep = "\t")
writeLines(c(
  paste0("train_samples\t", nrow(group_df)),
  paste0("validation_samples\t", nrow(val_group)),
  paste0("candidate_genes\t", length(candidate_genes)),
  paste0("train_transform_method\t", expr_info$transform_method),
  paste0("train_expr_hint\t", expr_info$expr_hint),
  paste0("validation_transform_method\t", val_expr_info$transform_method),
  paste0("validation_expr_hint\t", val_expr_info$expr_hint)
), con = file.path(input_dir, "input_summary.txt"))
if (length(candidate_genes) < min_vars) {
  log_msg("Skip 113 ML because shared candidate gene count is ", length(candidate_genes), " < ", min_vars)
  write_stage_state("input_check", "shared_genes", "failed_too_few_genes", paste0("shared_gene_count=", length(candidate_genes), "; min_vars=", min_vars))
  save_113_plots(reason = paste0("Too few shared genes for 113 ML: ", length(candidate_genes), " < ", min_vars))
  write_status_dual(data.frame(status = "skip", message = "Too few genes for 113 ML."))
  quit(save = "no")
}

train_expr <- as.data.frame(t(expr[candidate_genes, group_df$sample, drop = FALSE]))
val_expr2 <- as.data.frame(t(val_expr[candidate_genes, val_group$sample, drop = FALSE]))
train_expr <- train_expr[, intersect(colnames(train_expr), colnames(val_expr2)), drop = FALSE]
val_expr2 <- val_expr2[, colnames(train_expr), drop = FALSE]

Train_class <- data.frame(
  Sample = rownames(train_expr),
  outcome = ifelse(group_df$group[match(rownames(train_expr), group_df$sample)] == case_label, 1, 0),
  Cohort = "Train"
)
Train_class$outcome.group <- factor(Train_class$outcome, levels = c(0, 1))

Test_class <- data.frame(
  Sample = rownames(val_expr2),
  outcome = ifelse(val_group$group[match(rownames(val_expr2), val_group$sample)] == case_label, 1, 0),
  Cohort = "Test"
)
Test_class$outcome.group <- factor(Test_class$outcome, levels = c(0, 1))

Train_set <- scaleData(data = train_expr, centerFlags = FALSE, scaleFlags = FALSE)
Test_set <- scaleData(data = val_expr2, cohort = Test_class$Cohort, centerFlags = FALSE, scaleFlags = FALSE)

methods_xlsx <- file.path(assets_dir, "methods.xlsx")
catalog <- read_methods_catalog(methods_xlsx)
method_labels <- catalog$methods
methods <- compact_method_label(method_labels)
keep_methods <- !is.na(methods) & nzchar(methods)
method_labels <- method_labels[keep_methods]
methods <- methods[keep_methods]
method_map <- data.frame(
  method_label = method_labels,
  method_compact = methods,
  stringsAsFactors = FALSE
)

# 载入ML包--------
package_status <- load_method_packages(method_labels)
write_table_dual(method_map, "methods_catalog.tsv", file.path(methods_dir, "methods_catalog.tsv"))
if (!is.null(catalog$score_df)) {
  methods_score_df <- merge(method_map, catalog$score_df, by = "method_label", all.x = TRUE, sort = FALSE)
  write_table_dual(methods_score_df, "methods_catalog_scores.tsv", file.path(methods_dir, "methods_catalog_scores.tsv"))
}
write_table_dual(package_status$package_df, "methods_required_packages.tsv", file.path(methods_dir, "methods_required_packages.tsv"))

if (length(package_status$unsupported) > 0) {
  msg <- paste0("Unsupported algorithm component(s) in methods.xlsx: ", paste(package_status$unsupported, collapse = ", "))
  write_stage_state("package_check", "methods_catalog", "failed_unsupported_methods", msg)
  save_113_plots(reason = msg)
  write_status_dual(
    data.frame(
      status = "failed_unsupported_methods",
      message = msg,
      unsupported_components = paste(package_status$unsupported, collapse = ",")
    )
  )
  stop(msg, call. = FALSE)
}

if (nrow(package_status$missing_df) > 0) {
  missing_pkgs <- unique(package_status$missing_df$package)
  msg <- paste0("Missing required R package(s) for 113 ML: ", paste(missing_pkgs, collapse = ", "), ". Please install them and rerun.")
  write_stage_state("package_check", "methods_catalog", "failed_missing_packages", msg)
  save_113_plots(reason = msg)
  write_status_dual(
    data.frame(
      status = "failed_missing_packages",
      message = msg,
      missing_packages = paste(missing_pkgs, collapse = ",")
    )
  )
  stop(msg, call. = FALSE)
}

methods_ready_df <- method_map
write_table_dual(methods_ready_df, "methods_ready.tsv", file.path(methods_dir, "methods_ready.tsv"))
methods <- methods_ready_df$method_compact
write_stage_state("methods_ready", "methods_catalog", "ok", paste0("ready_method_count=", length(methods)))
log_msg("Method catalog sheet: ", catalog$sheet)
log_msg("Loaded method packages: ", ifelse(length(package_status$loaded) > 0, paste(package_status$loaded, collapse = ", "), "none"))
log_msg("Method catalog size: ", nrow(method_map), "; ready methods: ", length(methods))
log_msg("Method preview: ", paste(head(methods_ready_df$method_label, 20), collapse = ", "))


# Pretrain：保留 5 种 Pretrain 方式
Variable <- colnames(Train_set)
preTrain.method <- strsplit(methods, "\\+")
preTrain.method <- lapply(preTrain.method, function(x) rev(x)[-1])
preTrain.method <- unique(unlist(preTrain.method))
preTrain.method
package_status$components
# preTrain.method <- preTrain.method[preTrain.method %in% package_status$components]
preTrain.method <- preTrain.method[gsub("\\[.*\\]", "", preTrain.method) %in% package_status$components]
preTrain.method

preTrain.var <- list()
set.seed(777)
for (method in preTrain.method) {
  set.seed(777)
  pretrain_started_at <- Sys.time()
  write_stage_state("pretrain", method, "running", "Pretrain variable selection started.")
  log_msg("Pretrain variable selection started: ", method)
  pretrain_error <- NULL
  preTrain.var[[method]] <- tryCatch(
    RunML(method = method, Train_set = Train_set, Train_label = Train_class, mode = "Variable", classVar = "outcome.group"),
    error = function(e) {
      pretrain_error <<- conditionMessage(e)
      character()
    }
  )
  if (is.null(preTrain.var[[method]])) preTrain.var[[method]] <- character()
  status_now <- if (length(preTrain.var[[method]]) > 0) "ok" else if (!is.null(pretrain_error)) "failed" else "empty"
  pretrain_msg <- if (!is.null(pretrain_error) && nzchar(pretrain_error)) pretrain_error else if (identical(status_now, "empty")) "No variables returned by pretrain step." else NA_character_
  write_stage_state("pretrain", method, status_now, pretrain_msg, as.numeric(difftime(Sys.time(), pretrain_started_at, units = "secs")))
  record_pretrain_trace(method, status_now, length(preTrain.var[[method]]), pretrain_msg)
  log_msg("Pretrain variable selection finished: ", method, "; status=", status_now, "; variables=", length(preTrain.var[[method]]), if (!is.na(pretrain_msg)) paste0("; message=", pretrain_msg) else "")
}
preTrain.var[["simple"]] <- colnames(Train_set)
write_stage_state("pretrain", "simple", "ok", paste0("variables=", length(preTrain.var[["simple"]])), 0)
record_pretrain_trace("simple", "ok", length(preTrain.var[["simple"]]), "All candidate genes without pretrain filtering.")

pretrain_summary <- rbindlist(lapply(names(preTrain.var), function(method) {
  vars <- preTrain.var[[method]]
  var_file <- file.path(pretrain_dir, paste0("variables_", method, ".txt"))
  writeLines(as.character(vars), con = var_file)
  data.frame(method = method, variable_count = length(vars), variable_file = var_file, stringsAsFactors = FALSE)
}), fill = TRUE)
fwrite(pretrain_summary, file.path(pretrain_dir, "pretrain_variable_summary.tsv"), sep = "\t")
if (length(pretrain_trace) > 0) {
  fwrite(rbindlist(pretrain_trace, fill = TRUE), file.path(pretrain_dir, "pretrain_trace.tsv"), sep = "\t")
}

# Train--------
model <- list()
model_records <- list()
set.seed(777)
Train_set_bk <- Train_set
for (method_name in methods) {
  set.seed(777)
  parts <- strsplit(method_name, "\\+")[[1]]
  if (length(parts) == 1) parts <- c("simple", parts)
  variable_set <- preTrain.var[[parts[1]]]
  method_label_now <- method_map$method_label[match(method_name, method_map$method_compact)]
  if (length(variable_set) == 0) {
    log_msg("113 model skipped: ", method_name, "; reason=no variables from pretrain=", parts[1])
    write_stage_state("model", method_name, "skipped_no_variables", paste0("pretrain=", parts[1], "; method_label=", method_label_now))
    record_model_trace(method_name, method_label_now, parts[1], parts[2], "skipped_no_variables", 0L, paste0("No variables from pretrain method ", parts[1]), NA_character_)
    next
  }
  train_sub <- Train_set_bk[, variable_set, drop = FALSE]
  model_started_at <- Sys.time()
  write_stage_state("model", method_name, "running", paste0("pretrain=", parts[1], "; fit_method=", parts[2], "; train_vars=", ncol(train_sub), "; method_label=", method_label_now))
  log_msg("113 model training started: ", method_name, " using pretrain=", parts[1], ", train_vars=", ncol(train_sub))
  fit_error <- NULL
  fit <- tryCatch(
    RunML(method = parts[2], Train_set = train_sub, Train_label = Train_class, mode = "Model", classVar = "outcome.group"),
    error = function(e) {
      fit_error <<- conditionMessage(e)
      NULL
    }
  )
  if (is.null(fit)) {
    log_msg("113 model training failed: ", method_name, "; message=", fit_error)
    write_stage_state("model", method_name, "failed", fit_error, as.numeric(difftime(Sys.time(), model_started_at, units = "secs")))
    record_model_trace(method_name, method_label_now, parts[1], parts[2], "failed", NA_integer_, fit_error, NA_character_)
    next
  }
  n_vars <- length(ExtractVar(fit))
  selected_vars <- ExtractVar(fit)
  if (n_vars < min_vars || n_vars > max_vars) {
    log_msg("113 model filtered by feature-count rule: ", method_name, "; selected_feature_count=", n_vars, "; allowed_range=", min_vars, "-", max_vars, "; features=", paste(selected_vars, collapse = ";"))
    write_stage_state("model", method_name, "filtered_feature_count", paste0("selected_feature_count=", n_vars, "; allowed_range=", min_vars, "-", max_vars), as.numeric(difftime(Sys.time(), model_started_at, units = "secs")))
    record_model_trace(method_name, method_label_now, parts[1], parts[2], "filtered_feature_count", n_vars, paste0("Allowed range: ", min_vars, "-", max_vars), paste(selected_vars, collapse = ";"))
    next
  }
  model[[method_name]] <- fit
  model_records[[length(model_records) + 1L]] <- data.frame(
    method = method_name,
    method_label = method_label_now,
    pretrain_method = parts[1],
    fit_method = parts[2],
    selected_feature_count = n_vars,
    selected_features = paste(selected_vars, collapse = ";"),
    stringsAsFactors = FALSE
  )
  write_stage_state("model", method_name, "ok", paste0("selected_feature_count=", n_vars, "; method_label=", method_label_now), as.numeric(difftime(Sys.time(), model_started_at, units = "secs")))
  record_model_trace(method_name, method_label_now, parts[1], parts[2], "ok", n_vars, NA_character_, paste(selected_vars, collapse = ";"))
  log_msg("113 model training finished: ", method_name, "; selected_feature_count=", n_vars, "; features=", paste(selected_vars, collapse = ";"))
}
Train_set <- Train_set_bk

if (length(model_records) > 0) {
  model_record_df <- rbindlist(model_records, fill = TRUE)
  fwrite(model_record_df, file.path(model_dir, "model_feature_summary.tsv"), sep = "\t")
  feature_long_df <- rbindlist(lapply(seq_len(nrow(model_record_df)), function(i) {
    feats <- unlist(strsplit(as.character(model_record_df$selected_features[i]), ";", fixed = TRUE), use.names = FALSE)
    feats <- feats[nzchar(feats)]
    if (length(feats) == 0) return(NULL)
    data.frame(
      method = model_record_df$method[i],
      method_label = model_record_df$method_label[i],
      feature = feats,
      stringsAsFactors = FALSE
    )
  }), fill = TRUE)
  if (!is.null(feature_long_df) && nrow(feature_long_df) > 0) {
    fwrite(feature_long_df, file.path(model_dir, "feature_catalog.tsv"), sep = "\t")
  }
} else {
  model_record_df <- data.frame()
}
model_trace_df <- if (length(model_trace) > 0) rbindlist(model_trace, fill = TRUE) else NULL
if (!is.null(model_trace_df)) {
  fwrite(model_trace_df, file.path(model_dir, "model_training_trace.tsv"), sep = "\t")
}

if (length(model) == 0) {
  log_msg("No usable 113 models after variable-count filtering")
  write_stage_state("model", "all_methods", "failed_no_usable_models", "No usable 113 models after variable-count filtering.")
  save_113_plots(reason = "No usable 113 models after variable-count filtering.", model_trace_df = model_trace_df)
  write_status_dual(data.frame(status = "failed", message = "No usable 113 models after variable-count filtering."))
  quit(save = "no", status = 1L)
}

auc_list <- list()
for (method_name in names(model)) {
  evaluation_started_at <- Sys.time()
  write_stage_state("evaluation", method_name, "running", "AUC evaluation started.")
  log_msg("113 evaluation started: ", method_name)
  eval_error <- NULL
  auc_list[[method_name]] <- tryCatch(
    RunEval(
      fit = model[[method_name]],
      Test_set = Test_set,
      Test_label = Test_class,
      Train_set = Train_set,
      Train_label = Train_class,
      Train_name = "Train",
      cohortVar = "Cohort",
      classVar = "outcome.group"
    ),
    error = function(e) {
      eval_error <<- conditionMessage(e)
      c(Train = NA_real_, Test = NA_real_)
    }
  )
  eval_now <- auc_list[[method_name]]
  avg_now <- mean(as.numeric(eval_now[c("Train", "Test")]), na.rm = TRUE)
  if (!is.finite(avg_now)) avg_now <- NA_real_
  status_now <- if (all(is.na(eval_now[c("Train", "Test")]))) "failed" else "ok"
  write_stage_state("evaluation", method_name, status_now, eval_error, as.numeric(difftime(Sys.time(), evaluation_started_at, units = "secs")))
  record_evaluation_trace(method_name, method_map$method_label[match(method_name, method_map$method_compact)], status_now, as.numeric(eval_now["Train"]), as.numeric(eval_now["Test"]), avg_now, eval_error)
  log_msg("113 evaluation finished: ", method_name, "; status=", status_now, "; train_auc=", sprintf("%.4f", as.numeric(eval_now["Train"])), "; validation_auc=", sprintf("%.4f", as.numeric(eval_now["Test"])), "; avg_auc=", sprintf("%.4f", avg_now), if (!is.null(eval_error)) paste0("; message=", eval_error) else "")
}
evaluation_trace_df <- if (length(evaluation_trace) > 0) rbindlist(evaluation_trace, fill = TRUE) else NULL
if (!is.null(evaluation_trace_df)) {
  fwrite(evaluation_trace_df, file.path(model_dir, "evaluation_trace.tsv"), sep = "\t")
}

auc_mat <- do.call(rbind, auc_list)
auc_df <- as.data.frame(auc_mat)
auc_df$method <- rownames(auc_df)
auc_df$method_label <- method_map$method_label[match(auc_df$method, method_map$method_compact)]
auc_df$avg_auc <- rowMeans(auc_df[, c("Train", "Test"), drop = FALSE], na.rm = TRUE)
auc_df$evaluation_ok <- !is.na(auc_df$Train) & !is.na(auc_df$Test)
auc_df <- auc_df %>% dplyr::filter(Test!=1)
log_msg("113 evaluation summary: total_models=", nrow(auc_df), "; evaluable_models=", sum(auc_df$evaluation_ok))
if (!any(auc_df$evaluation_ok)) {
  log_msg("No evaluable 113 models after AUC calculation")
  write_stage_state("evaluation", "all_methods", "failed_no_evaluable_models", "No evaluable 113 models after AUC calculation.")
  save_113_plots(reason = "No evaluable 113 models after AUC calculation.", model_trace_df = model_trace_df, evaluation_trace_df = evaluation_trace_df)
  write_status_dual(data.frame(status = "failed", message = "No evaluable 113 models after AUC calculation."))
  quit(save = "no", status = 1L)
}
auc_df <- auc_df[auc_df$evaluation_ok, , drop = FALSE]
auc_df$meets_auc_threshold <- !is.na(auc_df$Train) & !is.na(auc_df$Test) & auc_df$Train > 0.7 & auc_df$Test > 0.7
auc_df$selection_priority <- ifelse(auc_df$meets_auc_threshold, 1L, 0L)
auc_df <- auc_df[order(-auc_df$selection_priority, -auc_df$avg_auc, -auc_df$Test, -auc_df$Train), , drop = FALSE]
qualified_df <- auc_df[auc_df$meets_auc_threshold, , drop = FALSE]
selection_pool <- if (nrow(qualified_df) > 0) "qualified_auc_pool" else "plot_only_no_qualified_model"
plot_rank_df <- auc_df[order(-auc_df$avg_auc, -auc_df$Test, -auc_df$Train), , drop = FALSE]
plot_model_method <- plot_rank_df$method[1]
plot_model_label <- method_label_or_name(plot_rank_df$method[1], plot_rank_df$method_label[1])
has_qualified_best <- nrow(qualified_df) > 0
if (has_qualified_best) {
  best_rank_df <- qualified_df[order(-qualified_df$avg_auc, -qualified_df$Test, -qualified_df$Train), , drop = FALSE]
  best_rank_df <- best_rank_df %>% dplyr::filter(Test!=1)
  best_method <- best_rank_df$method[1]
  best_method_label <- method_label_or_name(best_rank_df$method[1], best_rank_df$method_label[1])
  best_genes <- ExtractVar(model[[best_method]])
} else {
  best_method <- NA_character_
  best_method_label <- NA_character_
  best_genes <- character()
}
auc_matrix_df <- auc_df[, c("method", "method_label", "Train", "Test", "avg_auc", "selection_priority", "meets_auc_threshold"), drop = FALSE]
colnames(auc_matrix_df)[colnames(auc_matrix_df) == "Test"] <- "Validation"
fwrite(auc_matrix_df, file.path(output_dir, "AUC_mat.tsv"), sep = "\t")
method_detail_root <- ensure_dir(file.path(model_dir, "methods"))
for (i in seq_len(nrow(auc_df))) {
  method_name <- auc_df$method[i]
  method_subdir <- ensure_dir(file.path(method_detail_root, safe_component(method_name)))
  features_i <- ExtractVar(model[[method_name]])
  meta_i <- if (nrow(model_record_df) > 0) model_record_df[model_record_df$method == method_name, , drop = FALSE] else data.frame()
  summary_i <- data.frame(
    method = method_name,
    method_label = auc_df$method_label[i],
    train_auc = auc_df$Train[i],
    test_auc = auc_df$Test[i],
    avg_auc = auc_df$avg_auc[i],
    meets_auc_threshold = auc_df$meets_auc_threshold[i],
    selection_priority = auc_df$selection_priority[i],
    selected_best = has_qualified_best && identical(method_name, best_method),
    selected_for_plot = identical(method_name, plot_model_method),
    feature_count = length(features_i),
    pretrain_method = if (nrow(meta_i) > 0) meta_i$pretrain_method[1] else NA_character_,
    fit_method = if (nrow(meta_i) > 0) meta_i$fit_method[1] else NA_character_,
    stringsAsFactors = FALSE
  )
  fwrite(summary_i, file.path(method_subdir, "method_summary.tsv"), sep = "\t")
  fwrite(data.frame(symbol = features_i, stringsAsFactors = FALSE), file.path(method_subdir, "selected_features.tsv"), sep = "\t")
}
plot_fit <- model[[plot_model_method]]
plot_roc <- build_best_model_predictions(plot_fit, Train_set, Train_class, Test_set, Test_class, plot_model_label)
fwrite(plot_roc$prediction_df, file.path(output_dir, "best_model_prediction_scores.tsv"), sep = "\t")
fwrite(plot_roc$summary_df, file.path(output_dir, "best_model_roc_summary.tsv"), sep = "\t")
fwrite(plot_roc$curve_df, file.path(output_dir, "best_model_roc_curve_points.tsv"), sep = "\t")
plot_model_df <- data.frame(
  method = plot_model_method,
  method_label = plot_model_label,
  train_auc = plot_rank_df$Train[1],
  test_auc = plot_rank_df$Test[1],
  avg_auc = plot_rank_df$avg_auc[1],
  meets_auc_threshold = plot_rank_df$meets_auc_threshold[1],
  selection_pool = selection_pool,
  plot_only = !has_qualified_best,
  feature_count = length(ExtractVar(plot_fit)),
  stringsAsFactors = FALSE
)
fwrite(plot_model_df, file.path(output_dir, "plot_model.tsv"), sep = "\t")
log_msg("113 model used for plotting: method=", plot_model_label, "; avg_auc=", sprintf("%.4f", plot_rank_df$avg_auc[1]), "; train_auc=", sprintf("%.4f", plot_rank_df$Train[1]), "; test_auc=", sprintf("%.4f", plot_rank_df$Test[1]), "; meets_auc_threshold=", plot_rank_df$meets_auc_threshold[1])
log_msg("113 top-ranked methods: ", paste(utils::head(paste0(auc_df$method_label, "[Train=", sprintf("%.3f", auc_df$Train), ",Validation=", sprintf("%.3f", auc_df$Test), ",Avg=", sprintf("%.3f", auc_df$avg_auc), "]"), 8), collapse = " | "))

write_table_dual(auc_df, "auc_summary.tsv", file.path(model_dir, "auc_summary.tsv"))

save_113_plots(
  auc_df = auc_df,
  reason = if (has_qualified_best) "" else paste0("No model met Train>0.7 and Validation>0.7. Plotting top average-AUC model only: ", plot_model_label),
  model_trace_df = model_trace_df,
  evaluation_trace_df = evaluation_trace_df,
  best_roc = plot_roc,
  best_method_label = if (has_qualified_best) best_method_label else paste0(plot_model_label, " [plot-only]")
)
saveRDS(model, file.path(model_dir, "model_list.rds"))
saveRDS(model, file.path(outdir, "model_list.rds"))
if (has_qualified_best) {
  write_stage_state("finished", best_method_label, "ok", paste0("selection_pool=", selection_pool, "; selected_feature_count=", length(best_genes), "; train_auc=", sprintf("%.4f", best_rank_df$Train[1]), "; validation_auc=", sprintf("%.4f", best_rank_df$Test[1])))
  best_fit <- model[[best_method]]
  best_model_df <- data.frame(
    method = best_method,
    method_label = best_method_label,
    train_auc = best_rank_df$Train[1],
    test_auc = best_rank_df$Test[1],
    avg_auc = best_rank_df$avg_auc[1],
    meets_auc_threshold = best_rank_df$meets_auc_threshold[1],
    selection_pool = selection_pool,
    feature_count = length(best_genes),
    stringsAsFactors = FALSE
  )
  write_table_dual(best_model_df, "best_model.tsv", file.path(model_dir, "best_model.tsv"))
  write_table_dual(data.frame(symbol = best_genes), "selected_features.tsv", file.path(model_dir, "selected_features.tsv"))
  write.table(data.frame(symbol = best_genes, stringsAsFactors = FALSE), file.path(output_dir, "best_Genes.txt"), sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
  write_status_dual(
    data.frame(
      status = "ok",
      catalog_sheet = catalog$sheet,
      catalog_method_count = nrow(method_map),
      ready_method_count = nrow(methods_ready_df),
      qualified_model_count = nrow(qualified_df),
      auc_threshold = "Train>0.7 & Validation>0.7",
      best_method = best_method,
      best_method_label = best_method_label,
      selection_pool = selection_pool,
      selected_feature_count = length(best_genes)
    )
  )
} else {
  unlink(file.path(outdir, "best_model.tsv"), force = TRUE)
  unlink(file.path(model_dir, "best_model.tsv"), force = TRUE)
  unlink(file.path(outdir, "selected_features.tsv"), force = TRUE)
  unlink(file.path(model_dir, "selected_features.tsv"), force = TRUE)
  unlink(file.path(output_dir, "best_Genes.txt"), force = TRUE)
  write_stage_state("finished", plot_model_label, "no_qualified_model", paste0("selection_pool=", selection_pool, "; plot_method=", plot_model_label, "; train_auc=", sprintf("%.4f", plot_rank_df$Train[1]), "; validation_auc=", sprintf("%.4f", plot_rank_df$Test[1])))
  write_status_dual(
    data.frame(
      status = "no_qualified_model",
      message = "No 113 model met Train>0.7 and Validation>0.7; top average-AUC model was used for plotting only.",
      catalog_sheet = catalog$sheet,
      catalog_method_count = nrow(method_map),
      ready_method_count = nrow(methods_ready_df),
      qualified_model_count = 0L,
      auc_threshold = "Train>0.7 & Validation>0.7",
      plot_method = plot_model_method,
      plot_method_label = plot_model_label,
      plot_train_auc = plot_rank_df$Train[1],
      plot_test_auc = plot_rank_df$Test[1],
      plot_avg_auc = plot_rank_df$avg_auc[1],
      selection_pool = selection_pool,
      stringsAsFactors = FALSE
    )
  )
}
top_rank_lines <- paste(
  utils::head(auc_df$method_label, 5),
  sprintf("Train=%.4f", auc_df$Train[seq_len(min(5, nrow(auc_df)))]),
  sprintf("Validation=%.4f", auc_df$Test[seq_len(min(5, nrow(auc_df)))]),
  sprintf("Avg=%.4f", auc_df$avg_auc[seq_len(min(5, nrow(auc_df)))]),
  sep = "\t"
)
writeLines(
  c(
    paste0("selection_pool\t", selection_pool),
    paste0("auc_threshold\tTrain>0.7 & Validation>0.7"),
    paste0("qualified_model_count\t", nrow(qualified_df)),
    paste0("plot_method\t", plot_model_method),
    paste0("plot_method_label\t", plot_model_label),
    paste0("plot_train_auc\t", plot_rank_df$Train[1]),
    paste0("plot_test_auc\t", plot_rank_df$Test[1]),
    paste0("plot_avg_auc\t", plot_rank_df$avg_auc[1]),
    if (has_qualified_best) paste0("best_method\t", best_method) else "best_method\tNA",
    if (has_qualified_best) paste0("best_method_label\t", best_method_label) else "best_method_label\tNA",
    if (has_qualified_best) paste0("train_auc\t", best_rank_df$Train[1]) else "train_auc\tNA",
    if (has_qualified_best) paste0("test_auc\t", best_rank_df$Test[1]) else "test_auc\tNA",
    if (has_qualified_best) paste0("avg_auc\t", best_rank_df$avg_auc[1]) else "avg_auc\tNA",
    if (has_qualified_best) paste0("selected_feature_count\t", length(best_genes)) else "selected_feature_count\t0",
    if (has_qualified_best) paste0("selected_features\t", paste(best_genes, collapse = ";")) else "selected_features\t",
    paste0("best_train_auc\t", ifelse(any(plot_roc$summary_df$dataset == "Train"), plot_roc$summary_df$auc[plot_roc$summary_df$dataset == "Train"][1], NA_real_)),
    paste0("best_validation_auc\t", ifelse(any(plot_roc$summary_df$dataset == "Validation"), plot_roc$summary_df$auc[plot_roc$summary_df$dataset == "Validation"][1], NA_real_)),
    "top_ranked_methods",
    top_rank_lines
  ),
  con = file.path(output_dir, "result_summary.txt")
)

df <- read.table('./113ML/selected_features.tsv', header = TRUE, sep = '\t')
write.csv(df, 'Best_model_gene.csv')

log_msg("113 wrapper finished; plot_method=", plot_model_label, if (has_qualified_best) paste0(", best_method=", best_method_label, ", selected_feature_count=", length(best_genes)) else ", no qualified best model", ", best_train_auc=", sprintf("%.4f", plot_roc$summary_df$auc[plot_roc$summary_df$dataset == "Train"][1]), ", best_validation_auc=", sprintf("%.4f", plot_roc$summary_df$auc[plot_roc$summary_df$dataset == "Validation"][1]), ", pretrain_trace=", file.path(pretrain_dir, "pretrain_trace.tsv"), ", model_trace=", file.path(model_dir, "model_training_trace.tsv"), ", evaluation_trace=", file.path(model_dir, "evaluation_trace.tsv"))
message('Step 4 finished: ', outdir)


save.image("113ML.rda")
load("113ML.rda")
