# assets 目录说明

本目录用于存放脚本运行所需的**本地素材**（函数集与基因集，未随仓库提交）。
脚本一律使用相对路径，素材缺失时会在对应位置给出提示（部分素材可自动获取）。

## 一、目录结构约定

```
qxxncz/                        <- 项目根目录（含 00.rawdata、01.different_genes … 14.scRNA）
├── Rscript/                   <- 本套分析脚本（R.00 ~ R.14）
├── assets/                    <- 本目录：本地静态素材
│   ├── README.md
│   ├── 113MLtoolkit/          <- 113 种机器学习组合框架的本地函数集（R.05 使用）
│   ├── mycolor.R              <- 配色向量 color（R.14 可选）
│   ├── CellReports.txt        <- 基因集文件（R.14 可选）
│   └── geneset/
│       └── c2.cp.kegg.v7.5.1.symbols.gmt   <- KEGG 基因集（R.09 可选，缺失时自动获取）
└── 00.rawdata、01.different_genes … 14.scRNA   <- 输入数据与输出结果
```

脚本通过以下顺序自动定位项目根目录（即 `qxxncz/`），无需手工修改路径：

1. 环境变量 `PROJECT_ROOT`（如有设置，优先级最高）；
2. 脚本自身所在目录的上一级（脚本位于 `<项目根>/Rscript/`）；
3. 当前工作目录及其上一级中含 `00.rawdata` 的目录。

三种方式覆盖了 `Rscript Rscript/R.xx.R`、RStudio 中 Source、以及把脚本与数据放在同一目录等常见情形。
定位结果会在运行时打印为 `项目根目录：...`，可据此确认是否识别正确。

## 二、素材清单

| 路径 | 用途 | 是否必需 | 缺失时的行为 |
|---|---|---|---|
| `assets/113MLtoolkit/` | R.05 机器学习流程使用的函数集：`113ML_functions_noneed.R`、`113AUCplot_adjust.R`、`113ML_adjust.R`、`methods.xlsx` | 是 | R.05 会在开头停止并提示该目录缺失，其余脚本不受影响 |
| `assets/geneset/c2.cp.kegg.v7.5.1.symbols.gmt` | R.09 的 KEGG 通路基因集 | 可选 | 自动改用 `msigdbr`（C2:CP:KEGG）获取，并缓存到该路径；不同版本的基因集名称可能略有差异 |
| `assets/mycolor.R` | R.14 中 `color` 配色向量（`mycols <- color[51:70]`） | 可选 | 改用脚本内置备用配色（`wesanderson` 组合），并给出 warning；聚类配色与默认配色不同 |
| `assets/CellReports.txt` | R.14 中构建 `geneSetsList` | 可选 | 跳过 `geneSetsList` 构建（该对象在本脚本中未参与后续计算），不影响结果 |
| `00.rawdata/GPL/` | R.00 的 GEO 平台注释缓存目录 | 自动 | 运行时由 `getGEO()` 自动下载到 `00.rawdata/GPL/`，无需手工准备 |
| `14.scRNA/00.rawdata/GSE225948_RAW.tar` | R.14 单细胞原始数据（GEO 公开数据） | 必需（R.14） | 需下载后解压至 `14.scRNA/00.rawdata/GSE225948/` |

## 三、运行方式

在项目根目录 `qxxncz/` 下，按顺序执行：

```bash
Rscript Rscript/R.00.rawdata.R
Rscript Rscript/R.01.limma.R
Rscript Rscript/R.02.venn.R
Rscript Rscript/R.03.enrich.R
Rscript Rscript/R.04.PPI.R
Rscript Rscript/R.05.113ML.R      # 需先准备 assets/113MLtoolkit
Rscript Rscript/R.06.expressionROC.R
Rscript Rscript/R.07.nomogram.R
Rscript Rscript/R.08.location.R
Rscript Rscript/R.09.GSEA.R
Rscript Rscript/R.10.cibersort.R
Rscript Rscript/R.11.multiMiR.R
Rscript Rscript/R.13.enrichR.R
Rscript Rscript/R.14.scRNA.R      # 需先下载 14.scRNA/00.rawdata/GSE225948 原始数据
```

也可在 RStudio 中打开对应脚本后 Source。各脚本的输出写入 `<项目根>/<对应序号目录>/`，
读取输入统一使用 `../00.rawdata/...` 等相对路径，整体可随项目文件夹一起移动。
