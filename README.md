# MitoLipidStroke

Analysis code, result tables and figures for a study on **mitochondrial- and lipid-metabolism-related
genes in ischemic stroke**, combining bulk transcriptome discovery and validation with an
exploratory mouse single-cell transcriptome analysis.

This repository is intended to accompany the manuscript. It contains the full analysis pipeline
(`Rscript/R.00` – `Rscript/R.14`), the derived result tables and figures, and the input data used
for the analysis.

---

## 1. Repository layout

| Path | Content |
|---|---|
| `Rscript/` | Analysis pipeline, 14 standalone R scripts (`R.00.rawdata.R` … `R.14.scRNA.R`) |
| `assets/` | Local assets required by some scripts (function sets, gene sets) and `assets/README.md` with the path/asset conventions |
| `00.rawdata/` | Input data: GEO series matrices and processed expression matrices, group files, and the mitochondrial/lipid-metabolism gene list |
| `01.different_genes/` … `13.drugs/` | Result tables and figures of the bulk-transcriptome modules (differential expression, Venn, enrichment, PPI, machine learning, ROC, nomogram, subcellular localization, GSEA, immune deconvolution, ceRNA/TF networks, drug prediction) |
| `14.scRNA/` | Result tables and figures of the single-cell module (QC, integration, dimensionality reduction, annotation, key cells, cell–cell communication, metabolism, pseudotime, pathway activity) |

All scripts use **relative paths only** and locate the repository root automatically, so the whole
folder can be moved or cloned anywhere. See `assets/README.md` for details.

---

## 2. Data sources and provenance

### Public datasets re-used in this study

| Accession | Title as submitted to GEO | Platform | Samples used here | Species | Original study |
|---|---|---|---|---|---|
| [GSE16561](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16561) | Gene expression analysis of peripheral whole blood RNA following ischemic stroke | GPL6883 | 63 (24 control / 39 stroke), discovery cohort | *Homo sapiens* | Barr TL et al.; PMID [20837969](https://pubmed.ncbi.nlm.nih.gov/20837969/) (related: PMID 28446746, 29263821) |
| [GSE58294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58294) | Gene Expression Following Cardioembolic Stroke | GPL570 | 92 (23 control / 69 stroke), validation cohort | *Homo sapiens* | Stamova B et al.; PMID [25036109](https://pubmed.ncbi.nlm.nih.gov/25036109/) |
| [GSE225948](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225948) | Brain and blood single-cell transcriptomic analysis in acute and subacute phases after experimental stroke | GPL19057 | 29 samples, single-cell module | *Mus musculus* | PMID [38177281](https://pubmed.ncbi.nlm.nih.gov/38177281/) |

The single-cell dataset (GSE225948) is a **mouse** experimental-stroke dataset; human hub genes were
mapped to mouse orthologue symbols for the single-cell modules.

### Gene sets compiled from published studies

`00.rawdata/线粒体和脂质代谢相关基因.xlsx` contains two gene sets assembled from published work
(sheet names carry the PubMed IDs):

| Sheet | Source | Gene set |
|---|---|---|
| `PMID 39827204…` | Sci Rep 2025; PMID [39827204](https://pubmed.ncbi.nlm.nih.gov/39827204/) | lipid-metabolism-related genes |
| `线粒体PMID36915111` | J Transl Med 2023; PMID [36915111](https://pubmed.ncbi.nlm.nih.gov/36915111/) | mitochondrial-related genes |

Their intersection (253 genes, `02.venn/MLM related genes list.csv`) was used as the
mitochondrial- and lipid-metabolism-related gene list (MLM); 11 of these genes were differentially
expressed in the discovery cohort (`02.venn/01.cadi.genes.csv`).

### Third-party data notice

The files under `00.rawdata/` are **not** covered by the licenses granted for this repository. They
are redistributed here only for convenience of reproduction and remain subject to the terms of the
original data providers (NCBI GEO) and the copyright of the original submitters. When re-using
them, please cite the original studies listed above.

---

## 3. Requirements and how to run

* **R** 4.3.3 (version under which the deposited results were produced)
* R packages installed from CRAN / Bioconductor, including: `limma`, `GEOquery`, `Biobase`,
  `tinyarray`, `clusterProfiler`, `org.Hs.eg.db`, `DOSE`, `enrichplot`, `msigdbr`, `biomaRt`,
  `GSEABase`, `ReactomeGSA`, `pROC`, `rms`, `regplot`, `rmda`, `caret`, `glmnet`, `e1071`,
  `ComplexHeatmap`, `circlize`, `igraph`, `corrplot`, the `tidyverse` family, `ggpubr`, `ggvenn`,
  `VennDiagram`, `IOBR`, `multiMiR`, `enrichR`, `GseaVis`, `Seurat`, `harmony`, `SingleR`,
  `celldex`, `monocle3`, `monocle`, `CellChat`, `AUCell`, `scMetabolism` and `export`.
  The exact list of packages required by each script is given by its `library()` calls — e.g.
  `grep -h "library(" Rscript/*.R | sort -u`.
* Two script-specific requirements (see `assets/README.md`):
  * `R.05.113ML.R` uses the local machine-learning function set `assets/113MLtoolkit/`;
  * `R.14.scRNA.R` needs the GSE225948 raw count matrices placed in
    `14.scRNA/00.rawdata/GSE225948/` (public GEO data).

Run the scripts in order from the repository root:

```bash
Rscript Rscript/R.00.rawdata.R     # download/curate input data
Rscript Rscript/R.01.limma.R       # differential expression
Rscript Rscript/R.02.venn.R        # MLM ∩ DEGs
Rscript Rscript/R.03.enrich.R      # GO / KEGG enrichment
Rscript Rscript/R.04.PPI.R         # PPI network (STRING)
Rscript Rscript/R.05.113ML.R       # machine-learning model selection and validation
Rscript Rscript/R.06.expressionROC.R
Rscript Rscript/R.07.nomogram.R    # nomogram, calibration, DCA, ROC
Rscript Rscript/R.08.location.R    # subcellular localization
Rscript Rscript/R.09.GSEA.R        # single-gene GSEA
Rscript Rscript/R.10.cibersort.R   # immune cell deconvolution
Rscript Rscript/R.11.multiMiR.R    # miRNA / lncRNA / TF networks
Rscript Rscript/R.13.enrichR.R     # drug prediction
Rscript Rscript/R.14.scRNA.R       # single-cell analysis
```

Each script writes its outputs into the correspondingly numbered folder. Some modules query public
web services at run time (STRING, Enrichr, ENCORI, KEGG, Ensembl/BioMart); results may therefore
depend on the version of those resources at the time of the run.

> Note: package versions for the deposited results were not captured at the time of the original
> analysis. The package list above reflects the functions used, not a frozen environment.

---

## 4. Licensing

| Content | License |
|---|---|
| Analysis code (`Rscript/*.R`) | **MIT** — see [`LICENSE`](LICENSE) |
| Derived result tables and figures (files under `01.*`–`14.*` produced by this pipeline) | **CC BY 4.0** — https://creativecommons.org/licenses/by/4.0/ |
| Third-party data under `00.rawdata/` (GEO files) and gene sets compiled from published studies | **Not covered** by the licenses above; subject to the original providers' terms — please cite the sources in Section 2 |

The copyright holder is currently indicated by the repository account name; it will be updated to
the authors'/institution's name when the accompanying manuscript is published.

---

## 5. Citation

Citation details (authors, manuscript title, journal, DOI) will be added here once the accompanying
manuscript is published. In the meantime, please cite this repository by its URL and the GEO
accessions listed in Section 2.

Questions, bug reports and requests for the additional assets mentioned in Section 3 are welcome
via the repository's issue tracker.
