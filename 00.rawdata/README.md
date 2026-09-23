# 00.rawdata — input data and provenance

This folder holds the input data used by the analysis pipeline.

| Item | Source | Notes |
|---|---|---|
| `GSE16561/` | NCBI GEO [GSE16561](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16561) (GPL6883) | human whole-blood expression, discovery cohort (63 samples: 24 control / 39 stroke); original study PMID [20837969](https://pubmed.ncbi.nlm.nih.gov/20837969/) |
| `GSE58294/` | NCBI GEO [GSE58294](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58294) (GPL570) | human expression, validation cohort (92 samples: 23 control / 69 stroke); original study PMID [25036109](https://pubmed.ncbi.nlm.nih.gov/25036109/) |
| `01.GSE225948_group.csv` | NCBI GEO [GSE225948](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE225948) | sample grouping for the mouse single-cell module (experimental stroke); original study PMID [38177281](https://pubmed.ncbi.nlm.nih.gov/38177281/) |
| `线粒体和脂质代谢相关基因.xlsx` | gene sets compiled from published studies | two sheets: lipid-metabolism-related genes (PMID [39827204](https://pubmed.ncbi.nlm.nih.gov/39827204/)) and mitochondrial-related genes (PMID [36915111](https://pubmed.ncbi.nlm.nih.gov/36915111/)); see the root `README.md`, Section 2 |
| `GPL/` | created at run time by `Rscript/R.00.rawdata.R` | GEO platform annotation cache; not stored in this repository |

## License notice

The GEO-derived files listed above are **not** covered by the MIT license of the analysis code, nor
by the CC BY 4.0 license of the derived result tables and figures. They are redistributed here only
to make the analysis reproducible and remain subject to the terms of NCBI GEO and to the copyright
of the original data submitters. Please cite the original studies when re-using them.
