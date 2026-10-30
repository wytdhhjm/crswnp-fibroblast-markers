# CRSwNP multi-cohort transcriptomics and single-cell mapping

Analysis code for the manuscript:

> **Multi-cohort transcriptomic meta-analysis and single-cell mapping identify fibroblast-restricted *SCG2* and *COL10A1* as robust markers of chronic rhinosinusitis with nasal polyps**
> Yitao Wu, Xinjiang Medical University

All data are public (GEO). The scripts reproduce every table and figure in the manuscript from raw GEO downloads.

## Data

| Accession | Type | Use |
|---|---|---|
| GSE136825 | bulk RNA-seq | training / meta-analysis |
| GSE179265 | bulk RNA-seq | training / meta-analysis |
| GSE36830 | microarray (GPL570) | external validation / meta-analysis |
| GSE23552 | microarray (GPL5175) | external validation / meta-analysis |
| GSE194282 | microarray (GPL17692) | external validation / meta-analysis |
| GSE235711 | scRNA-seq (10x) | cell-type localisation |

## Pipeline

Run in order from the repository root:

| Script | Step |
|---|---|
| `scripts/00_install_packages.R` | install R/Bioconductor packages |
| `scripts/01_download.sh` | download bulk data, Hallmark gene sets and KEGG snapshot (chunked HTTPS, resumable) |
| `scripts/01b_download_scRNA.sh` | download GSE235711 filtered 10x matrices |
| `scripts/02_preprocess.R` | sample selection, gene annotation, per-dataset objects |
| `scripts/03_deg_meta.R` | per-cohort limma / limma-voom, Hedges' g, DerSimonian–Laird random-effects meta-analysis, leave-one-out and sensitivity analyses |
| `scripts/04_enrichment.R` | GO/KEGG over-representation and Hallmark GSEA |
| `scripts/05_ml_hub_genes.R` | LASSO stability selection, random forest, SVM-RFE on training cohorts only; external validation (AUC, DeLong CI) |
| `scripts/06_scRNA.R` | Seurat + Harmony integration, cell-type annotation, sample-level pseudobulk tests |
| `scripts/07_fetch_references.R` | builds the bibliography from PubMed IDs |
| `scripts/08_robustness.R` | robustness checks: duplicate-sample exclusion, threshold grid, bootstrap AUC CIs, fibroblast-content adjustment, single-cell specificity |

`results/tables/` contains all result tables produced by the scripts.
`sessionInfo.txt` records the R and package versions used (R 4.6.1).

## Design notes

- **No information leakage:** machine-learning feature selection uses only the two RNA-seq cohorts; the three microarray cohorts are used once, for final validation.
- **Cross-platform meta-analysis** pools standardised effect sizes rather than merged expression values.
- **Single-cell statistics** are computed at the sample level (pseudobulk) to avoid pseudoreplication.
- The meta-analysis inclusion rule was changed after inspection from "measured in all 5 cohorts" to "at least 4"; results under the original rule are written to `meta_sensitivity_k5.csv`.

## License

MIT (see `LICENSE`).
