# CRSwNP multi-cohort transcriptomics and single-cell analysis

Analysis code for the manuscript:

> **Multi-cohort transcriptomic and single-cell analysis identifies polyp-associated expression of *COL10A1* in fibroblasts in chronic rhinosinusitis with nasal polyps**
> Yitao Wu, Xinjiang Medical University

All data are public (GEO). The scripts reproduce every table and figure in the manuscript from raw GEO downloads.

## Data

Datasets were identified by a systematic GEO search (`scripts/16_geo_search.R`; screening decisions in `results/tables/rev_geo_screening.csv`).

| Accession | Type | Use |
|---|---|---|
| GSE136825 | bulk RNA-seq | training; meta-analysis; paired polyp vs own turbinate |
| GSE179265 | bulk RNA-seq | training; meta-analysis |
| GSE36830 | microarray (GPL570) | validation; meta-analysis; CRSsNP / CRSwNP uncinate specificity |
| GSE23552 | microarray (GPL5175) | validation; meta-analysis; non-polyp CRSwNP mucosa and allergic rhinitis specificity |
| GSE194282 | microarray (GPL17692) | validation; meta-analysis |
| GSE338426 | bulk RNA-seq | post hoc validation; eosinophilic vs non-eosinophilic polyps |
| GSE72713 | bulk RNA-seq (FPKM) | post hoc validation (n = 3 per group) |
| GSE235711 | scRNA-seq (10x) | single-cell discovery, all 20 samples |
| GSE276503 | scRNA-seq (DNBelab C4) | single-cell replication; 4 samples with corrupted GEO files excluded; analysed at inferred-donor level |

## Pipeline

Run in order from the repository root.

| Script | Step |
|---|---|
| `00_install_packages.R` | install R/Bioconductor packages |
| `01_download.sh`, `01b_download_scRNA.sh`, `01c_download_scRNA2.sh` | download bulk data, gene sets, KEGG snapshot and both single-cell datasets (resumable) |
| `02_preprocess.R` | sample selection, gene annotation |
| `03_deg_meta.R` | per-cohort limma / limma-voom, Hedges' g, DerSimonian–Laird meta-analysis |
| `04_enrichment.R` | GO/KEGG over-representation, Hallmark GSEA |
| `05_ml_hub_genes.R` | LASSO stability selection, random forest, SVM-RFE (training cohorts only); validation AUCs |
| `06_scRNA.R` | first single-cell analysis (superseded by 15 for the manuscript) |
| `07_fetch_references.R` | bibliography from verified PubMed IDs |
| `08_robustness.R` | duplicate samples, threshold grid, fibroblast-score adjustment |
| `14_specificity_baselines.R` | REML + Hartung–Knapp meta-analysis; random-gene and classic-marker AUC baselines; CRSsNP / non-polyp CRSwNP / allergic rhinitis specificity; paired polyp vs own turbinate |
| `15_scRNA_revision.R`, `15b_fibroblast_state.R`, `15c_sensitivity_eth3.R` | GSE235711, all 20 samples: stromal re-annotation, counts-sum pseudobulk (edgeR QL, exact Wilcoxon), co-expression, exploratory fibroblast subclusters, sensitivity analysis |
| `16_geo_search.R` | systematic GEO search |
| `17_new_bulk_cohorts.R` | GSE338426 and GSE72713 |
| `18_scRNA_replication_GSE276503.R`, `18b_replicate_check.R`, `18c_donor_level.R` | GSE276503 replication, technical-replicate detection, donor-level pseudobulk |
| `12_workflow_figure.R`, `19_revision_figures.R`, `09_plos_figures.R` | figures |
| `20_supporting_files.R` | Supporting information files |

`results/tables/` contains all result tables. `sessionInfo.txt` records package versions (R 4.6.1).

## Notes on the analysis

- Gene selection uses only the two RNA-seq cohorts; the three microarray cohorts are used once for validation.
- Separation of polyps from normal mucosa is compared against random gene sets and classic markers; the consensus genes are not distinctive in this respect.
- The DerSimonian–Laird meta-analysis is descriptive; REML with Hartung–Knapp adjustment is reported alongside it.
- Single-cell statistics use pseudobulk counts per sample, or per inferred donor in GSE276503.
- Post hoc changes and analyses (meta-analysis inclusion rule, systematic search, baselines, specificity, second single-cell dataset) are identified as such in the manuscript.
- An AI assistant (Claude, Anthropic) wrote and ran the code under the author's direction; see the manuscript Methods.

## License

MIT (see `LICENSE`).
