# 20_supporting_files.R
# 按修订稿（2026-11-02）的 Supporting information 编号生成补充文件
# 多张表的补充文件用 Excel（每个来源表一个工作表），单一表用 CSV
# 输出：manuscript_plos/supporting/S1_Fig.tif、S1_Table.csv ... S11_Table.xlsx

suppressPackageStartupMessages({ library(tidyverse); library(writexl) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
tab <- "results/tables"; out <- "manuscript_plos/supporting"
unlink(list.files(out, full.names = TRUE)); dir.create(out, showWarnings = FALSE)
rd <- function(f) read_csv(file.path(tab, f), show_col_types = FALSE)
file.copy("manuscript_plos/figures/S1_Fig.tif", file.path(out, "S1_Fig.tif"))

write_csv(rd("DEG_per_dataset_summary.csv"), file.path(out, "S1_Table.csv"))
write_csv(rd("rev_geo_screening.csv"), file.path(out, "S2_Table.csv"))
write_csv(rd("robust_A_dedup_GSE136825.csv"), file.path(out, "S3_Table.csv"))
write_csv(rd("robust_B_threshold_grid.csv"), file.path(out, "S4_Table.csv"))
write_csv(rd("robust_D_fibroblast_adjusted.csv"), file.path(out, "S5_Table.csv"))
write_csv(rd("rev_A_meta_REML_HK.csv") %>%
            left_join(rd("meta_all_genes.csv") %>% select(gene, tau2, Q, Q_P, I2), by = "gene") %>%
            left_join(rd("meta_robust_DEGs.csv") %>% select(gene, loo_stable = still_sig_all), by = "gene") %>%
            rename(g_DL = g_pooled, ci_low_DL = ci_low, ci_high_DL = ci_high, FDR_DL = FDR, direction_DL = direction),
          file.path(out, "S6_Table.csv"))
write_xlsx(list(reference_sets = rd("rev_B_reference_AUC.csv"), random4_summary = rd("rev_B_random4_summary.csv"),
                per_gene_AUC = rd("ML_hub_gene_AUC.csv") %>% filter(set == "External validation"),
                GSE338426_AUC = rd("rev_E_new_cohorts_AUC.csv")), file.path(out, "S7_Table.xlsx"))
write_xlsx(list(GSE36830_GSE23552 = rd("rev_C_specificity.csv"), GSE136825_paired_site = rd("rev_D_site_paired_GSE136825.csv"),
                GSE338426_GSE72713 = rd("rev_E_new_cohorts_DE.csv")), file.path(out, "S8_Table.xlsx"))
write_xlsx(list(pseudobulk_edgeR = rd("rev_sc1_pseudobulk_edgeR.csv"), without_CRSwNP_Eth_3 = rd("rev_sc1_pseudobulk_without_Eth3.csv"),
                cells_per_sample_type = rd("rev_sc1_cells_per_sample_type.csv"), umi_share = rd("rev_sc1_umi_share_by_type.csv"),
                expression_by_type = rd("rev_sc1_expression_by_type.csv"), coexpression = rd("rev_sc1_coexpression.csv")),
           file.path(out, "S9_Table.xlsx"))
write_xlsx(list(donor_map = rd("rev_sc2_donor_map.csv"), donor_pseudobulk = rd("rev_sc2_donor_level_pseudobulk.csv"),
                donor_positivity = rd("rev_sc2_donor_level_positivity.csv"), sample_QC = rd("rev_sc2_QC.csv"),
                excluded_corrupt = tibble(sample = readLines(file.path(tab, "rev_sc2_excluded_corrupt_samples.txt"))),
                stromal_expression = rd("rev_sc2_expression_stromal_types.csv"), coexpression = rd("rev_sc2_coexpression.csv")),
           file.path(out, "S10_Table.xlsx"))
write_xlsx(list(GSE235711_subclusters = rd("rev_sc1_fibroblast_subclusters.csv"), GSE235711_F1_markers = rd("rev_sc1_F1_markers.csv"),
                GSE235711_F1_fraction = rd("rev_sc1_F1_fraction_per_sample.csv"),
                GSE276503_per_sample = rd("rev_sc2_fibroblast_per_sample.csv"), GSE276503_tests = rd("rev_sc2_fibroblast_state_tests.csv")),
           file.path(out, "S11_Table.xlsx"))
print(list.files(out))
