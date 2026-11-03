# 15b_fibroblast_state.R
# 接 15_scRNA_revision.R：SCG2 标记的成纤维细胞亚群（F1）在各组中的比例，及其标志基因
# 输入：data/processed/scRNA_GSE235711_fibroblasts.rds
# 输出：results/tables/rev_sc1_F1_*.csv

suppressPackageStartupMessages({ library(tidyverse); library(Seurat) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
tab <- "results/tables"
fb <- readRDS("data/processed/scRNA_GSE235711_fibroblasts.rds")

# 1. 亚群比例：每个样本 F1 占成纤维细胞的比例，组间精确 Wilcoxon 检验
fr <- fb@meta.data %>% group_by(sample, group) %>% summarise(n = n(), F1 = mean(fsub == "F1"), .groups = "drop")
cmp <- list(c("NP", "Control"), c("NP", "CRSwNP-Eth"), c("NP", "CRSsNP"), c("CRSwNP-Eth", "Control"), c("CRSsNP", "Control"))
resF <- map_dfr(cmp, function(p) {
  a <- fr$F1[fr$group == p[1]]; b <- fr$F1[fr$group == p[2]]
  tibble(comparison = paste(p, collapse = " vs "), median_1 = median(a), median_0 = median(b),
         n1 = length(a), n0 = length(b), P_exact = wilcox.test(a, b, exact = TRUE)$p.value)
})
print(fr %>% arrange(group)); print(resF)
write_csv(fr, file.path(tab, "rev_sc1_F1_fraction_per_sample.csv"))
write_csv(resF, file.path(tab, "rev_sc1_F1_fraction_tests.csv"))

# 2. 去掉单个样本主导的 F0 后，F1 与其他成纤维细胞亚群相比的标志基因
Idents(fb) <- "fsub"
core <- subset(fb, fsub %in% c("F1", "F2", "F3", "F4"))
mk <- FindMarkers(core, ident.1 = "F1", min.pct = 0.2, logfc.threshold = 0.5, only.pos = TRUE, verbose = FALSE) %>%
  rownames_to_column("gene") %>% arrange(desc(avg_log2FC))
print(head(mk, 25))
write_csv(mk, file.path(tab, "rev_sc1_F1_markers.csv"))

# 3. F1 内部 SCG2/COL10A1 阳性率在各组的差异（每样本）
x <- GetAssayData(fb, layer = "data")
within <- fb@meta.data %>% mutate(S = x["SCG2", ] > 0, C = x["COL10A1", ] > 0) %>% filter(fsub == "F1") %>%
  group_by(sample, group) %>% summarise(n = n(), pct_SCG2 = mean(S) * 100, pct_COL10A1 = mean(C) * 100, .groups = "drop")
print(within %>% group_by(group) %>% summarise(across(c(n, pct_SCG2, pct_COL10A1), median)))
write_csv(within, file.path(tab, "rev_sc1_F1_positivity.csv"))
