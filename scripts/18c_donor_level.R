# 18c_donor_level.R
# GSE276503 样本合并到"推定供者"层面后重做成纤维细胞比较（18b 发现大量技术重复）
#   推定供者：成纤维细胞伪 bulk（2,000 高变基因）Spearman 相关 ≥ 阈值的样本合并（平均连锁聚类）
#   主分析阈值 0.94（高于全部样本对相关性的第 95 百分位 0.883 很多），敏感性分析 0.97、0.88
#   同一推定供者的 counts 求和后，再按供者做 edgeR QL 与精确 Wilcoxon；健康对照合并后可能只剩 1 个供者，不做统计
# 输出：results/tables/rev_sc2_donor_level_*.csv
suppressPackageStartupMessages({ library(tidyverse); library(edgeR) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
tab <- "results/tables"
st <- readRDS("data/processed/scRNA_GSE276503_stromal.rds")
fibcells <- st$sub_type == "Fibroblast"
cnt <- SeuratObject::GetAssayData(st, layer = "counts")[, fibcells]
md <- st@meta.data[fibcells, c("sample", "group")]
ln <- SeuratObject::GetAssayData(st, layer = "data")[c("SCG2", "COL10A1"), fibcells]
r <- read_csv(file.path(tab, "rev_sc2_sample_correlation.csv"), show_col_types = FALSE) %>% column_to_rownames("sample") %>% as.matrix()
hc <- hclust(as.dist(1 - r), "average")
genes <- c("SCG2", "COL10A1", "PTHLH", "PDE7A")

run <- function(th) {
  donor <- cutree(hc, h = 1 - th)
  dmap <- tibble(sample = names(donor), donor = paste0("D", donor)) %>%
    mutate(group = as.character(md$group[match(sample, md$sample)]))
  # 一个推定供者内若混有不同组（不同部位）则不合并：按 供者×组 定义单位
  dmap <- dmap %>% mutate(unit = paste(donor, group, sep = "|"))
  unit_of_cell <- dmap$unit[match(md$sample, dmap$sample)]
  units <- dmap %>% distinct(unit, group) %>% arrange(group)
  pb <- sapply(units$unit, function(u) Matrix::rowSums(cnt[, unit_of_cell == u, drop = FALSE]))
  grp <- factor(make.names(units$group))
  y <- DGEList(pb, group = grp); y <- normLibSizes(y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE])
  des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
  fit <- glmQLFit(estimateDisp(y, des), des, robust = TRUE)
  lc <- cpm(y, log = TRUE, prior.count = 1)
  cons <- c(`NP vs CRSwNP IT` = "NP - CRSwNP.IT", `NP vs CRSwNP MT` = "NP - CRSwNP.MT")
  res <- map_dfr(names(cons), function(nm) {
    p <- strsplit(cons[[nm]], " - ")[[1]]
    tt <- topTags(glmQLFTest(fit, contrast = makeContrasts(contrasts = cons[[nm]], levels = des)), n = Inf)$table
    a1 <- grp == p[1]; a0 <- grp == p[2]
    tibble(threshold = th, comparison = nm, gene = genes, n_donor_1 = sum(a1), n_donor_0 = sum(a0),
           logFC = tt[genes, "logFC"], P_QL = tt[genes, "PValue"], FDR_genomewide = tt[genes, "FDR"],
           P_exact_wilcoxon = map_dbl(genes, ~ wilcox.test(lc[.x, a1], lc[.x, a0], exact = TRUE)$p.value))
  })
  pos <- tibble(unit = unit_of_cell, S = ln["SCG2", ] > 0, C = ln["COL10A1", ] > 0) %>% group_by(unit) %>%
    summarise(n_cells = n(), pct_SCG2 = mean(S) * 100, pct_COL10A1 = mean(C) * 100, .groups = "drop") %>%
    left_join(units, by = "unit") %>% mutate(threshold = th)
  list(res = res, pos = pos, map = dmap %>% mutate(threshold = th))
}
out <- map(c(0.94, 0.97, 0.88), run)
res <- map_dfr(out, "res"); pos <- map_dfr(out, "pos"); dm <- map_dfr(out, "map")
write_csv(res, file.path(tab, "rev_sc2_donor_level_pseudobulk.csv"))
write_csv(pos, file.path(tab, "rev_sc2_donor_level_positivity.csv"))
write_csv(dm, file.path(tab, "rev_sc2_donor_map.csv"))
print(dm %>% filter(threshold == 0.94) %>% group_by(unit) %>% summarise(samples = paste(sample, collapse = ",")), n = 40)
print(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 40)
print(pos %>% filter(threshold == 0.94) %>% arrange(group) %>% mutate(across(where(is.numeric), ~ round(.x, 1))), n = 40)
