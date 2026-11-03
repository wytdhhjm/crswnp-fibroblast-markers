# 15c_sensitivity_eth3.R
# 敏感性分析：GSE235711 中 CRSwNP_Eth_3 贡献了 1,917 个成纤维细胞（92% 属于以线粒体假基因/组蛋白基因为标志的
# F0 亚群，提示低质量细胞），去掉该样本后重做成纤维细胞伪 bulk（counts 求和 + edgeR QL）
# 输出：results/tables/rev_sc1_pseudobulk_without_Eth3.csv
suppressPackageStartupMessages({ library(tidyverse); library(Seurat); library(edgeR) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
sc <- readRDS("data/processed/scRNA_GSE235711_all20.rds")
fib <- sc$cell_type == "Fibroblast" & sc$sample != "CRSwNP_Eth_3"
cn <- GetAssayData(sc, layer = "counts")[, fib]; md <- sc@meta.data[fib, ]
ok <- names(which(table(md$sample) >= 20))
pb <- sapply(ok, function(s) Matrix::rowSums(cn[, md$sample == s, drop = FALSE]))
grp <- factor(make.names(as.character(md$group[match(ok, md$sample)])))
y <- DGEList(pb, group = grp); y <- normLibSizes(y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE])
des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
fit <- glmQLFit(estimateDisp(y, des), des, robust = TRUE); lc <- cpm(y, log = TRUE, prior.count = 1)
cons <- c(`NP vs Control` = "NP - Control", `NP vs CRSwNP-Eth` = "NP - CRSwNP.Eth", `CRSwNP-Eth vs Control` = "CRSwNP.Eth - Control")
res <- map_dfr(names(cons), function(nm) {
  p <- strsplit(cons[[nm]], " - ")[[1]]
  tt <- topTags(glmQLFTest(fit, contrast = makeContrasts(contrasts = cons[[nm]], levels = des)), n = Inf)$table
  g <- c("SCG2", "COL10A1"); a1 <- grp == p[1]; a0 <- grp == p[2]
  tibble(comparison = nm, gene = g, n1 = sum(a1), n0 = sum(a0), logFC = tt[g, "logFC"], P_QL = tt[g, "PValue"],
         FDR_genomewide = tt[g, "FDR"], P_exact = map_dbl(g, ~ wilcox.test(lc[.x, a1], lc[.x, a0], exact = TRUE)$p.value))
})
print(res %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))
write_csv(res, "results/tables/rev_sc1_pseudobulk_without_Eth3.csv")
