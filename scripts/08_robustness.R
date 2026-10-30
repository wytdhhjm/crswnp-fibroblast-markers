# 08_robustness.R
# 自查与稳健性分析（回应审稿人最可能的质疑）
#   A. 同一患者重复取样：GSE136825 中 PY/TT 编号相同的息肉是否来自同一患者？
#      -> 敏感性分析：每个编号只保留一个息肉样本后重做训练集差异分析
#   B. 机器学习选择对阈值的依赖：RF/SVM-RFE 前 20/30/50、LASSO 频率 0.4/0.5/0.6 下核心基因是否稳定
#   C. 外部验证 AUC 的置信区间改用 bootstrap（AUC=1 时 DeLong 区间退化为 1-1，不可信）
#   D. 组织成分混杂：bulk 中 SCG2/COL10A1 升高是否只是成纤维细胞比例升高？
#      -> 在 bulk 中用成纤维细胞标志基因评分做协变量，检验核心基因是否仍显著
#   E. 单细胞：成纤维细胞特异性用 AUC（标志基因检验）量化；核心基因 vs 成纤维细胞比例
#
# 输出：results/tables/robust_*.csv

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(pROC)
  library(Seurat)
})

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
proc <- "data/processed"
tab  <- "results/tables"
set.seed(2026)
hub <- read_csv(file.path(tab, "ML_feature_selection.csv"), show_col_types = FALSE) %>% filter(hub)
hub_genes <- hub$gene
logcpm <- function(cnt) { y <- normLibSizes(DGEList(cnt)); cpm(y, log = TRUE, prior.count = 1) }

# ---------- A. 同一编号重复样本 ----------
d1 <- readRDS(file.path(proc, "GSE136825.rds"))
ph <- d1$pheno %>% mutate(prefix = sub(" .*", "", title), num = sub(".* ", "", title))
dup_ids <- ph %>% filter(group == "NP") %>% count(num) %>% filter(n > 1) %>% pull(num)
cat("A. GSE136825 中编号重复的息肉样本：", paste0(dup_ids, collapse = ", "), "\n")
keep <- !(ph$group == "NP" & ph$prefix == "TT" & ph$num %in% dup_ids)   # 重复编号只保留 PY
cnt <- d1$counts[, keep]; grp <- droplevels(ph$group[keep])
y <- DGEList(cnt, group = grp); y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
y <- normLibSizes(y); des <- model.matrix(~ grp)
tt <- topTable(eBayes(lmFit(voom(y, des), des)), coef = 2, number = Inf, sort.by = "none")
resA <- tibble(gene = hub_genes, logFC = tt[hub_genes, "logFC"], FDR = tt[hub_genes, "adj.P.Val"])
cat("   去重后 NP =", sum(grp == "NP"), "，对照 =", sum(grp == "Control"), "\n"); print(resA)
write_csv(resA %>% mutate(n_np = sum(grp == "NP"), n_ctrl = sum(grp == "Control")),
          file.path(tab, "robust_A_dedup_GSE136825.csv"))

# ---------- B. 阈值敏感性 ----------
ml <- read_csv(file.path(tab, "ML_feature_selection.csv"), show_col_types = FALSE)
grid <- expand_grid(lasso = c(0.4, 0.5, 0.6), top = c(20, 30, 50))
resB <- grid %>% mutate(hub = map2_chr(lasso, top, function(l, k)
  paste(ml$gene[ml$lasso_freq >= l & ml$rf_rank <= k & ml$svm_rank <= k], collapse = ", ")))
cat("\nB. 不同阈值下三算法交集：\n"); print(resB, n = 20)
write_csv(resB, file.path(tab, "robust_B_threshold_grid.csv"))

# ---------- C. bootstrap AUC 置信区间 ----------
zrow <- function(m) t(scale(t(m)))
val <- c("GSE36830", "GSE23552", "GSE194282")
score <- function(z) colMeans(z[intersect(hub_genes, rownames(z)), , drop = FALSE])
resC <- map_dfr(val, function(i) {
  d <- readRDS(file.path(proc, paste0(i, ".rds"))); z <- zrow(d$expr)
  r <- roc(d$pheno$group, score(z), levels = c("Control", "NP"), direction = "<", quiet = TRUE)
  ci <- ci.auc(r, method = "bootstrap", boot.n = 2000, boot.stratified = TRUE)
  tibble(dataset = i, AUC = ci[2], lo = ci[1], hi = ci[3])
})
cat("\nC. 四基因评分 AUC（分层 bootstrap 2000 次 95% CI）：\n"); print(resC)
write_csv(resC, file.path(tab, "robust_C_bootstrap_AUC.csv"))

# ---------- D. 校正成纤维细胞含量 ----------
fib_markers <- c("COL1A1", "COL1A2", "COL3A1", "DCN", "LUM", "PDGFRA")
resD <- map_dfr(c("GSE136825", "GSE179265", "GSE36830", "GSE23552", "GSE194282"), function(i) {
  d <- readRDS(file.path(proc, paste0(i, ".rds")))
  e <- if (!is.null(d$counts)) logcpm(d$counts) else d$expr
  fm <- intersect(fib_markers, rownames(e))
  fib <- colMeans(zrow(e[fm, , drop = FALSE]))
  grp <- d$pheno$group
  fib_diff <- t.test(fib[grp == "NP"], fib[grp == "Control"])$p.value
  map_dfr(intersect(hub_genes, rownames(e)), function(g) {
    m0 <- lm(e[g, ] ~ grp); m1 <- lm(e[g, ] ~ grp + fib)
    tibble(dataset = i, gene = g,
           beta_unadj = coef(m0)[2], p_unadj = summary(m0)$coef[2, 4],
           beta_adj_fib = coef(m1)[2], p_adj_fib = summary(m1)$coef[2, 4],
           fib_score_NP_vs_ctrl_p = fib_diff)
  })
})
cat("\nD. 校正成纤维细胞标志基因评分后（各数据集）：\n")
print(resD %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 30)
write_csv(resD, file.path(tab, "robust_D_fibroblast_adjusted.csv"))

# ---------- E. 单细胞：特异性与细胞比例 ----------
sc <- readRDS(file.path(proc, "scRNA_GSE235711.rds"))
Idents(sc) <- "cell_type"
resE1 <- map_dfr(hub_genes, function(g) {
  x <- GetAssayData(sc, layer = "data")[g, ]
  map_dfr(levels(sc$cell_type), function(ct) {
    inside <- sc$cell_type == ct
    if (sum(inside) < 50) return(NULL)
    tibble(gene = g, cell_type = ct, pct_in = mean(x[inside] > 0) * 100, pct_out = mean(x[!inside] > 0) * 100,
           auc = as.numeric(auc(roc(as.integer(inside), x, quiet = TRUE, direction = "<"))))
  }) %>% slice_max(auc, n = 1)
})
cat("\nE1. 每个核心基因区分度最高的细胞类型（AUC：该细胞类型 vs 其他所有细胞）：\n"); print(resE1)
write_csv(resE1, file.path(tab, "robust_E1_cell_specificity.csv"))

prop <- sc@meta.data %>% count(sample, group, cell_type) %>% group_by(sample) %>%
  mutate(frac = n / sum(n)) %>% ungroup() %>% filter(cell_type == "Fibroblast")
resE2 <- prop %>% group_by(group) %>% summarise(median_frac = median(frac), n = n(), .groups = "drop")
pE2 <- wilcox.test(frac ~ droplevels(group), data = prop, exact = FALSE)$p.value
cat("\nE2. 成纤维细胞占比（每个样本）：\n"); print(resE2); cat("   Wilcoxon P =", signif(pE2, 3), "\n")
write_csv(prop, file.path(tab, "robust_E2_fibroblast_fraction.csv"))

cat("\n完成\n")
