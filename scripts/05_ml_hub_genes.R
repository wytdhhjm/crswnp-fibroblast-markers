# 05_ml_hub_genes.R
# 第 2 步（下）：机器学习筛选核心基因 + 外部独立验证
#
# 避免数据泄漏（审稿重点）：
#   发现阶段只用 2 个 RNA-seq 训练集（GSE136825、GSE179265）；
#   3 个芯片数据集（GSE36830、GSE23552、GSE194282）不参与任何筛选，只在最后做一次验证。
#   （03 的 5 队列荟萃分析是描述性的全景结果，不用于本脚本的候选基因筛选）
#
# 流程：
#   1) 候选基因：两个训练集中均 FDR < 0.05、|log2FC| > 0.58（1.5 倍）、方向一致
#   2) 合并训练集：每个数据集内部对每个基因做 z-score 标准化
#      （不用 ComBat：ComBat 需要输入分组信息，会把分组信息带进后续交叉验证）
#   3) 三种算法，均基于重复抽样以保证稳定：
#      - LASSO 稳定性选择：200 次分层 bootstrap，统计每个基因被选中的频率
#      - 随机森林：20 个随机种子，置换重要性（MeanDecreaseAccuracy）取平均
#      - SVM-RFE（线性核）：每轮剔除权重最小的 10% 基因，得到排序
#   4) 核心基因 = LASSO 频率 >= 0.5 ∩ RF 前 30 ∩ SVM-RFE 前 30
#   5) 外部验证：每个基因在 3 个芯片数据集中的 AUC（DeLong 95% CI）；
#      再用"核心基因评分"（数据集内 z-score 后，上调基因均值 − 下调基因均值）
#      评估整体诊断效能——不需要跨平台迁移模型系数
#
# 输出：results/tables/ML_*.csv、results/figures/Fig3_*

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(glmnet)
  library(randomForest)
  library(e1071)
  library(pROC)
  library(patchwork)
})

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
proc <- "data/processed"
tab  <- "results/tables"
fig  <- "results/figures"
set.seed(2026)

col_group <- c(Control = "#2a78d6", NP = "#eb6834")
col_alg   <- c(LASSO = "#2a78d6", RF = "#eb6834", `SVM-RFE` = "#1baf7a")
ink2 <- "#52514e"
theme_pub <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.ticks = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.text = element_text(colour = ink2),
          plot.title = element_text(face = "bold", size = base + 1),
          legend.key.size = unit(3, "mm"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}

train_ids <- c("GSE136825", "GSE179265")
val_ids   <- c("GSE36830", "GSE23552", "GSE194282")

# ---------- 1. 候选基因（仅训练集） ----------
deg <- map(setNames(train_ids, train_ids),
           ~ read_csv(file.path(tab, paste0("DEG_", .x, ".csv")), show_col_types = FALSE))
cand <- inner_join(deg[[1]], deg[[2]], by = "gene", suffix = c(".a", ".b")) %>%
  filter(FDR.a < 0.05, FDR.b < 0.05, sign(logFC.a) == sign(logFC.b),
         abs(logFC.a) > 0.58, abs(logFC.b) > 0.58) %>%
  mutate(direction = ifelse(logFC.a > 0, "Up", "Down"))
cat("训练集候选基因：", nrow(cand), "（上调", sum(cand$direction == "Up"),
    "，下调", sum(cand$direction == "Down"), "）\n")

# ---------- 2. 合并训练集 ----------
zrow <- function(m) t(scale(t(m)))
logcpm <- function(d) {
  y <- DGEList(d$counts); y <- normLibSizes(y)
  cpm(y, log = TRUE, prior.count = 1)
}
tr <- map(train_ids, function(i) {
  d <- readRDS(file.path(proc, paste0(i, ".rds")))
  list(z = zrow(logcpm(d)[cand$gene, ]), group = d$pheno$group)
})
X <- t(do.call(cbind, map(tr, "z")))
y <- factor(unlist(map(tr, ~ as.character(.x$group))), c("Control", "NP"))
batch <- rep(train_ids, map_int(tr, ~ length(.x$group)))
stopifnot(!anyNA(X))
cat("合并训练集：", nrow(X), "个样本（NP", sum(y == "NP"), "，对照", sum(y == "Control"), "）\n")

# ---------- 3a. LASSO 稳定性选择 ----------
B <- 200
strat_boot <- function() unlist(lapply(split(seq_along(y), interaction(y, batch)),
                                       function(ix) sample(ix, length(ix), replace = TRUE)))
lasso_sel <- map(seq_len(B), function(b) {
  ix <- strat_boot()
  cv <- cv.glmnet(X[ix, ], y[ix], family = "binomial", alpha = 1, nfolds = 5)
  cf <- coef(cv, s = "lambda.1se")[-1, 1]
  names(cf)[cf != 0]
})
lasso_freq <- tibble(gene = colnames(X),
                     lasso_freq = map_dbl(colnames(X), ~ mean(map_lgl(lasso_sel, function(s) .x %in% s))))

# LASSO 全样本拟合（用于画系数路径和交叉验证曲线）
cv_full <- cv.glmnet(X, y, family = "binomial", alpha = 1, nfolds = 10)

# ---------- 3b. 随机森林 ----------
rf_imp <- map_dfc(1:20, function(s) {
  set.seed(s)
  rf <- randomForest(X, y, ntree = 1000, importance = TRUE)
  tibble(!!paste0("seed", s) := importance(rf, type = 1, scale = FALSE)[, 1])
})
rf_tab <- tibble(gene = colnames(X), rf_importance = rowMeans(rf_imp)) %>%
  mutate(rf_rank = rank(-rf_importance, ties.method = "first"))

# ---------- 3c. SVM-RFE（线性核） ----------
svm_rfe_rank <- function(X, y, drop_frac = 0.1) {
  remaining <- colnames(X); eliminated <- character(0)
  while (length(remaining) > 1) {
    m <- svm(X[, remaining, drop = FALSE], y, kernel = "linear", cost = 1, scale = FALSE)
    w <- as.vector(t(m$coefs) %*% m$SV)
    k <- max(1, floor(length(remaining) * drop_frac))
    if (length(remaining) <= 50) k <- 1               # 最后 50 个基因逐个剔除，排序更精细
    out <- remaining[order(w^2)][1:k]
    eliminated <- c(out, eliminated)                  # 先被剔除的排在后面
    remaining <- setdiff(remaining, out)
  }
  c(remaining, eliminated)
}
svm_order <- svm_rfe_rank(X, y)
svm_tab <- tibble(gene = svm_order, svm_rank = seq_along(svm_order))

# ---------- 4. 取交集 ----------
ml <- cand %>% select(gene, direction, logFC.a, logFC.b) %>%
  left_join(lasso_freq, by = "gene") %>% left_join(rf_tab, by = "gene") %>% left_join(svm_tab, by = "gene") %>%
  mutate(sel_lasso = lasso_freq >= 0.5, sel_rf = rf_rank <= 30, sel_svm = svm_rank <= 30,
         hub = sel_lasso & sel_rf & sel_svm) %>%
  arrange(desc(hub), desc(lasso_freq))
write_csv(ml, file.path(tab, "ML_feature_selection.csv"))
hub <- ml %>% filter(hub)
cat("\nLASSO 频率 >= 0.5：", sum(ml$sel_lasso), "个；RF 前 30；SVM-RFE 前 30\n")
cat("三种算法交集（核心基因）：", nrow(hub), "个\n")
print(hub %>% select(gene, direction, lasso_freq, rf_rank, svm_rank))

# ---------- 5. 外部验证 ----------
read_val <- function(i) {
  d <- readRDS(file.path(proc, paste0(i, ".rds")))
  list(z = zrow(d$expr), group = d$pheno$group)
}
sets <- c(setNames(tr, train_ids), map(setNames(val_ids, val_ids), read_val))

auc_ci <- function(score, grp, dir) {
  r <- roc(grp, score, levels = c("Control", "NP"), direction = dir, quiet = TRUE)
  ci <- ci.auc(r, method = "delong")
  tibble(AUC = as.numeric(auc(r)), lo = ci[1], hi = ci[3])
}

gene_auc <- map_dfr(names(sets), function(i) {
  s <- sets[[i]]
  map_dfr(hub$gene, function(g) {
    if (!g %in% rownames(s$z)) return(tibble(gene = g, dataset = i, AUC = NA, lo = NA, hi = NA))
    dir <- ifelse(hub$direction[hub$gene == g] == "Up", "<", ">")
    auc_ci(s$z[g, ], s$group, dir) %>% mutate(gene = g, dataset = i, .before = 1)
  })
}) %>% mutate(set = ifelse(dataset %in% train_ids, "Training", "External validation"))
write_csv(gene_auc, file.path(tab, "ML_hub_gene_AUC.csv"))

hub_score <- function(z) {
  upg <- intersect(hub$gene[hub$direction == "Up"], rownames(z))
  dng <- intersect(hub$gene[hub$direction == "Down"], rownames(z))
  (if (length(upg)) colMeans(z[upg, , drop = FALSE]) else 0) -
    (if (length(dng)) colMeans(z[dng, , drop = FALSE]) else 0)
}
score_auc <- map_dfr(names(sets), function(i) {
  s <- sets[[i]]
  auc_ci(hub_score(s$z), s$group, "<") %>%
    mutate(dataset = i, n_genes = sum(hub$gene %in% rownames(s$z)), .before = 1)
})
write_csv(score_auc, file.path(tab, "ML_hub_score_AUC.csv"))
cat("\n核心基因评分 AUC：\n"); print(score_auc)

# ---------- 6. 作图 ----------
if (nrow(hub) > 0) {
  # 6a LASSO 交叉验证曲线
  cvd <- tibble(loglam = log(cv_full$lambda), dev = cv_full$cvm, lo = cv_full$cvlo, hi = cv_full$cvup,
                nz = cv_full$nzero)
  p_cv <- ggplot(cvd, aes(loglam, dev)) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0, colour = "#c3c2b7", linewidth = 0.3) +
    geom_point(colour = "#e34948", size = 0.9) +
    geom_vline(xintercept = log(c(cv_full$lambda.min, cv_full$lambda.1se)), linetype = 2,
               linewidth = 0.3, colour = "#898781") +
    labs(x = expression(log(lambda)), y = "Binomial deviance", title = "LASSO (10-fold CV)") +
    theme_pub(8)

  # 6b LASSO 稳定性频率（前 20）
  p_lasso <- ml %>% slice_max(lasso_freq, n = 20, with_ties = FALSE) %>%
    ggplot(aes(lasso_freq, fct_reorder(gene, lasso_freq))) +
    geom_col(aes(fill = hub), width = 0.7) +
    geom_vline(xintercept = 0.5, linetype = 2, linewidth = 0.3, colour = "#898781") +
    scale_fill_manual(values = c(`TRUE` = "#2a78d6", `FALSE` = "#c3c2b7"), guide = "none") +
    labs(x = "Selection frequency (200 bootstraps)", y = NULL, title = "LASSO stability") +
    theme_pub(8)

  # 6c RF 重要性（前 20）
  p_rf <- ml %>% slice_min(rf_rank, n = 20, with_ties = FALSE) %>%
    ggplot(aes(rf_importance, fct_reorder(gene, rf_importance))) +
    geom_col(aes(fill = hub), width = 0.7) +
    scale_fill_manual(values = c(`TRUE` = "#eb6834", `FALSE` = "#c3c2b7"), guide = "none") +
    labs(x = "Permutation importance", y = NULL, title = "Random forest") +
    theme_pub(8)

  # 6d 三种算法韦恩图（用集合计数条形表示，避免手绘重叠圆）
  sets_sel <- list(LASSO = ml$gene[ml$sel_lasso], RF = ml$gene[ml$sel_rf], `SVM-RFE` = ml$gene[ml$sel_svm])
  combos <- tibble(
    combo = c("LASSO only", "RF only", "SVM-RFE only", "LASSO ∩ RF", "LASSO ∩ SVM-RFE", "RF ∩ SVM-RFE", "All three"),
    n = c(length(setdiff(sets_sel$LASSO, union(sets_sel$RF, sets_sel$`SVM-RFE`))),
          length(setdiff(sets_sel$RF, union(sets_sel$LASSO, sets_sel$`SVM-RFE`))),
          length(setdiff(sets_sel$`SVM-RFE`, union(sets_sel$LASSO, sets_sel$RF))),
          length(setdiff(intersect(sets_sel$LASSO, sets_sel$RF), sets_sel$`SVM-RFE`)),
          length(setdiff(intersect(sets_sel$LASSO, sets_sel$`SVM-RFE`), sets_sel$RF)),
          length(setdiff(intersect(sets_sel$RF, sets_sel$`SVM-RFE`), sets_sel$LASSO)),
          nrow(hub))) %>%
    mutate(combo = factor(combo, rev(combo)), all3 = combo == "All three")
  p_venn <- ggplot(combos, aes(n, combo)) +
    geom_col(aes(fill = all3), width = 0.7) +
    geom_text(aes(label = n), hjust = -0.3, size = 2.4, colour = ink2) +
    scale_fill_manual(values = c(`TRUE` = "#e34948", `FALSE` = "#c3c2b7"), guide = "none") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(x = "Genes", y = NULL, title = "Overlap of selected genes") +
    theme_pub(8)

  p_fig3a <- (p_cv | p_lasso) / (p_rf | p_venn)
  save_fig(p_fig3a, "Fig3A_ML_selection", 180, 150)

  # 6e 每个核心基因在各数据集的 AUC
  gplot <- gene_auc %>%
    mutate(dataset = factor(dataset, rev(c(train_ids, val_ids))),
           gene = factor(gene, hub$gene))
  p_auc <- ggplot(gplot, aes(AUC, dataset)) +
    geom_vline(xintercept = 0.5, colour = "#c3c2b7", linewidth = 0.3) +
    geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.4, colour = ink2, orientation = "y") +
    geom_point(aes(shape = set), size = 1.6, colour = "#0b0b0b") +
    scale_shape_manual(values = c(Training = 16, `External validation` = 17), name = NULL) +
    scale_x_continuous(limits = c(0.3, 1), breaks = c(0.5, 0.75, 1)) +
    facet_wrap(~ gene, ncol = min(5, nrow(hub))) +
    labs(x = "AUC (95% CI, DeLong)", y = NULL) +
    theme_pub(8) + theme(legend.position = "bottom")
  save_fig(p_auc, "Fig3B_hub_gene_AUC", 180, 40 + 35 * ceiling(nrow(hub) / 5))

  # 6f 核心基因评分 ROC 曲线
  roc_df <- map_dfr(names(sets), function(i) {
    s <- sets[[i]]
    r <- roc(s$group, hub_score(s$z), levels = c("Control", "NP"), direction = "<", quiet = TRUE)
    tibble(dataset = i, fpr = 1 - r$specificities, tpr = r$sensitivities)
  }) %>% left_join(score_auc, by = "dataset") %>%
    mutate(lab = sprintf("%s  AUC %.2f (%.2f-%.2f)", dataset, AUC, lo, hi),
           set = ifelse(dataset %in% train_ids, "Training", "External validation"))
  p_roc <- ggplot(roc_df %>% arrange(dataset, fpr, tpr), aes(fpr, tpr)) +
    geom_abline(colour = "#c3c2b7", linewidth = 0.3, linetype = 2) +
    geom_step(colour = "#2a78d6", linewidth = 0.5) +
    facet_wrap(~ factor(lab, unique(lab[order(match(dataset, c(train_ids, val_ids)))])), ncol = 5) +
    coord_equal() +
    labs(x = "1 - Specificity", y = "Sensitivity") +
    theme_pub(7)
  save_fig(p_roc, "Fig3C_hub_score_ROC", 190, 60)
}

cat("\n完成\n")
