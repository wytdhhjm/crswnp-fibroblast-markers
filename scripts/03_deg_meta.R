# 03_deg_meta.R
# 第 1 步：各数据集内部差异分析 + 5 队列随机效应荟萃分析
#
# 思路（写进论文方法部分）：
#   1) RNA-seq（GSE136825、GSE179265）：edgeR 过滤低表达基因 + TMM 标准化 + limma-voom
#      芯片（GSE36830、GSE23552、GSE194282）：limma
#   2) 不同平台的 log2FC 尺度不可比（芯片信号会被压缩），因此每个数据集把 t 统计量
#      换算成标准化均数差 Hedges' g（效应量），再用 DerSimonian-Laird 随机效应模型合并
#      （与 MetaIntegrator 的做法一致）
#   3) "稳健差异基因"标准：
#      至少 4 个数据集检测到该基因（个别平台缺探针或低表达被过滤时不至于丢掉真实信号）；
#      荟萃 FDR < 0.05、|合并 g| >= 0.8（Cohen 定义的大效应）、所有检测到的数据集方向一致
#      注意：最初的标准是"5 个数据集都检测到"，看到 CLC、PRR4 等经典标志基因因单一平台
#      缺探针被排除后才改为"至少 4 个"——这是看过结果后的调整，论文中须如实说明；
#      原标准的结果作为敏感性分析一并输出（meta_sensitivity_k5.csv）
#
# 输出：results/tables/、results/figures/

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(ggrepel)
  library(patchwork)
})

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
proc <- "data/processed"
tab  <- "results/tables"
fig  <- "results/figures"
dir.create(tab, showWarnings = FALSE, recursive = TRUE)
dir.create(fig, showWarnings = FALSE, recursive = TRUE)

# 配色（参考调色板中已验证的色值）
col_group <- c(Control = "#2a78d6", NP = "#eb6834")       # 分组：分类色 1、2
col_dir   <- c(Down = "#2a78d6", Up = "#e34948", NS = "#c3c2b7")  # 上/下调：发散色两端 + 灰色中性
ink2 <- "#52514e"; grid_col <- "#e1e0d9"

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

ds_ids <- c("GSE136825", "GSE179265", "GSE36830", "GSE23552", "GSE194282")
ds <- setNames(lapply(ds_ids, function(i) readRDS(file.path(proc, paste0(i, ".rds")))), ds_ids)

# ---------- 1. 各数据集差异分析 ----------
run_limma <- function(d) {
  grp <- d$pheno$group
  design <- model.matrix(~ grp)
  if (!is.null(d$counts)) {                           # RNA-seq
    y <- DGEList(d$counts, group = grp)
    y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
    y <- calcNormFactors(y)
    v <- voom(y, design)
    expr <- v$E
    fit <- lmFit(v, design)
  } else {                                            # 芯片
    expr <- d$expr
    # 去掉低表达探针：至少在最小组样本数个样本中高于全体下四分位数（背景噪声水平）
    keep <- rowSums(expr > quantile(expr, 0.25)) >= min(table(grp))
    expr <- expr[keep, ]
    fit <- lmFit(expr, design)
  }
  fit <- eBayes(fit)
  tt <- topTable(fit, coef = 2, number = Inf, sort.by = "none")
  n1 <- sum(grp == "NP"); n0 <- sum(grp == "Control")
  # t -> Cohen's d -> Hedges' g 及其方差
  d_ <- tt$t * sqrt(1 / n1 + 1 / n0)
  J  <- 1 - 3 / (4 * (n1 + n0) - 9)
  g  <- J * d_
  vg <- 1 / n1 + 1 / n0 + g^2 / (2 * (n1 + n0))
  list(expr = expr,
       res = tibble(gene = rownames(tt), logFC = tt$logFC, t = tt$t,
                    P = tt$P.Value, FDR = tt$adj.P.Val, g = g, var_g = vg),
       n = c(NP = n1, Control = n0))
}

deg <- lapply(ds, run_limma)

for (i in ds_ids) {
  write_csv(deg[[i]]$res %>% arrange(P), file.path(tab, paste0("DEG_", i, ".csv")))
}

per_ds_summary <- map_dfr(ds_ids, function(i) {
  r <- deg[[i]]$res
  tibble(dataset = i, NP = deg[[i]]$n["NP"], Control = deg[[i]]$n["Control"],
         genes_tested = nrow(r),
         up = sum(r$FDR < 0.05 & r$logFC > 1), down = sum(r$FDR < 0.05 & r$logFC < -1))
})
cat("\n各数据集差异基因数（FDR < 0.05 且 |log2FC| > 1）：\n"); print(per_ds_summary)
write_csv(per_ds_summary, file.path(tab, "DEG_per_dataset_summary.csv"))

# ---------- 2. 随机效应荟萃分析（DerSimonian-Laird） ----------
long <- map_dfr(ds_ids, ~ deg[[.x]]$res %>% mutate(dataset = .x))
common <- long %>% count(gene) %>% filter(n >= length(ds_ids) - 1) %>% pull(gene)
cat("\n至少 4 个数据集检测到的基因数：", length(common), "\n")

meta <- long %>%
  filter(gene %in% common) %>%
  group_by(gene) %>%
  summarise(
    k = n(),
    w = list(1 / var_g), gg = list(g), vv = list(var_g),
    n_up = sum(g > 0), n_down = sum(g < 0),
    .groups = "drop") %>%
  mutate(res = pmap(list(gg, vv), function(g, v) {
    w  <- 1 / v
    gf <- sum(w * g) / sum(w)                          # 固定效应合并值
    Q  <- sum(w * (g - gf)^2)
    df <- length(g) - 1
    tau2 <- max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))
    ws <- 1 / (v + tau2)
    est <- sum(ws * g) / sum(ws)
    se  <- sqrt(1 / sum(ws))
    z   <- est / se
    tibble(g_pooled = est, se = se, ci_low = est - 1.96 * se, ci_high = est + 1.96 * se,
           z = z, P = 2 * pnorm(-abs(z)),
           tau2 = tau2, Q = Q, Q_P = pchisq(Q, df, lower.tail = FALSE),
           I2 = if (Q > 0) max(0, (Q - df) / Q) * 100 else 0)
  })) %>%
  select(gene, k, n_up, n_down, res) %>%
  unnest(res) %>%
  mutate(FDR = p.adjust(P, "BH"),
         consistent = n_up == k | n_down == k,
         direction = case_when(FDR < 0.05 & abs(g_pooled) >= 0.8 & consistent & g_pooled > 0 ~ "Up",
                               FDR < 0.05 & abs(g_pooled) >= 0.8 & consistent & g_pooled < 0 ~ "Down",
                               TRUE ~ "NS")) %>%
  arrange(P)

write_csv(meta, file.path(tab, "meta_all_genes.csv"))
robust <- meta %>% filter(direction != "NS")
write_csv(robust, file.path(tab, "meta_robust_DEGs.csv"))
cat("稳健差异基因：上调", sum(robust$direction == "Up"), "个，下调", sum(robust$direction == "Down"), "个\n")

# 稳健性检验：依次去掉一个数据集重做荟萃分析，看稳健基因是否仍显著
loo <- map_dfr(ds_ids, function(drop) {
  sub <- long %>% filter(gene %in% robust$gene, dataset != drop)
  sub %>% group_by(gene) %>% summarise(
    sig = {
      w <- 1 / var_g; gf <- sum(w * g) / sum(w); Q <- sum(w * (g - gf)^2); df <- n() - 1
      tau2 <- max(0, (Q - df) / (sum(w) - sum(w^2) / sum(w)))
      ws <- 1 / (var_g + tau2); est <- sum(ws * g) / sum(ws)
      2 * pnorm(-abs(est / sqrt(1 / sum(ws))))
    }, .groups = "drop") %>%
    mutate(dropped = drop)
}) %>% group_by(dropped) %>% mutate(FDR = p.adjust(sig, "BH")) %>% ungroup()
loo_sum <- loo %>% group_by(gene) %>% summarise(still_sig_all = all(FDR < 0.05), .groups = "drop")
cat("留一法检验：", sum(loo_sum$still_sig_all), "/", nrow(loo_sum), "个稳健基因在去掉任一数据集后仍显著\n")
write_csv(robust %>% left_join(loo_sum, by = "gene"), file.path(tab, "meta_robust_DEGs.csv"))

# 敏感性分析：原始标准（5 个数据集均检测到），FDR 在该子集内重新校正
sens <- meta %>% filter(k == length(ds_ids)) %>%
  mutate(FDR = p.adjust(P, "BH"),
         robust_k5 = FDR < 0.05 & abs(g_pooled) >= 0.8 & consistent)
write_csv(sens, file.path(tab, "meta_sensitivity_k5.csv"))
cat("敏感性分析（k = 5）：稳健基因", sum(sens$robust_k5), "个，其中",
    sum(sens$gene[sens$robust_k5] %in% robust$gene), "个也在主分析结果中\n")

# ---------- 3. 作图 ----------
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}

# 3a 各数据集 PCA（样本质控：看对照与息肉是否分开、有无离群样本）
pca_plots <- map(ds_ids, function(i) {
  e <- deg[[i]]$expr
  e <- e[order(apply(e, 1, var), decreasing = TRUE)[1:2000], ]
  pc <- prcomp(t(e), scale. = TRUE)
  ve <- round(100 * pc$sdev^2 / sum(pc$sdev^2), 1)
  df <- tibble(PC1 = pc$x[, 1], PC2 = pc$x[, 2], group = ds[[i]]$pheno$group)
  ggplot(df, aes(PC1, PC2, colour = group)) +
    geom_point(size = 1.8, alpha = 0.9) +
    scale_colour_manual(values = col_group, name = NULL) +
    labs(title = i, x = paste0("PC1 (", ve[1], "%)"), y = paste0("PC2 (", ve[2], "%)")) +
    theme_pub()
})
p_pca <- wrap_plots(pca_plots, ncol = 3) + plot_layout(guides = "collect") &
  theme(legend.position = "bottom")
save_fig(p_pca, "FigS1_PCA_per_dataset", 180, 120)

# 3b 荟萃火山图
lab <- bind_rows(robust %>% filter(direction == "Up") %>% slice_max(abs(g_pooled), n = 10),
                 robust %>% filter(direction == "Down") %>% slice_max(abs(g_pooled), n = 10))
p_vol <- ggplot(meta, aes(g_pooled, -log10(FDR))) +
  geom_point(data = filter(meta, direction == "NS"), colour = col_dir["NS"], size = 0.6, alpha = 0.5) +
  geom_point(data = filter(meta, direction != "NS"), aes(colour = direction), size = 1.1) +
  geom_vline(xintercept = c(-0.8, 0.8), linetype = 2, linewidth = 0.3, colour = "#898781") +
  geom_hline(yintercept = -log10(0.05), linetype = 2, linewidth = 0.3, colour = "#898781") +
  geom_text_repel(data = lab, aes(label = gene), size = 2.4, colour = "#0b0b0b",
                  max.overlaps = 30, segment.size = 0.2, segment.colour = "#898781", min.segment.length = 0) +
  scale_colour_manual(values = col_dir, name = NULL,
                      labels = c(Down = paste0("Down (n=", sum(robust$direction == "Down"), ")"),
                                 Up = paste0("Up (n=", sum(robust$direction == "Up"), ")"))) +
  labs(x = "Pooled effect size (Hedges' g)", y = expression(-log[10]~"(meta FDR)"),
       title = NULL) +
  theme_pub() + theme(legend.position = c(0.15, 0.9))
n_np <- sum(per_ds_summary$NP); n_ct <- sum(per_ds_summary$Control)
p_vol <- p_vol + labs(title = paste0("5-cohort meta-analysis (NP n=", n_np, ", control n=", n_ct, ")"))
save_fig(p_vol, "Fig1_meta_volcano", 120, 100)

# 3c 森林图：效应最大的上调/下调各 6 个稳健基因
top_forest <- c(lab$gene[lab$direction == "Up"][1:6], lab$gene[lab$direction == "Down"][1:6])
fdat <- long %>% filter(gene %in% top_forest) %>%
  mutate(lo = g - 1.96 * sqrt(var_g), hi = g + 1.96 * sqrt(var_g), type = "Cohort") %>%
  select(gene, dataset, g, lo, hi, type) %>%
  bind_rows(meta %>% filter(gene %in% top_forest) %>%
              transmute(gene, dataset = "Pooled", g = g_pooled, lo = ci_low, hi = ci_high, type = "Pooled")) %>%
  mutate(dataset = factor(dataset, rev(c(ds_ids, "Pooled"))),
         gene = factor(gene, top_forest))
p_forest <- ggplot(fdat, aes(g, dataset)) +
  geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = hi), width = 0, linewidth = 0.4, colour = ink2, orientation = "y") +
  geom_point(aes(shape = type, size = type), colour = "#0b0b0b", fill = "#eb6834") +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 23), name = NULL) +
  scale_size_manual(values = c(Cohort = 1.4, Pooled = 2.4), name = NULL) +
  facet_wrap(~ gene, ncol = 6, scales = "free_x") +
  labs(x = "Hedges' g (95% CI)", y = NULL) +
  theme_pub(8) + theme(legend.position = "bottom")
save_fig(p_forest, "Fig1_forest_top_genes", 180, 110)

cat("\n完成。表格见 results/tables/，图见 results/figures/\n")
