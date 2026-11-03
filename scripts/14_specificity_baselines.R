# 14_specificity_baselines.R
# 回应审稿意见（2026-11-01）的 bulk 部分补充分析
#   A. 荟萃分析模型敏感性：REML + Hartung-Knapp（metafor），与 DerSimonian-Laird 结果比较
#   B. AUC 基线：4 个共识基因的评分 vs
#        (1) 从训练集上调候选基因中随机抽 4 个（1,000 次）
#        (2) 不用机器学习、直接取训练集 logFC 最大的 4 个上调基因
#        (3) 经典息肉标志物（POSTN、CLC、CCL26、CST1、ALOX15）
#   C. 疾病特异性：共识基因在"非息肉的慢性鼻窦炎黏膜"中是否也升高
#        GSE36830：CRSsNP 钩突、CRSwNP 患者钩突 vs 对照钩突；息肉 vs 同组患者钩突
#        GSE23552：哮喘型 CRSwNP 患者的炎症黏膜（非息肉）、过敏性鼻炎黏膜 vs 正常黏膜
#   D. 解剖部位：GSE136825 中同一编号患者的息肉 vs 下鼻甲（配对 limma）；
#      CRSwNP 患者下鼻甲 vs 对照下鼻甲（同一部位，只差疾病）
#
# 输出：results/tables/rev_*.csv、results/figures/FigR_*

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(limma)
  library(metafor)
  library(pROC)
  library(GEOquery)
  library(AnnotationDbi)
  library(patchwork)
})
select <- dplyr::select; rename <- dplyr::rename; first <- dplyr::first
slice <- dplyr::slice; desc <- dplyr::desc; reduce <- purrr::reduce

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
proc <- "data/processed"; raw <- "data/raw"
tab  <- "results/tables"; fig <- "results/figures"
set.seed(2026)

ds_ids    <- c("GSE136825", "GSE179265", "GSE36830", "GSE23552", "GSE194282")
train_ids <- c("GSE136825", "GSE179265")
val_ids   <- c("GSE36830", "GSE23552", "GSE194282")
ml  <- read_csv(file.path(tab, "ML_feature_selection.csv"), show_col_types = FALSE)
hub_genes <- ml$gene[ml$hub]
zrow <- function(m) t(scale(t(m)))

# ---------- A. REML + Hartung-Knapp ----------
long <- map_dfr(ds_ids, ~ read_csv(file.path(tab, paste0("DEG_", .x, ".csv")), show_col_types = FALSE) %>%
                  mutate(dataset = .x))
meta_dl <- read_csv(file.path(tab, "meta_all_genes.csv"), show_col_types = FALSE)
fit_reml <- function(g, v) {
  r <- tryCatch(rma(yi = g, vi = v, method = "REML", test = "knha", control = list(maxiter = 1000)),
                error = function(e) NULL)
  if (is.null(r)) return(tibble(g_reml = NA_real_, lo_reml = NA_real_, hi_reml = NA_real_,
                                P_reml = NA_real_, tau2_reml = NA_real_))
  tibble(g_reml = as.numeric(r$b), lo_reml = r$ci.lb, hi_reml = r$ci.ub, P_reml = r$pval, tau2_reml = r$tau2)
}
meta_reml <- long %>% filter(gene %in% meta_dl$gene) %>%
  group_by(gene) %>% summarise(res = list(fit_reml(g, var_g)), .groups = "drop") %>% unnest(res) %>%
  mutate(FDR_reml = p.adjust(P_reml, "BH"))
cmp <- meta_dl %>% select(gene, k, g_pooled, ci_low, ci_high, FDR, consistent, direction) %>%
  left_join(meta_reml, by = "gene") %>%
  mutate(robust_reml = FDR_reml < 0.05 & abs(g_reml) >= 0.8 & consistent)
write_csv(cmp, file.path(tab, "rev_A_meta_REML_HK.csv"))
nr <- sum(cmp$direction != "NS"); nb <- sum(cmp$direction != "NS" & cmp$robust_reml, na.rm = TRUE)
cat("A. 稳健差异基因：DL", nr, "个；其中在 REML+HK 下仍满足标准", nb, "个；REML+HK 共",
    sum(cmp$robust_reml, na.rm = TRUE), "个\n")
print(cmp %>% filter(gene %in% hub_genes) %>%
        select(gene, g_pooled, ci_low, ci_high, FDR, g_reml, lo_reml, hi_reml, FDR_reml) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))))

# ---------- 验证集表达（z 标准化，与 05 一致） ----------
val <- map(setNames(val_ids, val_ids), function(i) {
  d <- readRDS(file.path(proc, paste0(i, ".rds"))); list(z = zrow(d$expr), group = d$pheno$group)
})
score_auc <- function(genes, s) {
  g <- intersect(genes, rownames(s$z)); if (!length(g)) return(NA_real_)
  as.numeric(auc(roc(s$group, colMeans(s$z[g, , drop = FALSE]), levels = c("Control", "NP"),
                     direction = "<", quiet = TRUE)))
}

# ---------- B. AUC 基线 ----------
up_pool <- ml %>% filter(direction == "Up") %>% pull(gene)
up_pool <- up_pool[map_lgl(up_pool, function(g) all(map_lgl(val, ~ g %in% rownames(.x$z))))]
cat("\nB. 可用于随机抽样的训练集上调候选基因（3 个验证集均检测到）：", length(up_pool), "\n")
rand <- map_dfr(1:1000, function(b) {
  g <- sample(up_pool, 4)
  tibble(iter = b, dataset = val_ids, AUC = map_dbl(val, ~ score_auc(g, .x)))
})
top_fc <- ml %>% filter(direction == "Up", gene %in% up_pool) %>%
  mutate(mfc = (logFC.a + logFC.b) / 2) %>% slice_max(mfc, n = 4) %>% pull(gene)
classic <- c("POSTN", "CLC", "CCL26", "CST1", "ALOX15")
sets_ref <- c(list(`Consensus 4 genes` = hub_genes, `Top-4 training logFC (no ML)` = top_fc,
                   `Classic 5-gene score` = classic), setNames(as.list(classic), classic))
ref_auc <- imap_dfr(sets_ref, function(g, nm)
  tibble(set = nm, genes = paste(intersect(g, Reduce(union, map(val, ~ rownames(.x$z)))), collapse = ", "),
         dataset = val_ids, AUC = map_dbl(val, ~ score_auc(g, .x))))
rand_sum <- rand %>% group_by(dataset) %>%
  summarise(rand_median = median(AUC), rand_q05 = quantile(AUC, 0.05), rand_q95 = quantile(AUC, 0.95),
            rand_frac_ge_0.9 = mean(AUC >= 0.9), .groups = "drop")
hub_pct <- ref_auc %>% filter(set == "Consensus 4 genes") %>% select(dataset, AUC) %>%
  left_join(rand %>% rename(rAUC = AUC), by = "dataset") %>%
  group_by(dataset) %>% summarise(hub_AUC = first(AUC), pct_random_ge_hub = mean(rAUC >= first(AUC)), .groups = "drop")
rand_mean <- rand %>% group_by(iter) %>% summarise(mAUC = mean(AUC))
hub_mean <- mean(ref_auc$AUC[ref_auc$set == "Consensus 4 genes"])
cat("   训练集 logFC 最大的 4 个上调基因：", paste(top_fc, collapse = ", "), "\n")
print(ref_auc %>% mutate(AUC = round(AUC, 3)) %>% pivot_wider(names_from = dataset, values_from = AUC))
print(rand_sum %>% left_join(hub_pct, by = "dataset") %>% mutate(across(where(is.numeric), ~ round(.x, 3))))
cat("   三个验证集平均 AUC：共识基因", round(hub_mean, 3), "；随机 4 基因平均 AUC >= 共识基因的比例",
    mean(rand_mean$mAUC >= hub_mean), "\n")
write_csv(rand, file.path(tab, "rev_B_random4_AUC.csv"))
write_csv(ref_auc, file.path(tab, "rev_B_reference_AUC.csv"))
write_csv(rand_sum %>% left_join(hub_pct, by = "dataset"), file.path(tab, "rev_B_random4_summary.csv"))

# ---------- C. 疾病特异性：读取全部芯片样本 ----------
read_all <- function(id, pkg, grp_fun) {
  es <- getGEO(filename = file.path(raw, paste0(id, "_series_matrix.txt.gz")), getGPL = FALSE)
  t  <- pData(es)$title; grp <- grp_fun(t); keep <- !is.na(grp)
  ex <- exprs(es)[, keep]
  qx <- quantile(ex, c(0.99, 1, 0, 0.25), na.rm = TRUE)
  if (qx[1] > 100 || (qx[2] - qx[3] > 50 && qx[4] > 0)) { ex[ex <= 0] <- NA; ex <- log2(ex) }
  ex <- ex[complete.cases(ex), ]
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
  sym <- limma::alias2SymbolTable(mapIds(get(pkg), rownames(ex), "SYMBOL", "PROBEID", multiVals = "first"))
  ok <- !is.na(sym); ex <- ex[ok, ]; sym <- sym[ok]
  o <- order(rowMeans(ex), decreasing = TRUE); ex <- ex[o, ]; sym <- sym[o]
  ex <- ex[!duplicated(sym), ]; rownames(ex) <- sym[!duplicated(sym)]
  list(expr = ex, group = grp[keep], title = t[keep])
}
g36 <- read_all("GSE36830", "hgu133plus2.db", function(t) case_when(
  grepl("^Control UT", t) ~ "Control UT", grepl("^CRSsNP UT", t) ~ "CRSsNP UT",
  grepl("^CRSwNP UT", t) ~ "CRSwNP UT", grepl("^CRSwNP polyp", t) ~ "CRSwNP polyp"))
g23 <- read_all("GSE23552", "huex10sttranscriptcluster.db", function(t) case_when(
  grepl("^NasalMucosa_Normal", t) ~ "Normal mucosa", grepl("^NasalMucosa_AR", t) ~ "Allergic rhinitis",
  grepl("^NasalMucosa\\(inflamed\\)_aCRSm", t) & t != "NasalMucosa(inflamed)_aCRSm2" ~ "CRSwNP inflamed mucosa",
  grepl("^NasalPolyp_aCRSp", t) ~ "CRSwNP polyp"))
genes_c <- c(hub_genes, "POSTN", "CLC", "CCL26")
contrast_tab <- function(d, id, contrasts) {
  f <- factor(make.names(d$group)); des <- model.matrix(~ 0 + f); colnames(des) <- levels(f)
  fit <- lmFit(d$expr, des)
  map_dfr(names(contrasts), function(nm) {
    cm <- makeContrasts(contrasts = contrasts[[nm]], levels = des)
    tt <- topTable(eBayes(contrasts.fit(fit, cm)), number = Inf, sort.by = "none")
    g <- intersect(genes_c, rownames(tt))
    tibble(dataset = id, comparison = nm, gene = g, logFC = tt[g, "logFC"], P = tt[g, "P.Value"],
           n1 = sum(make.names(d$group) == sub(" - .*", "", contrasts[[nm]])),
           n0 = sum(make.names(d$group) == sub(".* - ", "", contrasts[[nm]])))
  })
}
resC <- bind_rows(
  contrast_tab(g36, "GSE36830", c(`CRSsNP UT vs control UT` = "CRSsNP.UT - Control.UT",
                                  `CRSwNP UT vs control UT` = "CRSwNP.UT - Control.UT",
                                  `Polyp vs control UT` = "CRSwNP.polyp - Control.UT",
                                  `Polyp vs CRSwNP UT` = "CRSwNP.polyp - CRSwNP.UT")),
  contrast_tab(g23, "GSE23552", c(`CRSwNP inflamed mucosa vs normal` = "CRSwNP.inflamed.mucosa - Normal.mucosa",
                                  `Allergic rhinitis vs normal` = "Allergic.rhinitis - Normal.mucosa",
                                  `Polyp vs normal` = "CRSwNP.polyp - Normal.mucosa",
                                  `Polyp vs CRSwNP inflamed mucosa` = "CRSwNP.polyp - CRSwNP.inflamed.mucosa"))) %>%
  group_by(dataset, comparison) %>% mutate(FDR_within = p.adjust(P, "BH")) %>% ungroup()
cat("\nC. 疾病特异性（logFC / P）：\n")
print(resC %>% mutate(v = sprintf("%.2f (%.2g)", logFC, P)) %>% select(dataset, comparison, gene, v) %>%
        pivot_wider(names_from = gene, values_from = v), width = 250)
write_csv(resC, file.path(tab, "rev_C_specificity.csv"))

# ---------- D. GSE136825：同一患者息肉 vs 下鼻甲；患者下鼻甲 vs 对照下鼻甲 ----------
x <- read.delim(file.path(raw, "GSE136825_counts.txt.gz"), check.names = FALSE)
cnt <- as.matrix(x[, -(1:6)]); colnames(cnt) <- sub("\\.bam$", "", colnames(cnt))
gsm <- local({
  l <- readLines(file.path(raw, "meta/GSE136825.gsm.txt"), warn = FALSE)
  tibble(title = sub("^!Sample_title = ", "", grep("^!Sample_title", l, value = TRUE)),
         desc  = sub("^!Sample_description = ", "", grep("^!Sample_description", l, value = TRUE)))
}) %>% mutate(prefix = sub(" .*", "", title), num = sub(".* ", "", title)) %>% filter(desc %in% colnames(cnt))
ens <- sub("\\..*$", "", x$Geneid)
sym_tab <- local({ y <- read.delim(file.path(raw, "GSE179265_counts.txt.gz")); setNames(y$Gene_name, y$geneID) })
rownames(cnt) <- ens
hub_ens <- names(sym_tab)[sym_tab %in% genes_c]
voom_test <- function(cols, grp, block = NULL) {
  y <- DGEList(cnt[, cols]); keep <- filterByExpr(y, group = grp)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  des <- if (is.null(block)) model.matrix(~ grp) else model.matrix(~ block + grp)
  tt <- topTable(eBayes(lmFit(voom(y, des), des)), coef = ncol(des), number = Inf, sort.by = "none")
  e <- intersect(hub_ens, rownames(tt))
  tibble(gene = unname(sym_tab[e]), logFC = tt[e, "logFC"], P = tt[e, "P.Value"], FDR = tt[e, "adj.P.Val"])
}
pairs <- gsm %>% mutate(pair = case_when(prefix %in% c("PY", "SCP") ~ paste0("P", num),
                                         prefix %in% c("TT", "SCT") ~ paste0("T", num)),
                        tissue = case_when(prefix %in% c("PY", "TT") ~ "Polyp", prefix %in% c("SCP", "SCT") ~ "IT")) %>%
  filter(!is.na(pair)) %>% group_by(pair) %>% filter(n_distinct(tissue) == 2) %>% ungroup()
cat("\nD. GSE136825 按编号配对的息肉/下鼻甲：", n_distinct(pairs$pair), "对\n")
d1 <- voom_test(pairs$desc, factor(pairs$tissue, c("IT", "Polyp")), block = factor(pairs$pair)) %>%
  mutate(comparison = sprintf("Polyp vs own IT (paired, n=%d pairs)", n_distinct(pairs$pair)))
it <- gsm %>% filter(prefix %in% c("SCP", "SCT", "CC", "CCS", "CCY")) %>%
  mutate(grp = factor(ifelse(prefix %in% c("SCP", "SCT"), "CRSwNP_IT", "Control_IT"), c("Control_IT", "CRSwNP_IT")))
d2 <- voom_test(it$desc, it$grp) %>%
  mutate(comparison = sprintf("CRSwNP IT vs control IT (n=%d vs %d)", sum(it$grp == "CRSwNP_IT"), sum(it$grp == "Control_IT")))
resD <- bind_rows(d1, d2)
print(resD %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 30)
write_csv(resD, file.path(tab, "rev_D_site_paired_GSE136825.csv"))

# ---------- 作图：随机基因集 AUC 分布 ----------
p_rand <- ggplot(rand, aes(AUC)) +
  geom_histogram(bins = 40, fill = "#c3c2b7", colour = NA) +
  geom_vline(data = ref_auc %>% filter(set %in% c("Consensus 4 genes", "Top-4 training logFC (no ML)", "Classic 5-gene score")),
             aes(xintercept = AUC, colour = set, linetype = set), linewidth = 0.6) +
  scale_colour_manual(values = c("#e34948", "#2a78d6", "#1baf7a"), name = NULL) +
  scale_linetype_manual(values = c(1, 2, 4), name = NULL) +
  facet_wrap(~ factor(dataset, val_ids), nrow = 1) +
  labs(x = "AUC in withheld validation cohort", y = "Random 4-gene sets (n = 1,000)") +
  theme_classic(base_size = 8) + theme(legend.position = "bottom", strip.background = element_blank(),
                                       strip.text = element_text(face = "bold"))
ggsave(file.path(fig, "FigR_random_AUC.png"), p_rand, width = 180, height = 70, units = "mm", dpi = 300, bg = "white")
ggsave(file.path(fig, "FigR_random_AUC.pdf"), p_rand, width = 180, height = 70, units = "mm")
cat("\n完成\n")
