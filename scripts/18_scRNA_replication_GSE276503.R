# 18_scRNA_replication_GSE276503.R
# 在独立单细胞数据集 GSE276503（Xie 2024 Nat Commun；BGI DNBelab C4 平台）中复现成纤维细胞结果
# 样本：健康对照下鼻甲 2、CRSwNP 下鼻甲 9、CRSwNP 中鼻甲 7、CRSwNP 鼻息肉 15（约 25 万个细胞）
# 电脑内存有限，分两步：
#   1) 逐个样本读入、质控、按标志基因均值给每个细胞粗分大类；保存每个样本×大类的 counts 求和（伪 bulk），
#      只保留基质候选细胞（成纤维/周细胞/内皮/神经胶质）用于下一步
#   2) 合并基质候选细胞 → Harmony → 聚类 → 精细注释（与 15_scRNA_revision.R 同一套标志基因）
#   3) 成纤维细胞伪 bulk（counts 求和 + edgeR QL）：息肉 vs CRSwNP 下鼻甲 / 中鼻甲 / 健康下鼻甲
#   4) 把 GSE235711 中 SCG2 主导亚群（F1）的标志基因（去掉即刻早期基因）投射到本数据集，
#      看该状态是否在息肉成纤维细胞中富集；并检查解离应激（即刻早期基因评分）是否能解释 SCG2 的差异
#   5) SCG2 与 COL10A1 共表达
#
# 输出：results/tables/rev_sc2_*.csv、results/figures/FigR_sc2_*

suppressPackageStartupMessages({
  library(tidyverse); library(Seurat); library(harmony); library(edgeR); library(Matrix); library(patchwork)
})
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
raw <- "data/raw/GSE276503"; tab <- "results/tables"; fig <- "results/figures"
set.seed(2026)
options(future.globals.maxSize = 8 * 1024^3)
genes_focus <- c("SCG2", "COL10A1", "PTHLH", "PDE7A")
grp_levels <- c("Healthy IT", "CRSwNP IT", "CRSwNP MT", "NP")

broad <- list(
  Epithelial = c("EPCAM", "KRT5", "KRT8", "KRT19", "SCGB1A1", "FOXJ1"),
  Stromal = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "RGS5", "NOTCH3"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP"),
  Glandular = c("LTF", "BPIFB2", "ZG16B", "PIP"),
  `T/NK` = c("CD3E", "CD3D", "NKG7", "GNLY"),
  `B/Plasma` = c("MS4A1", "CD79A", "JCHAIN", "MZB1"),
  Myeloid = c("CD68", "CD14", "LYZ", "C1QA"),
  Mast = c("TPSAB1", "CPA3", "KIT"),
  Granulocyte = c("CLC", "CSF3R", "FCGR3B"),
  Glia = c("PLP1", "S100B", "SOX10", "CDH19"))

# ---------- 1. 逐个样本处理 ----------
smp <- tibble(f = list.files(raw, "_barcodes\\.tsv\\.gz$")) %>%
  mutate(gsm = sub("_.*", "", f),
         title = sub("^GSM[0-9]+_(.*)_barcodes\\.tsv\\.gz$", "\\1", f),
         title = ifelse(title == f, NA, title))
meta_gsm <- local({
  l <- readLines("data/raw/meta/GSE276503.gsm.txt", warn = FALSE)
  tibble(gsm = sub("^\\^SAMPLE = ", "", grep("^\\^SAMPLE", l, value = TRUE)),
         sample = sub("^!Sample_title = ", "", grep("^!Sample_title", l, value = TRUE)))
})
smp <- smp %>% left_join(meta_gsm, by = "gsm") %>%
  mutate(group = factor(case_when(grepl("^HC_IT", sample) ~ "Healthy IT", grepl("^CRSwNP_IT", sample) ~ "CRSwNP IT",
                                  grepl("^CRSwNP_MT", sample) ~ "CRSwNP MT", grepl("^CRSwNP_NP", sample) ~ "NP"), grp_levels))
print(count(smp, group))

read_one <- function(gsm) {
  fs <- list.files(raw, paste0("^", gsm, "_"), full.names = TRUE)
  m <- readMM(grep("matrix", fs, value = TRUE))
  hdr <- scan(gzfile(grep("matrix", fs, value = TRUE)), what = integer(), skip = 2, nmax = 3, quiet = TRUE)
  if (length(m@x) != hdr[3]) return(NULL)        # GEO 原始文件损坏（见下方说明），跳过
  g <- read.delim(grep("genes", fs, value = TRUE), header = FALSE)$V1
  b <- read.delim(grep("barcodes", fs, value = TRUE), header = FALSE)$V1
  m <- as(m, "CsparseMatrix"); rownames(m) <- make.unique(g); colnames(m) <- paste0(gsm, "_", b)
  m
}
qc_tab <- list(); pb_list <- list(); keep_list <- list(); excluded <- character(0)
for (i in seq_len(nrow(smp))) {
  m <- read_one(smp$gsm[i])
  # GSE276503 中 4 个样本（CRSwNP_IT11、IT14、IT15、MT22）的 GEO 原始矩阵文件内部有整段零字节，
  # 条目数少于文件头声明的数目；两次独立下载结果完全相同（文件大小与服务器一致），属于数据提交时的损坏，排除
  if (is.null(m)) { excluded <- c(excluded, smp$sample[i]); cat(smp$sample[i], ": 原始文件损坏，排除\n"); next }
  nf <- Matrix::colSums(m > 0); mt <- Matrix::colSums(m[grep("^MT-", rownames(m)), , drop = FALSE]) / Matrix::colSums(m) * 100
  ok <- nf > 200 & nf < 7500 & mt < 20
  m <- m[, ok]
  ln <- log1p(t(t(m) / Matrix::colSums(m)) * 1e4)
  sc_mat <- sapply(broad, function(g) Matrix::colMeans(ln[intersect(g, rownames(ln)), , drop = FALSE]))
  lab <- colnames(sc_mat)[max.col(sc_mat, ties.method = "first")]
  lab[apply(sc_mat, 1, max) == 0] <- "Unassigned"
  qc_tab[[i]] <- tibble(sample = smp$sample[i], group = smp$group[i], cells_raw = length(ok), cells_qc = sum(ok))
  pb_list[[i]] <- map(setNames(nm = setdiff(unique(lab), "Unassigned")), ~ Matrix::rowSums(m[, lab == .x, drop = FALSE]))
  keep <- lab %in% c("Stromal", "Endothelial", "Glia")
  s <- CreateSeuratObject(m[, keep], min.features = 0)
  s$sample <- smp$sample[i]; s$group <- smp$group[i]; s$broad <- lab[keep]
  keep_list[[i]] <- s
  cat(smp$sample[i], ": QC", sum(ok), "个细胞；基质候选", sum(keep), "\n")
  rm(m, ln); gc(verbose = FALSE)
}
qc_tab <- bind_rows(qc_tab); write_csv(qc_tab, file.path(tab, "rev_sc2_QC.csv"))
writeLines(excluded, file.path(tab, "rev_sc2_excluded_corrupt_samples.txt"))
used <- !smp$sample %in% excluded
pb_list <- pb_list[used]; keep_list <- compact(keep_list); smp <- smp[used, ]
cat("质控后共", sum(qc_tab$cells_qc), "个细胞；排除损坏样本", length(excluded), "个；纳入：\n"); print(count(smp, group))

# ---------- 2. 基质细胞整合与精细注释 ----------
st <- merge(keep_list[[1]], keep_list[-1]); rm(keep_list); gc()
st <- JoinLayers(st)
st <- NormalizeData(st, verbose = FALSE) %>% FindVariableFeatures(nfeatures = 2000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>% RunPCA(npcs = 20, verbose = FALSE)
st <- RunHarmony(st, "sample", verbose = FALSE)
st <- RunUMAP(st, reduction = "harmony", dims = 1:20, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:20, verbose = FALSE) %>% FindClusters(resolution = 0.8, verbose = FALSE)
st_markers <- list(
  Fibroblast = c("DCN", "LUM", "PDGFRA", "COL1A2", "MFAP5", "FBLN1"),
  `Pericyte/SMC` = c("RGS5", "NOTCH3", "KCNJ8", "MYH11", "ACTA2", "MCAM", "HIGD1B"),
  `Schwann/glia` = c("PLP1", "S100B", "SOX10", "MPZ", "NRXN1", "CDH19"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP", "CLDN5"),
  Lymphatic = c("PROX1", "LYVE1", "CCL21", "TFF3"),
  Contaminant = c("PTPRC", "EPCAM", "KRT5", "LYZ", "CD3E"))
st <- AddModuleScore(st, map(st_markers, ~ intersect(.x, rownames(st))), name = "st_", seed = 2026)
cl <- st@meta.data %>% group_by(seurat_clusters) %>% summarise(across(starts_with("st_"), mean), n = n(), .groups = "drop")
cl$label <- names(st_markers)[max.col(as.matrix(cl[, paste0("st_", seq_along(st_markers))]))]
st$sub_type <- cl$label[match(st$seurat_clusters, cl$seurat_clusters)]
write_csv(cl, file.path(tab, "rev_sc2_stromal_clusters.csv"))
cat("\n基质细胞精细注释 × 分组：\n"); print(table(st$sub_type, st$group))

ln <- GetAssayData(st, layer = "data")
expr_type <- map_dfr(intersect(c(genes_focus, "CHGA", "RGS5"), rownames(st)), function(g)
  tibble(gene = g, cell_type = st$sub_type, v = ln[g, ]) %>% group_by(gene, cell_type) %>%
    summarise(n = n(), pct = mean(v > 0) * 100, .groups = "drop"))
write_csv(expr_type, file.path(tab, "rev_sc2_expression_stromal_types.csv"))
print(expr_type %>% filter(gene %in% c("SCG2", "COL10A1")) %>% pivot_wider(names_from = gene, values_from = pct))
# 所有细胞大类中的检测比例（来自伪 bulk 的 counts 无法得到比例，这里只报告基质内细分；
# 非基质细胞中 SCG2/COL10A1 的总量见 rev_sc2_pseudobulk_alltypes.csv）

# ---------- 3. 伪 bulk：成纤维细胞（精细注释）+ 其他大类（粗注释） ----------
cnt <- GetAssayData(st, layer = "counts")
fib_pb <- sapply(smp$sample, function(s) {
  idx <- which(st$sample == s & st$sub_type == "Fibroblast")
  if (length(idx) < 20) return(rep(NA_real_, nrow(cnt)))
  Matrix::rowSums(cnt[, idx, drop = FALSE])
})
rownames(fib_pb) <- rownames(cnt)
n_fib <- table(factor(st$sample[st$sub_type == "Fibroblast"], smp$sample))
ok <- colSums(is.na(fib_pb)) == 0
cat("\n纳入伪 bulk 的样本（成纤维细胞 ≥ 20 个）：\n"); print(table(smp$group[ok]))
grp <- droplevels(smp$group[ok])
y <- DGEList(fib_pb[, ok], group = grp); y <- normLibSizes(y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE])
des <- model.matrix(~ 0 + grp); colnames(des) <- make.names(levels(grp))
fit <- glmQLFit(estimateDisp(y, des), des, robust = TRUE)
lc <- cpm(y, log = TRUE, prior.count = 1)
cons <- c(`NP vs CRSwNP IT` = "NP - CRSwNP.IT", `NP vs CRSwNP MT` = "NP - CRSwNP.MT",
          `NP vs healthy IT` = "NP - Healthy.IT", `CRSwNP IT vs healthy IT` = "CRSwNP.IT - Healthy.IT")
pb <- map_dfr(names(cons), function(nm) {
  parts <- strsplit(cons[[nm]], " - ")[[1]]
  if (!all(parts %in% colnames(des))) return(NULL)
  tt <- topTags(glmQLFTest(fit, contrast = makeContrasts(contrasts = cons[[nm]], levels = des)), n = Inf)$table
  g <- intersect(genes_focus, rownames(tt))
  a1 <- make.names(as.character(grp)) == parts[1]; a0 <- make.names(as.character(grp)) == parts[2]
  tibble(comparison = nm, gene = g, n1 = sum(a1), n0 = sum(a0), logFC = tt[g, "logFC"], P_QL = tt[g, "PValue"],
         FDR_genomewide = tt[g, "FDR"],
         P_exact_wilcoxon = map_dbl(g, ~ wilcox.test(lc[.x, a1], lc[.x, a0], exact = TRUE)$p.value))
})
write_csv(pb, file.path(tab, "rev_sc2_fibroblast_pseudobulk.csv"))
cat("\n成纤维细胞伪 bulk（edgeR QL）：\n"); print(pb %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 30)

# 粗注释大类的伪 bulk（看 SCG2/COL10A1/PTHLH 在哪类细胞中升高）
all_pb <- map_dfr(seq_len(nrow(smp)), function(i) imap_dfr(pb_list[[i]], function(v, ct)
  tibble(sample = smp$sample[i], group = smp$group[i], cell_type = ct, lib = sum(v),
         !!!setNames(as.list(v[intersect(genes_focus, names(v))]), intersect(genes_focus, names(v))))))
write_csv(all_pb, file.path(tab, "rev_sc2_pseudobulk_alltypes.csv"))
share <- all_pb %>% group_by(cell_type) %>% summarise(across(any_of(genes_focus), sum)) %>%
  mutate(across(any_of(genes_focus), ~ .x / sum(.x) * 100))
cat("\n各大类细胞贡献的 UMI 比例（%，全部样本）：\n"); print(share %>% mutate(across(where(is.numeric), ~ round(.x, 1))))

# ---------- 4. 投射 GSE235711 的 F1 状态；即刻早期基因（解离应激）检查 ----------
f1_mk <- read_csv(file.path(tab, "rev_sc1_F1_markers.csv"), show_col_types = FALSE)
ieg <- c("FOS", "FOSB", "JUN", "JUNB", "JUND", "ATF3", "EGR1", "EGR2", "EGR3", "IER2", "IER3", "ZFP36", "NR4A1",
         "NR4A2", "NR4A3", "DUSP1", "KLF2", "KLF4", "KLF6", "GADD45B", "PPP1R15A", "CCNL1", "SOCS3", "HSPA1A", "HSPA1B")
f1_sig <- f1_mk %>% filter(!gene %in% c(ieg, "SCG2", "COL10A1"), p_val_adj < 0.05) %>% slice_max(avg_log2FC, n = 30) %>% pull(gene)
cat("\nF1 状态标志（去除即刻早期基因及 SCG2/COL10A1，", length(f1_sig), "个）：", paste(f1_sig, collapse = ", "), "\n")
fb <- subset(st, sub_type == "Fibroblast")
fb <- AddModuleScore(fb, list(intersect(f1_sig, rownames(fb)), intersect(ieg, rownames(fb))), name = c("sig_"), seed = 2026)
fb$F1_score <- fb$sig_1; fb$IEG_score <- fb$sig_2
fl <- GetAssayData(fb, layer = "data")
fb$SCG2 <- fl["SCG2", ]; fb$COL10A1 <- if ("COL10A1" %in% rownames(fl)) fl["COL10A1", ] else 0
per_s <- fb@meta.data %>% group_by(sample, group) %>%
  summarise(n_fib = n(), F1_score = mean(F1_score), IEG_score = mean(IEG_score),
            pct_SCG2 = mean(SCG2 > 0) * 100, pct_COL10A1 = mean(COL10A1 > 0) * 100, .groups = "drop")
write_csv(per_s, file.path(tab, "rev_sc2_fibroblast_per_sample.csv"))
cat("\n每样本成纤维细胞：\n"); print(per_s %>% arrange(group) %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = 40)
tests <- map_dfr(list(c("NP", "CRSwNP IT"), c("NP", "CRSwNP MT"), c("NP", "Healthy IT")), function(p) {
  a <- per_s %>% filter(group == p[1], n_fib >= 20); b <- per_s %>% filter(group == p[2], n_fib >= 20)
  map_dfr(c("F1_score", "IEG_score", "pct_SCG2", "pct_COL10A1"), ~ tibble(
    comparison = paste(p, collapse = " vs "), metric = .x, median_1 = median(a[[.x]]), median_0 = median(b[[.x]]),
    n1 = nrow(a), n0 = nrow(b), P_exact = wilcox.test(a[[.x]], b[[.x]], exact = TRUE)$p.value))
})
write_csv(tests, file.path(tab, "rev_sc2_fibroblast_state_tests.csv"))
print(tests %>% mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 20)
# 细胞水平：SCG2 与 F1 评分、IEG 评分的关系（混合效应近似：样本内中心化后 Spearman）
cc <- fb@meta.data %>% group_by(sample) %>%
  mutate(across(c(SCG2, F1_score, IEG_score), ~ .x - mean(.x))) %>% ungroup()
cat("\n样本内中心化后 Spearman：SCG2~F1", round(cor(cc$SCG2, cc$F1_score, method = "spearman"), 3),
    "；SCG2~IEG", round(cor(cc$SCG2, cc$IEG_score, method = "spearman"), 3), "\n")
# 伪 bulk 中把 IEG 评分作为协变量，看 SCG2 的组间差异是否保留
cov_ieg <- per_s$IEG_score[match(colnames(y), per_s$sample)]
des2 <- model.matrix(~ 0 + grp + cov_ieg); colnames(des2) <- c(make.names(levels(grp)), "IEG")
fit2 <- glmQLFit(estimateDisp(y, des2), des2, robust = TRUE)
pb_ieg <- map_dfr(names(cons), function(nm) {
  parts <- strsplit(cons[[nm]], " - ")[[1]]; if (!all(parts %in% colnames(des2))) return(NULL)
  tt <- topTags(glmQLFTest(fit2, contrast = makeContrasts(contrasts = cons[[nm]], levels = des2)), n = Inf)$table
  g <- intersect(genes_focus, rownames(tt))
  tibble(comparison = nm, gene = g, logFC_IEGadj = tt[g, "logFC"], P_IEGadj = tt[g, "PValue"], FDR_IEGadj = tt[g, "FDR"])
})
write_csv(pb_ieg, file.path(tab, "rev_sc2_fibroblast_pseudobulk_IEGadjusted.csv"))
cat("\n校正 IEG 评分后：\n"); print(pb_ieg %>% mutate(across(where(is.numeric), ~ signif(.x, 3))))

# ---------- 5. 共表达 ----------
co <- fb@meta.data %>% mutate(S = SCG2 > 0, C = COL10A1 > 0) %>% group_by(sample, group) %>%
  summarise(both = sum(S & C), S_only = sum(S & !C), C_only = sum(!S & C), none = sum(!S & !C), .groups = "drop")
arr <- array(unlist(map(seq_len(nrow(co)), ~ matrix(c(co$both[.x], co$S_only[.x], co$C_only[.x], co$none[.x]), 2) + 0.5)),
             dim = c(2, 2, nrow(co)))
mh <- mantelhaen.test(arr)
cat("\n共表达 MH OR =", round(mh$estimate, 2), "95% CI", round(mh$conf.int, 2), "P =", signif(mh$p.value, 3), "\n")
write_csv(co %>% mutate(MH_OR = mh$estimate, MH_lo = mh$conf.int[1], MH_hi = mh$conf.int[2], MH_P = mh$p.value),
          file.path(tab, "rev_sc2_coexpression.csv"))

saveRDS(st, "data/processed/scRNA_GSE276503_stromal.rds")
cat("\n完成\n")
