# 15_scRNA_revision.R
# 回应审稿意见：GSE235711 单细胞部分重做（全部 20 个样本）
#   1) 纳入全部分组：对照筛窦 5、CRSsNP 筛窦 5、CRSwNP 筛窦（非息肉黏膜）4、CRSwNP 鼻息肉 6
#   2) 基质细胞重新分群：在 成纤维/内皮/周细胞 中高分辨率再聚类，单独识别周细胞/平滑肌、
#      施万/神经胶质细胞，确认 SCG2 不是来自周细胞或神经来源
#   3) 伪 bulk 改为原始 counts 求和 + edgeR 准似然检验（Squair 2021 推荐做法）；
#      另用 log-CPM 的精确 Wilcoxon 检验做敏感性分析；报告每个细胞类型实际纳入的样本数
#   4) 成纤维细胞亚群：再聚类，并用已发表的鼻息肉成纤维细胞状态标志打分
#      （IDO1+ 成纤维细胞 Xie 2024；Cao 2026 的 IL-4 相关 C8 亚群；CXCL12+ 成纤维细胞 Guo 2024；
#       通用成纤维细胞图谱中的 PI16+ 与肌成纤维细胞状态）
#   5) SCG2 与 COL10A1 在同一成纤维细胞中的共表达（逐样本 2x2 表，Mantel-Haenszel 合并）
#
# 输出：data/processed/scRNA_GSE235711_all20.rds、results/tables/rev_sc1_*.csv、results/figures/FigR_sc1_*

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(edgeR)
  library(patchwork)
})

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
raw <- "data/raw/GSE235711"
tab <- "results/tables"; fig <- "results/figures"
set.seed(2026)
options(future.globals.maxSize = 8 * 1024^3)
genes_focus <- c("SCG2", "COL10A1", "PTHLH", "PDE7A")
grp_levels <- c("Control", "CRSsNP", "CRSwNP-Eth", "NP")
col_group <- c(Control = "#2a78d6", CRSsNP = "#1baf7a", `CRSwNP-Eth` = "#eda100", NP = "#eb6834")
theme_pub <- function(base = 8) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.ticks = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.text = element_text(colour = "#52514e"),
          plot.title = element_text(face = "bold", size = base + 1),
          legend.key.size = unit(3, "mm"),
          strip.background = element_blank(), strip.text = element_text(face = "bold"))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}
cluster_harmony <- function(s, nfeat = 3000, dims = 30, res = 0.5) {
  s <- NormalizeData(s, verbose = FALSE) %>% FindVariableFeatures(nfeatures = nfeat, verbose = FALSE) %>%
    ScaleData(verbose = FALSE) %>% RunPCA(npcs = dims, verbose = FALSE)
  s <- RunHarmony(s, "sample", verbose = FALSE)
  RunUMAP(s, reduction = "harmony", dims = 1:dims, verbose = FALSE) %>%
    FindNeighbors(reduction = "harmony", dims = 1:dims, verbose = FALSE) %>%
    FindClusters(resolution = res, verbose = FALSE)
}
annotate_by_score <- function(s, markers, prefix) {
  markers <- map(markers, ~ intersect(.x, rownames(s)))
  s <- AddModuleScore(s, markers, name = prefix, seed = 2026)
  cols <- paste0(prefix, seq_along(markers))
  cl <- s@meta.data %>% group_by(seurat_clusters) %>% summarise(across(all_of(cols), mean), .groups = "drop")
  lab <- names(markers)[max.col(as.matrix(cl[, cols]))]
  list(obj = s, label = lab[match(s$seurat_clusters, cl$seurat_clusters)],
       table = tibble(cluster = cl$seurat_clusters, label = lab))
}

# ---------- 1. 读取全部 20 个样本 ----------
files <- list.files(raw, pattern = "_filtered\\.tar\\.gz$", full.names = TRUE)
tmp <- file.path(tempdir(), "sc10x_all")
objs <- map(files, function(f) {
  id <- sub("_filtered\\.tar\\.gz$", "", basename(f))
  d <- file.path(tmp, id); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  untar(f, exdir = d)
  m <- suppressMessages(Read10X(d)); if (is.list(m)) m <- m[["Gene Expression"]]
  s <- CreateSeuratObject(m, project = id, min.cells = 3, min.features = 200)
  s$sample <- sub("^GSM[0-9]+_", "", id); s
})
sc <- merge(objs[[1]], objs[-1], add.cell.ids = map_chr(objs, ~ .x$sample[1])); rm(objs); gc()
sc$group <- factor(case_when(grepl("^Control", sc$sample) ~ "Control", grepl("^CRSsNP", sc$sample) ~ "CRSsNP",
                             grepl("^CRSwNP_Eth", sc$sample) ~ "CRSwNP-Eth", grepl("^CRSwNP_NP", sc$sample) ~ "NP"),
                   grp_levels)
sc[["mt"]] <- PercentageFeatureSet(sc, "^MT-")
sc <- subset(sc, nFeature_RNA > 200 & nFeature_RNA < 7500 & mt < 20)
sc <- JoinLayers(sc)
cat("质控后：", ncol(sc), "个细胞，", n_distinct(sc$sample), "个样本\n"); print(table(sc$group))
sc <- cluster_harmony(sc)

markers <- list(
  Epithelial = c("EPCAM", "KRT5", "KRT8", "KRT19", "SCGB1A1", "FOXJ1"),
  Stromal = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA", "RGS5", "NOTCH3", "PDGFRB"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP"),
  Glandular = c("LTF", "BPIFB2", "ZG16B", "PIP"),
  `T/NK` = c("CD3E", "CD3D", "TRAC", "NKG7", "GNLY"),
  `B/Plasma` = c("MS4A1", "CD79A", "JCHAIN", "MZB1"),
  Myeloid = c("CD68", "CD14", "LYZ", "C1QA", "CD1C"),
  Mast = c("TPSAB1", "CPA3", "KIT", "HDC"),
  Granulocyte = c("CLC", "CSF3R", "FCGR3B", "SIGLEC8"))
a <- annotate_by_score(sc, markers, "ct_"); sc <- a$obj; sc$cell_type <- a$label

# ---------- 2. 基质+内皮细胞再聚类：区分成纤维、周细胞/平滑肌、施万/胶质 ----------
st <- subset(sc, cell_type %in% c("Stromal", "Endothelial"))
st <- cluster_harmony(st, nfeat = 2000, dims = 20, res = 0.8)
st_markers <- list(
  Fibroblast = c("DCN", "LUM", "PDGFRA", "COL1A2", "MFAP5", "FBLN1"),
  `Pericyte/SMC` = c("RGS5", "NOTCH3", "KCNJ8", "MYH11", "ACTA2", "MCAM", "HIGD1B"),
  `Schwann/glia` = c("PLP1", "S100B", "SOX10", "MPZ", "NRXN1", "CDH19"),
  Endothelial = c("PECAM1", "VWF", "CDH5", "PLVAP", "CLDN5"),
  Lymphatic = c("PROX1", "LYVE1", "CCL21", "TFF3"))
b <- annotate_by_score(st, st_markers, "st_"); st <- b$obj; st$sub_type <- b$label
write_csv(b$table, file.path(tab, "rev_sc1_stromal_subclusters.csv"))
sc$cell_type[colnames(st)] <- st$sub_type
sc$cell_type[sc$cell_type == "Stromal"] <- "Fibroblast"     # 未进入再聚类的基质细胞（不应存在）
ct_levels <- c("Epithelial", "Glandular", "Fibroblast", "Pericyte/SMC", "Schwann/glia", "Endothelial",
               "Lymphatic", "T/NK", "B/Plasma", "Myeloid", "Mast", "Granulocyte")
sc$cell_type <- factor(sc$cell_type, intersect(ct_levels, unique(sc$cell_type)))
cnt_tab <- sc@meta.data %>% count(sample, group, cell_type) %>% complete(nesting(sample, group), cell_type, fill = list(n = 0))
write_csv(cnt_tab, file.path(tab, "rev_sc1_cells_per_sample_type.csv"))
cat("\n细胞类型 × 分组：\n"); print(table(sc$cell_type, sc$group))

# 各细胞类型中的表达（检测比例与平均表达），含周细胞和神经胶质
lognorm <- GetAssayData(sc, layer = "data")
check_genes <- c(genes_focus, "RGS5", "NOTCH3", "PLP1", "S100B", "CHGA", "CHGB", "SCG3", "SCG5")
expr_by_type <- map_dfr(intersect(check_genes, rownames(sc)), function(g) {
  x <- lognorm[g, ]
  tibble(gene = g, cell_type = sc$cell_type, v = x) %>% group_by(gene, cell_type) %>%
    summarise(n_cells = n(), pct_expr = mean(v > 0) * 100, mean_expr = mean(v), .groups = "drop")
})
write_csv(expr_by_type, file.path(tab, "rev_sc1_expression_by_type.csv"))
cat("\nSCG2 / COL10A1 在各细胞类型的检测比例（%）：\n")
print(expr_by_type %>% filter(gene %in% c("SCG2", "COL10A1", "CHGA")) %>%
        select(gene, cell_type, n_cells, pct_expr) %>% pivot_wider(names_from = gene, values_from = pct_expr) %>%
        mutate(across(where(is.numeric), ~ round(.x, 1))))

# ---------- 3. 伪 bulk：counts 求和 + edgeR QL；精确 Wilcoxon 敏感性分析 ----------
raw_counts <- GetAssayData(sc, layer = "counts")
pb_test <- function(ct, min_cells = 20) {
  cells <- which(sc$cell_type == ct)
  md <- sc@meta.data[cells, ] %>% mutate(cell = colnames(sc)[cells])
  ok <- md %>% count(sample, group) %>% filter(n >= min_cells)
  if (n_distinct(ok$group) < 2 || any(table(ok$group)[c("Control", "NP")] < 3, na.rm = TRUE)) return(NULL)
  agg <- sapply(ok$sample, function(s) Matrix::rowSums(raw_counts[, md$cell[md$sample == s], drop = FALSE]))
  grp <- droplevels(ok$group)
  y <- DGEList(agg, group = grp); keep <- filterByExpr(y, group = grp)
  y <- normLibSizes(y[keep, , keep.lib.sizes = FALSE])
  des <- model.matrix(~ 0 + grp); colnames(des) <- make.names(levels(grp))
  fit <- glmQLFit(estimateDisp(y, des), des, robust = TRUE)
  lc <- cpm(y, log = TRUE, prior.count = 1)
  cons <- c(`NP vs Control` = "NP - Control", `CRSwNP-Eth vs Control` = "CRSwNP.Eth - Control",
            `CRSsNP vs Control` = "CRSsNP - Control", `NP vs CRSwNP-Eth` = "NP - CRSwNP.Eth")
  map_dfr(names(cons), function(nm) {
    parts <- strsplit(cons[[nm]], " - ")[[1]]
    if (!all(parts %in% colnames(des)) || any(colSums(des[, parts, drop = FALSE]) < 3)) return(NULL)
    tt <- topTags(glmQLFTest(fit, contrast = makeContrasts(contrasts = cons[[nm]], levels = des)), n = Inf)$table
    g <- intersect(genes_focus, rownames(tt))
    a1 <- make.names(as.character(grp)) == parts[1]; a0 <- make.names(as.character(grp)) == parts[2]
    tibble(cell_type = ct, comparison = nm, gene = g, n1 = sum(a1), n0 = sum(a0),
           logFC = tt[g, "logFC"], P_QL = tt[g, "PValue"], FDR_genomewide = tt[g, "FDR"],
           P_exact_wilcoxon = map_dbl(g, ~ wilcox.test(lc[.x, a1], lc[.x, a0], exact = TRUE)$p.value),
           genes_tested = nrow(tt))
  })
}
pb <- map_dfr(setdiff(levels(sc$cell_type), c("Schwann/glia", "Lymphatic")), pb_test)
pb <- pb %>% group_by(comparison) %>% mutate(FDR_focus = p.adjust(P_QL, "BH")) %>% ungroup()
write_csv(pb, file.path(tab, "rev_sc1_pseudobulk_edgeR.csv"))
cat("\n伪 bulk（counts 求和 + edgeR QL）：\n")
print(pb %>% filter(gene %in% c("SCG2", "COL10A1", "PTHLH"), cell_type %in% c("Fibroblast", "Epithelial", "Pericyte/SMC")) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))), n = 60, width = 200)

# ---------- 4. 成纤维细胞亚群 ----------
fb <- subset(sc, cell_type == "Fibroblast")
fb <- cluster_harmony(fb, nfeat = 2000, dims = 15, res = 0.4)
fb$fsub <- paste0("F", fb$seurat_clusters)
pub_sets <- list(
  `IDO1+ (Xie 2024)` = c("IDO1", "CXCL1", "CXCL2", "CXCL3", "CXCL5", "CXCL6", "CXCL8"),
  `IL-4-associated C8 (Cao 2026)` = c("CCL26", "POSTN", "CISH", "CCL13"),
  `CXCL12+ (Guo 2024)` = c("CXCL12"),
  `PI16+ universal` = c("PI16", "DPP4", "CD34"),
  Myofibroblast = c("ACTA2", "TAGLN", "MYL9", "POSTN"))
fb <- AddModuleScore(fb, map(pub_sets, ~ intersect(.x, rownames(fb))), name = "pub_", seed = 2026)
pub_cols <- paste0("pub_", seq_along(pub_sets))
fl <- GetAssayData(fb, layer = "data")
fsum <- fb@meta.data %>% mutate(SCG2 = fl["SCG2", ] > 0, COL10A1 = fl["COL10A1", ] > 0) %>%
  group_by(fsub) %>%
  summarise(n_cells = n(), n_samples = n_distinct(sample), pct_from_NP = mean(group == "NP") * 100,
            pct_SCG2 = mean(SCG2) * 100, pct_COL10A1 = mean(COL10A1) * 100,
            across(all_of(pub_cols), mean), .groups = "drop") %>%
  rename_with(~ names(pub_sets)[as.integer(sub("pub_", "", .x))], all_of(pub_cols))
fmk <- FindAllMarkers(fb, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.5, verbose = FALSE) %>%
  group_by(cluster) %>% slice_max(avg_log2FC, n = 8) %>% summarise(top_markers = paste(gene, collapse = ", "))
fsum <- fsum %>% left_join(fmk %>% mutate(fsub = paste0("F", cluster)) %>% select(fsub, top_markers), by = "fsub")
write_csv(fsum, file.path(tab, "rev_sc1_fibroblast_subclusters.csv"))
cat("\n成纤维细胞亚群：\n"); print(fsum %>% mutate(across(where(is.numeric), ~ round(.x, 2))), width = 250)
# 每个样本中各亚群的比例
fprop <- fb@meta.data %>% count(sample, group, fsub) %>% group_by(sample) %>% mutate(frac = n / sum(n)) %>% ungroup()
write_csv(fprop, file.path(tab, "rev_sc1_fibroblast_subcluster_fraction.csv"))
# SCG2+/COL10A1+ 细胞与已发表状态评分的相关
cor_pub <- map_dfr(names(pub_sets), function(nm) {
  sc_col <- fb@meta.data[[pub_cols[match(nm, names(pub_sets))]]]
  tibble(signature = nm, rho_SCG2 = cor(sc_col, fl["SCG2", ], method = "spearman"),
         rho_COL10A1 = cor(sc_col, fl["COL10A1", ], method = "spearman"))
})
write_csv(cor_pub, file.path(tab, "rev_sc1_fibroblast_signature_correlation.csv"))
print(cor_pub)

# ---------- 5. SCG2 与 COL10A1 共表达 ----------
co <- fb@meta.data %>% mutate(S = fl["SCG2", ] > 0, C = fl["COL10A1", ] > 0) %>%
  group_by(sample, group) %>%
  summarise(n = n(), both = sum(S & C), S_only = sum(S & !C), C_only = sum(!S & C), none = sum(!S & !C),
            pct_SCG2 = mean(S) * 100, pct_COL10A1 = mean(C) * 100,
            pct_COL10A1_in_SCG2pos = if (sum(S)) mean(C[S]) * 100 else NA_real_,
            pct_COL10A1_in_SCG2neg = if (sum(!S)) mean(C[!S]) * 100 else NA_real_, .groups = "drop")
arr <- array(unlist(map(seq_len(nrow(co)), ~ matrix(c(co$both[.x], co$S_only[.x], co$C_only[.x], co$none[.x]), 2) + 0.5)),
             dim = c(2, 2, nrow(co)))
mh <- mantelhaen.test(arr)
cat("\n共表达（成纤维细胞内）：Mantel-Haenszel 合并 OR =", round(mh$estimate, 2),
    "95% CI", round(mh$conf.int, 2), "P =", signif(mh$p.value, 3), "\n")
print(co)
write_csv(co %>% mutate(MH_OR = mh$estimate, MH_lo = mh$conf.int[1], MH_hi = mh$conf.int[2], MH_P = mh$p.value),
          file.path(tab, "rev_sc1_coexpression.csv"))

# ---------- 6. 作图 ----------
fb_umap <- bind_cols(as_tibble(Embeddings(fb, "umap"), .name_repair = ~ c("UMAP1", "UMAP2")),
                     fsub = fb$fsub, group = fb$group, SCG2 = fl["SCG2", ], COL10A1 = fl["COL10A1", ])
lab <- fb_umap %>% group_by(fsub) %>% summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
p1 <- ggplot(fb_umap, aes(UMAP1, UMAP2, colour = fsub)) + geom_point(size = 0.4) +
  geom_text(data = lab, aes(label = fsub), colour = "black", size = 2.6, fontface = "bold") +
  labs(title = "Fibroblast subclusters") + theme_pub() + theme(legend.position = "none")
p2 <- ggplot(fb_umap, aes(UMAP1, UMAP2, colour = group)) + geom_point(size = 0.4) +
  scale_colour_manual(values = col_group, name = NULL) + labs(title = "Group") + theme_pub()
pf <- map(c("SCG2", "COL10A1"), function(g) ggplot(arrange(fb_umap, .data[[g]]), aes(UMAP1, UMAP2, colour = .data[[g]])) +
  geom_point(size = 0.4) + scale_colour_gradient(low = "#e1e0d9", high = "#e34948", name = NULL) +
  labs(title = g) + theme_pub() + theme(plot.title = element_text(face = "bold.italic")))
heat <- fsum %>% select(fsub, all_of(names(pub_sets))) %>% pivot_longer(-fsub) %>%
  group_by(name) %>% mutate(z = as.numeric(scale(value))) %>% ungroup()
p3 <- ggplot(heat, aes(fsub, name, fill = z)) + geom_tile() +
  scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#e34948", name = "Scaled\nscore") +
  labs(x = NULL, y = NULL, title = "Published fibroblast states") + theme_pub()
pb_plot <- pb %>% filter(gene %in% c("SCG2", "COL10A1"), cell_type == "Fibroblast")
save_fig((p1 | p2) / (pf[[1]] | pf[[2]]) / p3 + plot_layout(heights = c(1, 1, 0.7)), "FigR_sc1_fibroblast_subsets", 180, 200)

saveRDS(sc, "data/processed/scRNA_GSE235711_all20.rds")
saveRDS(fb, "data/processed/scRNA_GSE235711_fibroblasts.rds")
cat("\n完成\n")
