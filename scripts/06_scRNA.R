# 06_scRNA.R
# 第 3 步：单细胞定位核心基因（GSE235711，Immunity 2025）
#
# 样本：健康对照筛窦黏膜 5、CRSwNP 鼻息肉 6（主分析）；
#       CRSsNP 筛窦 5、CRSwNP 筛窦 6（补充，下载完成则一并纳入）
# 流程：
#   1) 读取 10x 过滤矩阵（仅取 Gene Expression，忽略 CITE-seq 抗体标签）
#   2) 质控：200 < 基因数 < 7500，线粒体比例 < 20%（预设标准，按常规阈值）
#   3) 标准化 → 高变基因 → PCA → Harmony 按样本去批次 → 聚类 → UMAP
#   4) 用经典标志基因注释大类细胞（上皮、成纤维、内皮、平滑肌/周细胞、T/NK、B/浆细胞、
#      髓系、肥大细胞、嗜酸/中性粒）
#   5) 核心基因在各细胞类型中的表达；按样本做"伪 bulk"比较息肉 vs 对照
#      （以样本而非细胞为单位做统计，避免把成千上万个细胞当独立样本导致 P 值虚低）
#
# 输出：data/processed/scRNA_GSE235711.rds、results/tables/sc_*.csv、results/figures/Fig4_*

suppressPackageStartupMessages({
  library(tidyverse)
  library(Seurat)
  library(harmony)
  library(patchwork)
})

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
raw <- "data/raw/GSE235711"
tab <- "results/tables"
fig <- "results/figures"
set.seed(2026)
options(future.globals.maxSize = 8 * 1024^3)

hub_genes <- read_csv(file.path(tab, "ML_feature_selection.csv"), show_col_types = FALSE) %>%
  filter(hub) %>% pull(gene)
cat("核心基因：", hub_genes, "\n")

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
col_group <- c(Control = "#2a78d6", CRSsNP = "#1baf7a", `CRSwNP-Eth` = "#eda100", NP = "#eb6834")

# ---------- 1. 读取 ----------
files <- list.files(raw, pattern = "_filtered\\.tar\\.gz$", full.names = TRUE)
cat("已下载样本：", length(files), "\n")
tmp <- file.path(tempdir(), "sc10x")
objs <- map(files, function(f) {
  id <- sub("_filtered\\.tar\\.gz$", "", basename(f))
  d <- file.path(tmp, id); dir.create(d, recursive = TRUE, showWarnings = FALSE)
  untar(f, exdir = d)
  m <- suppressMessages(Read10X(d))
  if (is.list(m)) m <- m[["Gene Expression"]]
  s <- CreateSeuratObject(m, project = id, min.cells = 3, min.features = 200)
  s$sample <- sub("^GSM[0-9]+_", "", id)
  s
})
sc <- merge(objs[[1]], objs[-1], add.cell.ids = map_chr(objs, ~ .x$sample[1]))
sc$group <- case_when(grepl("^Control", sc$sample) ~ "Control",
                      grepl("^CRSsNP", sc$sample) ~ "CRSsNP",
                      grepl("^CRSwNP_Eth", sc$sample) ~ "CRSwNP-Eth",
                      grepl("^CRSwNP_NP", sc$sample) ~ "NP")
sc$group <- factor(sc$group, names(col_group))

# ---------- 2. 质控 ----------
sc[["mt"]] <- PercentageFeatureSet(sc, "^MT-")
n0 <- ncol(sc)
sc <- subset(sc, nFeature_RNA > 200 & nFeature_RNA < 7500 & mt < 20)
qc_tab <- sc@meta.data %>% count(sample, group, name = "cells_after_QC")
write_csv(qc_tab, file.path(tab, "sc_QC_cells_per_sample.csv"))
cat("质控：", n0, "->", ncol(sc), "个细胞\n")

# ---------- 3. 降维、去批次、聚类 ----------
sc <- JoinLayers(sc)
sc <- NormalizeData(sc, verbose = FALSE) %>%
  FindVariableFeatures(nfeatures = 3000, verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(npcs = 30, verbose = FALSE)
sc <- RunHarmony(sc, "sample", verbose = FALSE)
sc <- RunUMAP(sc, reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindNeighbors(reduction = "harmony", dims = 1:30, verbose = FALSE) %>%
  FindClusters(resolution = 0.5, verbose = FALSE)

# ---------- 4. 细胞类型注释（标志基因评分，按聚类取最高者） ----------
markers <- list(
  Epithelial       = c("EPCAM", "KRT5", "KRT8", "KRT19", "SCGB1A1", "FOXJ1"),
  Fibroblast       = c("COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"),
  Endothelial      = c("PECAM1", "VWF", "CDH5", "PLVAP"),
  `Smooth muscle/Pericyte` = c("ACTA2", "RGS5", "MYH11", "TAGLN"),
  `Glandular`      = c("LTF", "BPIFB2", "ZG16B", "PIP"),
  `T/NK`           = c("CD3E", "CD3D", "TRAC", "NKG7", "GNLY"),
  `B/Plasma`       = c("MS4A1", "CD79A", "JCHAIN", "MZB1"),
  Myeloid          = c("CD68", "CD14", "LYZ", "C1QA", "CD1C"),
  Mast             = c("TPSAB1", "CPA3", "KIT", "HDC"),
  Granulocyte      = c("CLC", "CSF3R", "FCGR3B", "SIGLEC8")
)
markers <- map(markers, ~ intersect(.x, rownames(sc)))
sc <- AddModuleScore(sc, markers, name = "ct_", seed = 2026)
score_cols <- paste0("ct_", seq_along(markers))
cl_scores <- sc@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean), .groups = "drop")
cl_type <- tibble(seurat_clusters = cl_scores$seurat_clusters,
                  cell_type = names(markers)[max.col(as.matrix(cl_scores[, score_cols]))])
sc$cell_type <- cl_type$cell_type[match(sc$seurat_clusters, cl_type$seurat_clusters)]
sc$cell_type <- factor(sc$cell_type, names(markers))
write_csv(cl_type, file.path(tab, "sc_cluster_annotation.csv"))
print(table(sc$cell_type, sc$group))

saveRDS(sc, "data/processed/scRNA_GSE235711.rds")

# ---------- 5. 作图 ----------
ct_cols <- setNames(c("#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4",
                      "#008300", "#4a3aa7", "#e34948", "#86b6ef", "#898781"),
                    names(markers))

umap_df <- bind_cols(as_tibble(Embeddings(sc, "umap")), sc@meta.data %>% select(cell_type, group))
names(umap_df)[1:2] <- c("UMAP1", "UMAP2")
lab_pos <- umap_df %>% group_by(cell_type) %>% summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
p_umap <- ggplot(umap_df, aes(UMAP1, UMAP2)) +
  geom_point(aes(colour = cell_type), size = 0.05, alpha = 0.5) +
  geom_text(data = lab_pos, aes(label = cell_type), size = 2.4, fontface = "bold", colour = "#0b0b0b") +
  scale_colour_manual(values = ct_cols, guide = "none") +
  labs(title = paste0("GSE235711 (", ncol(sc), " cells, ", n_distinct(sc$sample), " samples)")) +
  theme_pub(8) + theme(axis.text = element_blank(), axis.ticks = element_blank())

# 标志基因点图（证明注释可靠）
mk_show <- map_chr(markers, 1) %>% c(map_chr(markers, 2)) %>% unique()
p_mk <- DotPlot(sc, features = mk_show, group.by = "cell_type", cols = c("#e1e0d9", "#2a78d6")) +
  RotatedAxis() + labs(x = NULL, y = NULL) + theme_pub(7) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))

# 核心基因点图：各细胞类型 × 分组
sc$ct_group <- paste(sc$cell_type, sc$group, sep = " | ")
p_hub <- DotPlot(sc, features = hub_genes, group.by = "cell_type", split.by = NULL,
                 cols = c("#e1e0d9", "#e34948")) +
  labs(x = NULL, y = NULL, title = "Hub genes across cell types") + theme_pub(8) +
  theme(axis.text.x = element_text(face = "italic"))

save_fig((p_umap | p_hub) / p_mk + plot_layout(heights = c(1.3, 1)), "Fig4A_scRNA_overview", 190, 170)

# 核心基因 UMAP 表达图
p_feat <- FeaturePlot(sc, hub_genes, ncol = length(hub_genes), order = TRUE, pt.size = 0.05,
                      cols = c("#e1e0d9", "#e34948"), combine = FALSE) %>%
  map(~ .x + theme_pub(8) + theme(axis.text = element_blank(), axis.ticks = element_blank(),
                                  plot.title = element_text(face = "bold.italic"),
                                  legend.position = "none") + labs(x = NULL, y = NULL)) %>%
  wrap_plots(nrow = 1)
save_fig(p_feat, "Fig4B_hub_featureplot", 190, 55)

# ---------- 6. 伪 bulk：每个样本 × 细胞类型的平均表达，比较 NP vs Control ----------
expr <- GetAssayData(sc, layer = "data")[hub_genes, , drop = FALSE]
pb <- map_dfr(hub_genes, function(g) {
  tibble(gene = g, value = expr[g, ], pos = expr[g, ] > 0,
         sample = sc$sample, group = sc$group, cell_type = sc$cell_type)
}) %>%
  group_by(gene, sample, group, cell_type) %>%
  summarise(n_cells = n(), mean_expr = mean(value), pct_expr = mean(pos) * 100, .groups = "drop") %>%
  filter(n_cells >= 20)                                     # 某样本中该细胞类型不足 20 个细胞时不计
write_csv(pb, file.path(tab, "sc_pseudobulk_hub_genes.csv"))

pb_test <- pb %>% filter(group %in% c("Control", "NP")) %>%
  group_by(gene, cell_type) %>%
  filter(n_distinct(group) == 2, sum(group == "Control") >= 3, sum(group == "NP") >= 3) %>%
  summarise(mean_ctrl = mean(mean_expr[group == "Control"]), mean_np = mean(mean_expr[group == "NP"]),
            P = wilcox.test(mean_expr ~ droplevels(group), exact = FALSE)$p.value, .groups = "drop") %>%
  mutate(FDR = p.adjust(P, "BH")) %>% arrange(gene, P)
write_csv(pb_test, file.path(tab, "sc_pseudobulk_test.csv"))
cat("\n伪 bulk 比较（NP vs Control，以样本为单位）：\n"); print(pb_test, n = 50)

p_pb <- pb %>% filter(group %in% c("Control", "NP")) %>%
  ggplot(aes(cell_type, mean_expr, colour = group)) +
  geom_boxplot(outlier.shape = NA, linewidth = 0.3, width = 0.6, position = position_dodge(0.75)) +
  geom_point(size = 0.8, position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.75)) +
  scale_colour_manual(values = col_group[c("Control", "NP")], name = NULL) +
  facet_wrap(~ gene, ncol = 1, scales = "free_y") +
  labs(x = NULL, y = "Mean log-normalized expression per sample") +
  theme_pub(8) + theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
save_fig(p_pb, "Fig4C_hub_pseudobulk", 170, 50 * length(hub_genes) + 20)

cat("\n完成\n")
