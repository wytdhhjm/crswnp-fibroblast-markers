# 19_revision_figures.R
# 修订稿新增/重做的图（全部从 results/tables 中已保存的分析结果和已保存的单细胞对象生成）
#   FigR_auc_withheld   ：每个共识基因在 3 个预留验证集中的 AUC（去掉训练集，避免循环论证）+ 评分 ROC
#   FigR_specificity    ：bulk 中疾病/部位特异性（CRSsNP、非息肉炎症黏膜、过敏性鼻炎、同一患者下鼻甲、嗜酸/非嗜酸）
#   FigR_sc1_*          ：GSE235711 全部 20 个样本：UMAP、各细胞类型表达、成纤维细胞伪 bulk（每样本）
#   FigR_sc2_*          ：GSE276503 推定供者层面的成纤维细胞伪 bulk 与阳性率
# 同时输出各细胞类型对 SCG2/COL10A1 UMI 的贡献比例（正文引用）

suppressPackageStartupMessages({ library(tidyverse); library(Seurat); library(edgeR); library(pROC); library(patchwork); library(ggrepel) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
tab <- "results/tables"; fig <- "results/figures"
ink2 <- "#52514e"
theme_pub <- function(base = 8) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.ticks = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.text = element_text(colour = ink2), plot.title = element_text(face = "bold", size = base + 1),
          legend.key.size = unit(3, "mm"), strip.background = element_blank(), strip.text = element_text(face = "bold"))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}
genes4 <- c("PTHLH", "SCG2", "COL10A1", "PDE7A")
val_ids <- c("GSE36830", "GSE23552", "GSE194282")

# ---------- 1. 预留验证集 AUC（不含训练集） ----------
ga <- read_csv(file.path(tab, "ML_hub_gene_AUC.csv"), show_col_types = FALSE) %>% filter(dataset %in% val_ids) %>%
  mutate(dataset = factor(dataset, rev(val_ids)), gene = factor(gene, genes4))
p_auc <- ggplot(ga, aes(AUC, dataset)) +
  geom_vline(xintercept = 0.5, colour = "#c3c2b7", linewidth = 0.3) +
  geom_errorbar(aes(xmin = lo, xmax = pmin(hi, 1)), width = 0, linewidth = 0.4, colour = ink2, orientation = "y") +
  geom_point(size = 1.6, shape = 17, colour = "#0b0b0b") +
  scale_x_continuous(limits = c(0.2, 1), breaks = c(0.5, 0.75, 1)) +
  facet_wrap(~ gene, nrow = 1) + labs(x = "AUC (95% CI, DeLong)", y = NULL) + theme_pub() +
  theme(strip.text = element_text(face = "bold.italic"))
save_fig(p_auc, "FigR_auc_withheld", 180, 45)

# ---------- 2. 疾病/部位特异性 ----------
c1 <- read_csv(file.path(tab, "rev_C_specificity.csv"), show_col_types = FALSE) %>%
  transmute(dataset, comparison, gene, logFC, P)
d1 <- read_csv(file.path(tab, "rev_D_site_paired_GSE136825.csv"), show_col_types = FALSE) %>%
  transmute(dataset = "GSE136825", comparison = sub(" \\(.*", "", comparison), gene, logFC, P)
e1 <- read_csv(file.path(tab, "rev_E_new_cohorts_DE.csv"), show_col_types = FALSE) %>%
  filter(dataset == "GSE338426") %>% transmute(dataset, comparison, gene, logFC, P)
sp <- bind_rows(c1, d1, e1) %>% filter(gene %in% c(genes4, "POSTN", "CLC")) %>%
  mutate(row = paste0(comparison, "  [", dataset, "]"),
         sig = case_when(P < 0.001 ~ "***", P < 0.01 ~ "**", P < 0.05 ~ "*", TRUE ~ ""))
row_order <- c("Polyp vs control UT  [GSE36830]", "Polyp vs normal  [GSE23552]", "All NP vs control  [GSE338426]",
               "Polyp vs own IT  [GSE136825]", "Polyp vs CRSwNP UT  [GSE36830]", "Polyp vs CRSwNP inflamed mucosa  [GSE23552]",
               "CRSwNP inflamed mucosa vs normal  [GSE23552]", "CRSwNP UT vs control UT  [GSE36830]",
               "CRSwNP IT vs control IT  [GSE136825]", "CRSsNP UT vs control UT  [GSE36830]",
               "Allergic rhinitis vs normal  [GSE23552]", "Eos vs non-eos NP  [GSE338426]")
sp <- sp %>% filter(row %in% row_order) %>% mutate(row = factor(row, rev(row_order)), gene = factor(gene, c(genes4, "POSTN", "CLC")))
lim <- max(abs(sp$logFC)); lim <- ceiling(min(lim, 4))
p_sp <- ggplot(sp, aes(gene, row, fill = pmax(pmin(logFC, lim), -lim))) + geom_tile(colour = "white") +
  geom_text(aes(label = sig), size = 2.6, vjust = 0.75) +
  scale_fill_gradient2(low = "#2a78d6", mid = "white", high = "#e34948", limits = c(-lim, lim), name = expression(log[2]*FC)) +
  geom_vline(xintercept = 4.5, colour = ink2, linewidth = 0.3) +
  geom_hline(yintercept = c(1.5, 3.5, 6.5, 9.5) - 0, colour = "#c3c2b7", linewidth = 0.3) +
  labs(x = NULL, y = NULL) + theme_pub() + theme(axis.text.x = element_text(face = "italic"), axis.line = element_blank())
save_fig(p_sp, "FigR_specificity", 180, 95)

# ---------- 3. GSE235711（20 个样本） ----------
sc <- readRDS("data/processed/scRNA_GSE235711_all20.rds")
ct_cols <- c(Epithelial = "#2a78d6", Glandular = "#86b6ef", Fibroblast = "#eb6834", `Pericyte/SMC` = "#eda100",
             Endothelial = "#1baf7a", Lymphatic = "#008300", `T/NK` = "#4a3aa7", `B/Plasma` = "#e87ba4",
             Myeloid = "#e34948", Mast = "#898781", Granulocyte = "#52514e", `Schwann/glia` = "#c3c2b7")
ud <- bind_cols(as_tibble(Embeddings(sc, "umap"), .name_repair = ~ c("UMAP1", "UMAP2")), cell_type = sc$cell_type)
lp <- ud %>% group_by(cell_type) %>% summarise(UMAP1 = median(UMAP1), UMAP2 = median(UMAP2))
p_umap <- ggplot(ud, aes(UMAP1, UMAP2, colour = cell_type)) + geom_point(size = 0.03, alpha = 0.4) +
  geom_text_repel(data = lp, aes(label = cell_type), colour = "black", size = 2.2, fontface = "bold",
                  bg.color = "white", bg.r = 0.15, min.segment.length = 0, segment.size = 0.2, seed = 2026) +
  scale_colour_manual(values = ct_cols, guide = "none") +
  labs(title = sprintf("GSE235711 (%s cells, %d samples)", format(ncol(sc), big.mark = ","), n_distinct(sc$sample))) +
  theme_pub() + theme(axis.text = element_blank(), axis.ticks = element_blank())
p_dot <- DotPlot(sc, features = genes4, group.by = "cell_type", cols = c("#e1e0d9", "#e34948")) +
  labs(x = NULL, y = NULL) + theme_pub() + theme(axis.text.x = element_text(face = "italic"))
umi_share <- function(obj, ctcol) {
  cn <- GetAssayData(obj, layer = "counts")[intersect(genes4, rownames(obj)), ]
  as_tibble(t(as.matrix(sapply(split(seq_len(ncol(obj)), obj@meta.data[[ctcol]]), function(i) Matrix::rowSums(cn[, i, drop = FALSE])))),
            rownames = "cell_type") %>% mutate(across(-cell_type, ~ .x / sum(.x) * 100))
}
sh1 <- umi_share(sc, "cell_type"); write_csv(sh1, file.path(tab, "rev_sc1_umi_share_by_type.csv"))
cat("GSE235711 各细胞类型 UMI 贡献（%）：\n"); print(sh1 %>% mutate(across(where(is.numeric), ~ round(.x, 1))))
# 成纤维细胞每样本 counts 求和 → logCPM
fib <- sc$cell_type == "Fibroblast"; cn <- GetAssayData(sc, layer = "counts")[, fib]; md <- sc@meta.data[fib, ]
ok <- names(which(table(md$sample) >= 20))
pbm <- sapply(ok, function(s) Matrix::rowSums(cn[, md$sample == s, drop = FALSE]))
lc <- cpm(normLibSizes(DGEList(pbm)), log = TRUE, prior.count = 1)
grp_lab <- c(Control = "Control\nethmoid", CRSsNP = "CRSsNP\nethmoid", `CRSwNP-Eth` = "CRSwNP\nethmoid", NP = "Nasal\npolyp")
pd1 <- expand_grid(gene = c("SCG2", "COL10A1"), sample = ok) %>%
  mutate(logCPM = map2_dbl(gene, sample, ~ lc[.x, .y]), group = as.character(md$group[match(sample, md$sample)])) %>%
  mutate(group = factor(grp_lab[group], grp_lab), gene = factor(gene, c("SCG2", "COL10A1")))
pt1 <- read_csv(file.path(tab, "rev_sc1_pseudobulk_edgeR.csv"), show_col_types = FALSE) %>%
  filter(cell_type == "Fibroblast", gene %in% c("SCG2", "COL10A1"))
p_pb1 <- ggplot(pd1, aes(group, logCPM, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "#898781", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.4) +
  scale_colour_manual(values = c("#2a78d6", "#1baf7a", "#eda100", "#eb6834"), guide = "none") +
  facet_wrap(~ gene, scales = "free_y") + labs(x = NULL, y = "Fibroblast pseudobulk log2 CPM") + theme_pub() +
  theme(strip.text = element_text(face = "bold.italic"))
save_fig((p_umap | p_dot) / p_pb1 + plot_layout(heights = c(1.2, 1)), "FigR_sc1_overview", 180, 160)

# ---------- 4. GSE276503（推定供者层面） ----------
st <- readRDS("data/processed/scRNA_GSE276503_stromal.rds")
dm <- read_csv(file.path(tab, "rev_sc2_donor_map.csv"), show_col_types = FALSE) %>% filter(threshold == 0.94)
fib2 <- st$sub_type == "Fibroblast"; cn2 <- GetAssayData(st, layer = "counts")[, fib2]
unit <- dm$unit[match(st$sample[fib2], dm$sample)]
pbm2 <- sapply(sort(unique(unit)), function(u) Matrix::rowSums(cn2[, unit == u, drop = FALSE]))
lc2 <- cpm(normLibSizes(DGEList(pbm2)), log = TRUE, prior.count = 1)
grp2_lab <- c(`Healthy IT` = "Healthy\nIT", `CRSwNP IT` = "CRSwNP\nIT", `CRSwNP MT` = "CRSwNP\nMT", NP = "Nasal\npolyp")
pd2 <- expand_grid(gene = c("SCG2", "COL10A1"), unit = colnames(lc2)) %>%
  mutate(logCPM = map2_dbl(gene, unit, ~ lc2[.x, .y]), group = sub(".*\\|", "", unit)) %>%
  mutate(group = factor(grp2_lab[group], grp2_lab), gene = factor(gene, c("SCG2", "COL10A1")))
p_pb2 <- ggplot(pd2, aes(group, logCPM, colour = group)) +
  geom_boxplot(outlier.shape = NA, width = 0.55, colour = "#898781", linewidth = 0.3) +
  geom_jitter(width = 0.12, size = 1.4) +
  scale_colour_manual(values = c("#2a78d6", "#1baf7a", "#eda100", "#eb6834"), guide = "none") +
  facet_wrap(~ gene, scales = "free_y") + labs(x = NULL, y = "Fibroblast pseudobulk log2 CPM\n(one point per inferred donor)") +
  theme_pub() + theme(strip.text = element_text(face = "bold.italic"))
sh2 <- read_csv(file.path(tab, "rev_sc2_pseudobulk_alltypes.csv"), show_col_types = FALSE) %>%
  group_by(cell_type) %>% summarise(across(any_of(genes4), sum)) %>% mutate(across(-cell_type, ~ .x / sum(.x) * 100)) %>%
  pivot_longer(-cell_type, names_to = "gene", values_to = "pct") %>%
  mutate(cell_type = ifelse(cell_type == "Stromal", "Stromal (fibroblast,\npericyte)", cell_type),
         gene = factor(gene, genes4))
p_sh2 <- ggplot(sh2, aes(pct, gene, fill = cell_type)) + geom_col(width = 0.7) +
  scale_fill_manual(values = c(Epithelial = "#2a78d6", Glandular = "#86b6ef", `Stromal (fibroblast,\npericyte)` = "#eb6834",
                               Endothelial = "#1baf7a", `T/NK` = "#4a3aa7", `B/Plasma` = "#e87ba4", Myeloid = "#e34948",
                               Mast = "#898781", Granulocyte = "#52514e", Glia = "#c3c2b7"), name = NULL) +
  scale_y_discrete(limits = rev) + labs(x = "Share of UMIs (%)", y = NULL,
                                        title = "GSE276503: cellular source (all 29 samples)") +
  theme_pub() + theme(axis.text.y = element_text(face = "italic"))
save_fig(p_sh2 / p_pb2 + plot_layout(heights = c(0.8, 1)), "FigR_sc2_replication", 180, 130)
cat("\n完成\n")
