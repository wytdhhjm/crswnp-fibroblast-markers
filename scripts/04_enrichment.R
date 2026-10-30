# 04_enrichment.R
# 第 2 步（上）：稳健差异基因的功能富集
#   1) GO 生物过程 + KEGG 过表达分析（ORA），上调、下调基因分开做
#      背景基因 = 5 个数据集共同检测到的基因（不是全基因组，避免富集结果虚高）
#   2) GSEA：所有基因按荟萃 z 值排序，用 MSigDB Hallmark 50 个基因集
#
# 输入：results/tables/meta_all_genes.csv、data/genesets/（Hallmark GMT、KEGG 缓存）
# 输出：results/tables/enrich_*.csv、results/figures/Fig2_*

suppressPackageStartupMessages({
  library(tidyverse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(patchwork)
})

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
tab <- "results/tables"
fig <- "results/figures"
set.seed(2026)

col_dir <- c(Down = "#2a78d6", Up = "#e34948")
ink2 <- "#52514e"
theme_pub <- function(base = 9) {
  theme_classic(base_size = base) +
    theme(axis.line = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.ticks = element_line(colour = "#c3c2b7", linewidth = 0.4),
          axis.text = element_text(colour = ink2),
          plot.title = element_text(face = "bold", size = base + 1),
          legend.key.size = unit(3, "mm"))
}
save_fig <- function(p, name, w, h) {
  ggsave(file.path(fig, paste0(name, ".pdf")), p, width = w, height = h, units = "mm")
  ggsave(file.path(fig, paste0(name, ".png")), p, width = w, height = h, units = "mm", dpi = 300, bg = "white")
}

meta <- read_csv(file.path(tab, "meta_all_genes.csv"), show_col_types = FALSE)
universe <- meta$gene
up   <- meta$gene[meta$direction == "Up"]
down <- meta$gene[meta$direction == "Down"]

sym2eg <- bitr(universe, "SYMBOL", "ENTREZID", org.Hs.eg.db) %>% distinct(SYMBOL, .keep_all = TRUE)
eg <- function(s) sym2eg$ENTREZID[match(s, sym2eg$SYMBOL)] %>% na.omit() %>% as.character()

# KEGG 通路注释：使用本地缓存（rest.kegg.jp 下载，日期见 kegg_download_date.txt），保证可复现
kegg_t2g <- read_tsv("data/genesets/kegg_hsa_link.txt", col_names = c("path", "gene"), show_col_types = FALSE) %>%
  transmute(path = sub("^path:", "", path), gene = sub("^hsa:", "", gene))
kegg_t2n <- read_tsv("data/genesets/kegg_hsa_list.txt", col_names = c("path", "name"), show_col_types = FALSE) %>%
  mutate(name = sub(" - Homo sapiens \\(human\\)$", "", name))

# ---------- 1. ORA ----------
run_ora <- function(genes, label) {
  go <- enrichGO(genes, OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont = "BP",
                 universe = universe, pvalueCutoff = 0.05, qvalueCutoff = 0.05)
  go <- simplify(go, cutoff = 0.7)                    # 合并语义高度重复的 GO 条目
  kg <- enricher(eg(genes), TERM2GENE = kegg_t2g, TERM2NAME = kegg_t2n,
                 universe = eg(universe), pvalueCutoff = 0.05)
  kg <- setReadable(kg, org.Hs.eg.db, "ENTREZID")
  bind_rows(as_tibble(go@result) %>% mutate(db = "GO BP"),
            as_tibble(kg@result) %>% mutate(db = "KEGG")) %>%
    filter(p.adjust < 0.05) %>%
    mutate(set = label)
}
ora <- bind_rows(run_ora(up, "Up"), run_ora(down, "Down"))
write_csv(ora, file.path(tab, "enrich_ORA_GO_KEGG.csv"))
cat("ORA 显著条目数：\n"); print(count(ora, set, db))

# 图：每组各取 GO、KEGG 前 8 条
top_ora <- ora %>%
  group_by(set, db) %>% slice_min(p.adjust, n = 8, with_ties = FALSE) %>% ungroup() %>%
  mutate(ratio = map_dbl(GeneRatio, ~ eval(parse(text = .x))),
         Description = str_trunc(Description, 55),
         label = fct_reorder(paste(Description, set, db, sep = "|"), -log10(p.adjust)))
p_ora <- ggplot(top_ora, aes(-log10(p.adjust), label)) +
  geom_col(aes(fill = set), width = 0.7) +
  geom_text(aes(label = Count), hjust = -0.2, size = 2.3, colour = ink2) +
  scale_fill_manual(values = col_dir, guide = "none") +
  scale_y_discrete(labels = function(x) sub("\\|.*", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
  facet_grid(set + db ~ ., scales = "free_y", space = "free_y",
             labeller = labeller(.multi_line = FALSE)) +
  labs(x = expression(-log[10]~"(adjusted P)"), y = NULL) +
  theme_pub(8) +
  theme(strip.background = element_blank(), strip.text.y = element_text(angle = 0, face = "bold", hjust = 0))
save_fig(p_ora, "Fig2A_ORA_GO_KEGG", 170, 190)

# ---------- 2. GSEA（Hallmark） ----------
gmt <- read.gmt("data/genesets/h.all.v2024.1.Hs.symbols.gmt")
ranks <- meta %>% arrange(desc(z)) %>% { setNames(.$z, .$gene) }
gs <- GSEA(ranks, TERM2GENE = gmt, pvalueCutoff = 1, minGSSize = 15, maxGSSize = 500,
           eps = 0, seed = TRUE)
gsea_res <- as_tibble(gs@result) %>%
  mutate(Description = str_remove(Description, "HALLMARK_") %>% str_replace_all("_", " ")) %>%
  arrange(desc(NES))
write_csv(gsea_res, file.path(tab, "enrich_GSEA_hallmark.csv"))
cat("Hallmark GSEA 显著（FDR < 0.05）：", sum(gsea_res$p.adjust < 0.05), "个\n")

gplot <- gsea_res %>% filter(p.adjust < 0.05) %>%
  mutate(dir = ifelse(NES > 0, "Up", "Down"), Description = fct_reorder(Description, NES))
p_gsea <- ggplot(gplot, aes(NES, Description)) +
  geom_vline(xintercept = 0, colour = "#c3c2b7", linewidth = 0.3) +
  geom_segment(aes(x = 0, xend = NES, yend = Description), colour = "#c3c2b7", linewidth = 0.4) +
  geom_point(aes(colour = dir, size = -log10(p.adjust))) +
  scale_colour_manual(values = col_dir, name = NULL,
                      labels = c(Down = "Enriched in control", Up = "Enriched in NP")) +
  scale_size_continuous(range = c(1.5, 4), name = expression(-log[10]~"FDR")) +
  labs(x = "Normalized enrichment score (NES)", y = NULL) +
  theme_pub(8)
save_fig(p_gsea, "Fig2B_GSEA_hallmark", 150, max(60, 5 * nrow(gplot) + 25))

cat("\n完成\n")
