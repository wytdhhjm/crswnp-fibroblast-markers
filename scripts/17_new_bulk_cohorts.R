# 17_new_bulk_cohorts.R
# 系统检索（16_geo_search.R）新发现的两个 bulk 组织数据集，作为独立验证与内型比较
#   GSE338426（RNA-seq read counts）：对照下鼻甲 8、嗜酸性 CRSwNP 息肉 9、非嗜酸性 CRSwNP 息肉 5
#   GSE72713 （RNA-seq FPKM）       ：对照蝶窦黏膜 3、嗜酸性息肉 3、非嗜酸性息肉 3（样本极少，仅作补充）
# 两个数据集均未参与任何基因筛选，结果是真正的事后外部验证
#   1) 共识基因与经典 2 型标志物在 息肉 vs 对照、嗜酸性 vs 非嗜酸性 中的差异（limma-voom / limma-trend）
#   2) 共识基因评分与随机 4 基因评分的 AUC（与 14 相同的基线比较）
#
# 输出：results/tables/rev_E_new_cohorts_*.csv

suppressPackageStartupMessages({ library(tidyverse); library(edgeR); library(limma); library(pROC) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
raw <- "data/raw"; tab <- "results/tables"
set.seed(2026)
ml <- read_csv(file.path(tab, "ML_feature_selection.csv"), show_col_types = FALSE)
hub_genes <- ml$gene[ml$hub]
genes_c <- c(hub_genes, "POSTN", "CLC", "CCL26", "CST1", "ALOX15")
sym_tab <- local({ y <- read.delim(file.path(raw, "GSE179265_counts.txt.gz")); setNames(y$Gene_name, y$geneID) })

# ---------- GSE338426 ----------
x <- read_csv(file.path(raw, "GSE338426_Expression.csv.gz"), show_col_types = FALSE)
cnt <- as.matrix(x[, grep(":read count$", names(x))]); colnames(cnt) <- sub(":read count$", "", colnames(cnt))
cnt <- cnt[, colnames(cnt) != "read_count"]                 # 最后一列是重复的 Control1，丢弃
sym <- ifelse(!is.na(sym_tab[x$Gene_ID]), sym_tab[x$Gene_ID], x$Name)
sym <- limma::alias2SymbolTable(sym); keep <- !is.na(sym)
cnt <- cnt[keep, ]; sym <- sym[keep]; o <- order(rowSums(cnt), decreasing = TRUE)
cnt <- cnt[o, ]; sym <- sym[o]; cnt <- cnt[!duplicated(sym), ]; rownames(cnt) <- sym[!duplicated(sym)]
grp <- factor(case_when(grepl("^Control", colnames(cnt)) ~ "Control", grepl("^EosNP", colnames(cnt)) ~ "eNP",
                        grepl("^Non_EP", colnames(cnt)) ~ "neNP"), c("Control", "eNP", "neNP"))
print(table(grp))
y <- DGEList(cnt, group = grp); y <- normLibSizes(y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE])
des <- model.matrix(~ 0 + grp); colnames(des) <- levels(grp)
v <- voom(y, des); fit <- lmFit(v, des)
cons <- c(`All NP vs control` = "(eNP + neNP)/2 - Control", `Eos NP vs control` = "eNP - Control",
          `Non-eos NP vs control` = "neNP - Control", `Eos vs non-eos NP` = "eNP - neNP")
resE <- map_dfr(names(cons), function(nm) {
  tt <- topTable(eBayes(contrasts.fit(fit, makeContrasts(contrasts = cons[[nm]], levels = des))), number = Inf, sort.by = "none")
  g <- intersect(genes_c, rownames(tt))
  tibble(dataset = "GSE338426", comparison = nm, gene = g, logFC = tt[g, "logFC"], P = tt[g, "P.Value"], FDR = tt[g, "adj.P.Val"])
})

# 评分 AUC：息肉（全部）vs 对照；随机 4 基因基线
z <- t(scale(t(v$E))); np <- factor(ifelse(grp == "Control", "Control", "NP"), c("Control", "NP"))
auc_of <- function(g) { g <- intersect(g, rownames(z)); as.numeric(auc(roc(np, colMeans(z[g, , drop = FALSE]), levels = c("Control", "NP"), direction = "<", quiet = TRUE))) }
up_pool <- intersect(ml$gene[ml$direction == "Up"], rownames(z))
rand <- replicate(1000, auc_of(sample(up_pool, 4)))
aucE <- tibble(set = c("Consensus 4 genes", "POSTN", "Classic 5-gene score", "Random 4 genes (median)"),
               AUC = c(auc_of(hub_genes), auc_of("POSTN"), auc_of(c("POSTN", "CLC", "CCL26", "CST1", "ALOX15")), median(rand)),
               pct_random_ge = c(mean(rand >= auc_of(hub_genes)), NA, NA, NA))
cat("\nGSE338426 AUC（息肉 vs 对照）：\n"); print(aucE)

# ---------- GSE72713（FPKM，limma-trend） ----------
f <- read.delim(file.path(raw, "GSE72713_RNA-seq_processed_data.txt.gz"), check.names = FALSE)
m <- as.matrix(f[, -(1:2)]); g72 <- limma::alias2SymbolTable(f$`Gene name`)
ok <- !is.na(g72); m <- m[ok, ]; g72 <- g72[ok]
m <- rowsum(m, g72)                                       # 同一基因的多个转录本 FPKM 相加
e <- log2(m + 1); e <- e[rowSums(m >= 1) >= 3, ]
g2 <- factor(case_when(grepl("^Ctrl", colnames(e)) ~ "Control", grepl("^ECRS", colnames(e)) ~ "eNP", TRUE ~ "neNP"),
             c("Control", "eNP", "neNP"))
d2 <- model.matrix(~ 0 + g2); colnames(d2) <- levels(g2); f2 <- lmFit(e, d2)
res72 <- map_dfr(names(cons), function(nm) {
  tt <- topTable(eBayes(contrasts.fit(f2, makeContrasts(contrasts = cons[[nm]], levels = d2)), trend = TRUE), number = Inf, sort.by = "none")
  g <- intersect(genes_c, rownames(tt))
  tibble(dataset = "GSE72713", comparison = nm, gene = g, logFC = tt[g, "logFC"], P = tt[g, "P.Value"], FDR = tt[g, "adj.P.Val"])
})

res <- bind_rows(resE, res72)
write_csv(res, file.path(tab, "rev_E_new_cohorts_DE.csv")); write_csv(aucE, file.path(tab, "rev_E_new_cohorts_AUC.csv"))
cat("\n共识基因与经典标志物（logFC / P）：\n")
print(res %>% mutate(v = sprintf("%.2f (%.2g)", logFC, P)) %>% select(dataset, comparison, gene, v) %>%
        pivot_wider(names_from = gene, values_from = v), width = 250)
