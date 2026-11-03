# 18b_replicate_check.R
# GSE276503 中多组样本的成纤维细胞指标几乎相同（如 NP31/NP32、NP41/NP42、HC_IT01/HC_IT02），
# 怀疑为同一患者的技术重复（同一份组织分成多个文库）。用全转录组伪 bulk 相关性判断：
# 同一患者的技术重复相关性应明显高于不同患者之间。
# 输出：results/tables/rev_sc2_sample_correlation.csv、rev_sc2_replicate_groups.csv
suppressPackageStartupMessages({ library(tidyverse); library(edgeR) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
tab <- "results/tables"
st <- readRDS("data/processed/scRNA_GSE276503_stromal.rds")
cnt <- SeuratObject::GetAssayData(st, layer = "counts")
smp <- sort(unique(st$sample))
pb <- sapply(smp, function(s) Matrix::rowSums(cnt[, st$sample == s & st$sub_type == "Fibroblast", drop = FALSE]))
y <- DGEList(pb); y <- normLibSizes(y[filterByExpr(y, min.count = 20), , keep.lib.sizes = FALSE])
lc <- cpm(y, log = TRUE, prior.count = 2)
lc <- lc[order(apply(lc, 1, var), decreasing = TRUE)[1:2000], ]
r <- cor(lc, method = "spearman")
write_csv(as_tibble(r, rownames = "sample"), file.path(tab, "rev_sc2_sample_correlation.csv"))
# 每个样本最相似的另一个样本
best <- map_dfr(smp, function(s) { v <- r[s, setdiff(smp, s)]; tibble(sample = s, nearest = names(which.max(v)), rho = max(v),
                                    second = sort(v, decreasing = TRUE)[2]) })
print(best %>% mutate(across(where(is.numeric), ~ round(.x, 3))), n = 40)
# 背景：同组不同样本之间相关性的分布
cat("\n全部样本对 rho 分位数：\n"); print(round(quantile(r[upper.tri(r)], c(0, .25, .5, .75, .9, .95, 1)), 3))
# 层次聚类，rho 高于阈值的合并为一个"患者"
hc <- hclust(as.dist(1 - r), "average")
th <- 1 - as.numeric(commandArgs(TRUE)[1] %||% "0.97")
grp <- cutree(hc, h = th)
rep_tab <- tibble(sample = names(grp), donor_cluster = paste0("D", grp)) %>%
  mutate(group = st$group[match(sample, st$sample)]) %>% arrange(donor_cluster)
print(rep_tab %>% group_by(donor_cluster) %>% summarise(samples = paste(sample, collapse = ","), group = first(group)), n = 40)
write_csv(rep_tab, file.path(tab, "rev_sc2_replicate_groups.csv"))
