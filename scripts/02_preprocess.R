# 02_preprocess.R
# 数据预处理：读取原始计数/芯片矩阵 -> 统一样本分组 -> 探针/ID 转基因名 -> 保存
#
# 训练集（RNA-seq）：GSE136825、GSE179265
# 验证集（芯片）  ：GSE23552、GSE36830、GSE194282
# 分组统一为 NP（鼻息肉）与 Control（非 CRS 对照黏膜）
#
# 输出：data/processed/*.rds

suppressPackageStartupMessages({
  library(tidyverse)
  library(edgeR)
  library(GEOquery)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
})

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
raw <- "data/raw"
out <- "data/processed"
dir.create(out, showWarnings = FALSE, recursive = TRUE)

# ---------- 工具函数 ----------

# Ensembl -> 基因名的主对照表：GSE179265 原文件自带的 Gene_name 列（GENCODE 注释）
# org.Hs.eg.db 只用于补缺：它对一个 ID 对应多个基因名时取第一个，会出错（如 PRR4 被标为 PRH1）
ens_backup <- local({
  x <- read.delim(file.path(raw, "GSE179265_counts.txt.gz"))
  setNames(x$Gene_name, x$geneID)
})

# 把以 Ensembl ID 为行名的计数矩阵转为基因名；多个 ID 对应同一基因名时，保留总计数最高的那一个
ensembl_to_symbol <- function(counts, ens) {
  ens <- sub("\\..*$", "", ens)                       # 去掉版本号 ENSG...\.5
  sym <- unname(ens_backup[ens])
  miss <- is.na(sym)
  ok <- miss & ens %in% keys(org.Hs.eg.db, "ENSEMBL")
  if (any(ok)) sym[ok] <- mapIds(org.Hs.eg.db, ens[ok], "SYMBOL", "ENSEMBL", multiVals = "first")
  sym <- limma::alias2SymbolTable(sym)                # 统一为 HGNC 官方基因名（与芯片注释一致）
  keep <- !is.na(sym)
  counts <- counts[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(rowSums(counts), decreasing = TRUE)
  counts <- counts[o, , drop = FALSE]; sym <- sym[o]
  counts <- counts[!duplicated(sym), , drop = FALSE]
  rownames(counts) <- sym[!duplicated(sym)]
  counts
}

# 芯片：探针转基因名；多个探针对应同一基因时，保留平均表达最高的探针
probe_to_symbol <- function(expr, annot_pkg) {
  suppressPackageStartupMessages(library(annot_pkg, character.only = TRUE))
  db <- get(annot_pkg)
  sym <- mapIds(db, rownames(expr), "SYMBOL", "PROBEID", multiVals = "first")
  sym <- limma::alias2SymbolTable(sym)
  keep <- !is.na(sym)
  expr <- expr[keep, , drop = FALSE]; sym <- sym[keep]
  o <- order(rowMeans(expr), decreasing = TRUE)
  expr <- expr[o, , drop = FALSE]; sym <- sym[o]
  expr <- expr[!duplicated(sym), , drop = FALSE]
  rownames(expr) <- sym[!duplicated(sym)]
  expr
}

# 芯片数据若未取 log，则做 log2 转换（GEO2R 的判断规则）
log2_if_needed <- function(x) {
  qx <- quantile(x, c(0, 0.25, 0.5, 0.75, 0.99, 1), na.rm = TRUE)
  need <- (qx[5] > 100) || (qx[6] - qx[1] > 50 && qx[2] > 0)
  if (need) { x[x <= 0] <- NA; x <- log2(x) }
  attr(x, "logged_by_us") <- need
  x
}

# 读取 GEO 的 sample 文本元数据（01_download 时保存在 data/raw/meta/）
read_gsm_meta <- function(file) {
  l <- readLines(file, warn = FALSE)
  id <- NA; rows <- list()
  for (x in l) {
    if (startsWith(x, "^SAMPLE")) id <- sub("^\\^SAMPLE = ", "", x)
    else if (startsWith(x, "!Sample_title")) rows[[id]]$title <- sub("^!Sample_title = ", "", x)
    else if (startsWith(x, "!Sample_description")) rows[[id]]$desc <- sub("^!Sample_description = ", "", x)
  }
  tibble(gsm = names(rows),
         title = map_chr(rows, ~ .x$title %||% NA),
         desc  = map_chr(rows, ~ .x$desc %||% NA))
}

# ---------- 训练集 1：GSE136825 ----------
cat("== GSE136825\n")
x <- read.delim(file.path(raw, "GSE136825_counts.txt.gz"), check.names = FALSE)
cnt <- as.matrix(x[, -(1:6)]); colnames(cnt) <- sub("\\.bam$", "", colnames(cnt))
meta <- read_gsm_meta(file.path(raw, "meta/GSE136825.gsm.txt")) %>%
  mutate(prefix = sub(" .*", "", title),
         group = case_when(prefix %in% c("PY", "TT") ~ "NP",            # 鼻息肉组织
                           prefix %in% c("CC", "CCS", "CCY") ~ "Control", # 非 CRS 对照下鼻甲
                           TRUE ~ NA_character_)) %>%                     # SCP/SCT：息肉患者下鼻甲，剔除
  filter(!is.na(group), desc %in% colnames(cnt))
cnt <- cnt[, meta$desc]; colnames(cnt) <- meta$gsm
cnt <- ensembl_to_symbol(cnt, x$Geneid)
gse136825 <- list(counts = cnt,
                  pheno = data.frame(sample = meta$gsm, title = meta$title,
                                     group = factor(meta$group, c("Control", "NP")),
                                     dataset = "GSE136825"))
print(table(gse136825$pheno$group))

# ---------- 训练集 2：GSE179265 ----------
cat("== GSE179265\n")
x <- read.delim(file.path(raw, "GSE179265_counts.txt.gz"), check.names = FALSE)
smp <- grep("^(Control|Jpn_NP|US_NP)", colnames(x), value = TRUE)
cnt <- as.matrix(x[, smp]); rownames(cnt) <- x$geneID
cnt <- ensembl_to_symbol(cnt, x$geneID)
gse179265 <- list(counts = cnt,
                  pheno = data.frame(sample = smp, title = smp,
                                     group = factor(ifelse(grepl("^Control", smp), "Control", "NP"),
                                                    c("Control", "NP")),
                                     population = case_when(grepl("^Jpn", smp) ~ "Japanese",
                                                            grepl("^US", smp) ~ "Caucasian",
                                                            TRUE ~ NA_character_),
                                     dataset = "GSE179265"))
print(table(gse179265$pheno$group))

saveRDS(gse136825, file.path(out, "GSE136825.rds"))
saveRDS(gse179265, file.path(out, "GSE179265.rds"))

# ---------- 验证集：芯片 ----------
read_array <- function(id, annot_pkg, group_fun) {
  cat("==", id, "\n")
  f <- file.path(raw, paste0(id, "_series_matrix.txt.gz"))
  if (!file.exists(f)) { cat("  文件不存在，跳过\n"); return(NULL) }
  es <- getGEO(filename = f, getGPL = FALSE)
  ph <- pData(es)
  grp <- group_fun(ph$title)
  keep <- !is.na(grp)
  ex <- exprs(es)[, keep]
  ex <- log2_if_needed(ex)
  cat("  已做 log2 转换：", attr(ex, "logged_by_us"), "\n")
  ex <- ex[complete.cases(ex), ]
  ex <- probe_to_symbol(ex, annot_pkg)
  res <- list(expr = ex,
              pheno = data.frame(sample = colnames(ex), title = ph$title[keep],
                                 group = factor(grp[keep], c("Control", "NP")), dataset = id))
  print(table(res$pheno$group))
  saveRDS(res, file.path(out, paste0(id, ".rds")))
  res
}

# GSE36830：仅取对照者钩突 vs 鼻息肉（CRSsNP/CRSwNP 患者钩突剔除）
invisible(read_array("GSE36830", "hgu133plus2.db", function(t)
  case_when(grepl("^Control UT", t) ~ "Control", grepl("^CRSwNP polyp", t) ~ "NP")))

# GSE23552：仅取正常黏膜 vs 鼻息肉（炎症黏膜、过敏性鼻炎、混合样本剔除）
invisible(read_array("GSE23552", "huex10sttranscriptcluster.db", function(t)
  case_when(grepl("^NasalMucosa_Normal", t) ~ "Control", grepl("^NasalPolyp_aCRSp", t) ~ "NP")))

# GSE194282：对照钩突 vs 鼻息肉
invisible(read_array("GSE194282", "hugene21sttranscriptcluster.db", function(t)
  case_when(grepl("^Con", t) ~ "Control", grepl("^NP", t) ~ "NP")))

cat("\n预处理完成\n")
