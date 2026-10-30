# 00_install_packages.R
# 安装本项目所需的全部 R 包（使用清华镜像）

options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"),
        BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor",
        timeout = 600)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c(
  "tidyverse", "data.table", "ggplot2", "ggpubr", "pheatmap", "ggrepel",
  "glmnet", "randomForest", "e1071", "pROC", "caret",   # 机器学习与 ROC
  "WGCNA", "rms",                                      # 共表达网络、列线图
  "Seurat"                                             # 单细胞
)
bioc_pkgs <- c(
  "GEOquery", "limma", "sva", "clusterProfiler", "org.Hs.eg.db",
  "enrichplot", "GSVA", "GO.db", "impute", "preprocessCore",
  "DESeq2", "edgeR", "AnnotationDbi"
)

need_bioc <- setdiff(bioc_pkgs, rownames(installed.packages()))
if (length(need_bioc)) BiocManager::install(need_bioc, ask = FALSE, update = FALSE)

need_cran <- setdiff(cran_pkgs, rownames(installed.packages()))
if (length(need_cran)) BiocManager::install(need_cran, ask = FALSE, update = FALSE)

all_pkgs <- c(cran_pkgs, bioc_pkgs)
ok <- sapply(all_pkgs, requireNamespace, quietly = TRUE)
cat("\n===== 安装结果 =====\n")
print(ok)
if (all(ok)) cat("全部安装成功\n") else cat("以下未成功：", names(ok)[!ok], "\n")
