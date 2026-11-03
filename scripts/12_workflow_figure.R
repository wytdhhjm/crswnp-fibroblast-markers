# 12_workflow_figure.R
# 研究流程图（PLOS ONE 同类文章普遍以 Fig 1 给出研究设计概览；审稿人亦常要求）
# 三列：左 = 发现（训练集），中 = 描述性荟萃分析，右 = 外部验证；下方 = 单细胞定位与稳健性检验
# 箭头只走横/竖方向，且每条箭头都对应真实的数据流向；内容全部来自已完成分析的真实数字
# 输出：results/figures/Fig0_workflow.(pdf|png)

suppressPackageStartupMessages({ library(ggplot2) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))

ink <- "#0b0b0b"; line <- "#898781"
fill_data <- "#e8f1fc"; fill_step <- "#fcfcfb"; fill_val <- "#fdeee6"; fill_sc <- "#e6f6ef"; fill_out <- "#f0efec"

box <- function(x, y, w, h, label, fill, size = 2.5, face = "plain") list(
  annotate("rect", xmin = x - w / 2, xmax = x + w / 2, ymin = y - h / 2, ymax = y + h / 2,
           fill = fill, colour = line, linewidth = 0.35),
  annotate("text", x = x, y = y, label = label, size = size, colour = ink, lineheight = 0.95, fontface = face))
seg <- function(x1, y1, x2, y2, head = TRUE) {
  if (head) annotate("segment", x = x1, y = y1, xend = x2, yend = y2, colour = line, linewidth = 0.45,
                     arrow = grid::arrow(length = grid::unit(1.5, "mm"), type = "closed"))
  else annotate("segment", x = x1, y = y1, xend = x2, yend = y2, colour = line, linewidth = 0.45)
}

xL <- 17; xM <- 50; xR <- 83; W <- 31
y1 <- 88; y2 <- 67; y3 <- 46; y4 <- 26; yC <- 7

p <- ggplot() + coord_cartesian(xlim = c(0, 100), ylim = c(0, 100), expand = FALSE) + theme_void() +
  annotate("text", x = xL, y = 98, label = "Discovery", fontface = "bold", size = 3) +
  annotate("text", x = xM, y = 98, label = "Descriptive meta-analysis", fontface = "bold", size = 3) +
  annotate("text", x = xR, y = 98, label = "External validation", fontface = "bold", size = 3) +
  # 第 1 行：数据
  box(xL, y1, W, 10, "Training cohorts (RNA-seq)\nGSE136825: 42 NP / 28 control\nGSE179265: 17 NP / 7 control", fill_data) +
  box(xM, y1, W, 10, "All five bulk cohorts\n82 NP / 61 control", fill_step) +
  box(xR, y1, W, 10, "Validation cohorts (microarray)\nGSE36830 6/6, GSE23552 10/13\nGSE194282 7/7", fill_val) +
  # 第 2 行
  box(xL, y2, W, 18, "Within-cohort limma-voom\n1,481 candidate genes\n(FDR < 0.05, |log2FC| > 0.58,\nsame direction in both)\n\nLASSO stability selection\nrandom forest | SVM-RFE", fill_data, size = 2.4) +
  box(xM, y2, W, 18, "Within-cohort limma / limma-voom\nt converted to Hedges' g\nDerSimonian-Laird pooling\n(REML + Hartung-Knapp sensitivity)\n\n2,816 DEGs (DL); 190 (REML-HK)\nGO, KEGG, Hallmark GSEA", fill_step, size = 2.4) +
  box(xR, y2, W, 18, "Not used in any\nselection step", fill_val, size = 2.4) +
  # 第 3 行
  box(xL, y3, W, 9, "4 consensus genes\nPTHLH, SCG2, COL10A1, PDE7A", fill_out, face = "bold") +
  box(xR, y3, W, 9, "Per-gene AUC; random-gene and\nclassic-marker baselines", fill_val) +
  # 第 4 行
  box(xL, y4, W, 10, "Specificity and site checks\nCRSsNP, AR, non-polyp CRSwNP mucosa,\nown turbinate; GSE338426, GSE72713", fill_step, size = 2.3) +
  box(xM, y4, W, 10, "scRNA-seq discovery (GSE235711)\n20 samples, 63,059 cells\ncontrol, CRSsNP, CRSwNP ethmoid, NP", fill_sc, size = 2.3) +
  box(xR, y4, W, 10, "scRNA-seq replication (GSE276503)\n29 samples, 195,690 cells\ndonor-level pseudobulk", fill_sc, size = 2.3) +
  # 结论
  box(50, yC, 97, 6, "COL10A1 (and SCG2): fibroblast-derived, higher in polyp fibroblasts than in non-polyp CRSwNP mucosa  |  PTHLH: epithelial",
      fill_out, size = 2.45, face = "bold") +
  # 竖向箭头
  seg(xL, y1 - 5, xL, y2 + 9) + seg(xL, y2 - 9, xL, y3 + 4.5) +
  seg(xM, y1 - 5, xM, y2 + 9) +
  seg(xR, y1 - 5, xR, y2 + 9) + seg(xR, y2 - 9, xR, y3 + 4.5) +
  # 核心基因 → 外部验证（中间列此行为空，横线不穿过任何方框）
  seg(xL + W / 2, y3, xR - W / 2, y3) +
  # 核心基因 → 稳健性检验；核心基因 → 单细胞（经两行之间的空隙）
  seg(xL, y3 - 4.5, xL, y4 + 5) +
  seg(xL, (y3 + y4) / 2 + 0.5, xM, (y3 + y4) / 2 + 0.5, head = FALSE) +
  seg(xM, (y3 + y4) / 2 + 0.5, xM, y4 + 5) +
  # 单细胞 → 定位
  seg(xM + W / 2, y4, xR - W / 2, y4) +
  # → 结论
  seg(xL, y4 - 5, xL, yC + 3) + seg(xR, y4 - 5, xR, yC + 3)

ggsave("results/figures/Fig0_workflow.pdf", p, width = 180, height = 140, units = "mm")
ggsave("results/figures/Fig0_workflow.png", p, width = 180, height = 140, units = "mm", dpi = 300, bg = "white")
cat("完成\n")
