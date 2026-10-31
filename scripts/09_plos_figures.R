# 09_plos_figures.R
# 把各子图按 PLOS ONE 要求组合成单个文件（每个编号图一个文件，含 A/B/C 标签）
# 输出：manuscript_plos/figures/Fig1.tif ... Fig4.tif（300 dpi，RGB，LZW 压缩，宽 ≤ 7.5 inch）及同名 PDF 预览
# PLOS 要求：Fig1.tif、Fig2.tif 命名；多面板合并为一张图；图中不含图注

suppressPackageStartupMessages(library(magick))

setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
src <- "results/figures"
out <- "manuscript_plos/figures"
dir.create(out, recursive = TRUE, showWarnings = FALSE)

dpi <- 300
max_w <- round(7.5 * dpi)            # 2250 px：PLOS 最大宽度

# 用 ggsave 输出的 300 dpi PNG（PDF 栅格化会丢失 Symbol 字体，例如坐标轴上的负号）
rd <- function(f) image_read(file.path(src, paste0(f, ".png"))) |> image_background("white") |> image_flatten()
lab <- function(img, l) {
  pad <- image_blank(image_info(img)$width, 70, "white")
  image_annotate(image_append(c(pad, img), stack = TRUE), l, size = 60, weight = 700,
                 gravity = "northwest", location = "+20+0")
}
fitw <- function(img, w = max_w) image_resize(img, paste0(w, "x"))
stackv <- function(...) image_append(do.call(c, list(...)), stack = TRUE)

save_fig <- function(img, name) {
  img <- fitw(img)
  if (image_info(img)$height > round(8.75 * dpi)) img <- image_resize(img, paste0("x", round(8.75 * dpi)))  # 不超过 8.75 inch 高
  image_write(img, file.path(out, paste0(name, ".tif")), format = "tiff", compression = "LZW", density = dpi)
  image_write(img, file.path(out, paste0(name, ".png")), format = "png", density = dpi)
  i <- image_info(img); cat(name, i$width, "x", i$height, "px\n")
}

save_fig(stackv(lab(fitw(rd("Fig1_meta_volcano"), round(max_w * 0.62)), "A"),
                lab(fitw(rd("Fig1_forest_top_genes")), "B")), "Fig1")
save_fig(image_append(c(lab(fitw(rd("Fig2A_ORA_GO_KEGG"), round(max_w * 0.56)), "A"),
                        lab(fitw(rd("Fig2B_GSEA_hallmark"), round(max_w * 0.44)), "B"))), "Fig2")
save_fig(stackv(lab(fitw(rd("Fig3A_ML_selection")), "A"),
                lab(fitw(rd("Fig3B_hub_gene_AUC")), "B"),
                lab(fitw(rd("Fig3C_hub_score_ROC")), "C")), "Fig3")
# Fig4：左列 A、B 上下排列，右列 C，避免整体过高被压窄
left  <- stackv(lab(fitw(rd("Fig4A_scRNA_overview"), round(max_w * 0.6)), "A"),
                lab(fitw(rd("Fig4B_hub_featureplot"), round(max_w * 0.6)), "B"))
right <- lab(fitw(rd("Fig4C_hub_pseudobulk"), round(max_w * 0.4)), "C")
save_fig(image_append(c(left, right), stack = FALSE) |> image_background("white") |> image_flatten(), "Fig4")
save_fig(rd("FigS1_PCA_per_dataset"), "S1_Fig")
cat("完成\n")
