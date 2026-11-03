# 16_geo_search.R
# 回应审稿意见：交代数据集是如何检索和筛选出来的（PRISMA 式流程）
# 用 NCBI E-utilities 检索 GEO DataSets（db=gds），条件写死在脚本里，任何人可重复
# 自动筛选只做客观条件（人类、GSE、表达谱、样本数）；其余逐条按纳入/排除标准人工判定，
# 判定结果和理由写在 results/tables/rev_geo_screening.csv 的 decision / reason 两列
#
# 输出：results/tables/rev_geo_search_raw.csv、rev_geo_screening.csv

suppressPackageStartupMessages({ library(tidyverse); library(xml2) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
tab <- "results/tables"

query <- paste0('("nasal polyps"[MeSH Terms] OR "nasal polyp"[All Fields] OR "nasal polyps"[All Fields] OR ',
                '"chronic rhinosinusitis"[All Fields] OR CRSwNP[All Fields]) AND "Homo sapiens"[Organism] AND gse[Entry Type] AND ',
                '("expression profiling by array"[DataSet Type] OR "expression profiling by high throughput sequencing"[DataSet Type])')
eu <- "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
get <- function(url) {
  f <- tempfile()
  stopifnot(system2("curl", c("-s", "--ssl-no-revoke", "--max-time", "120", "--retry", "3", "-o", f, shQuote(url))) == 0)
  f
}
s <- read_xml(get(paste0(eu, "/esearch.fcgi?db=gds&retmax=1000&term=", URLencode(query, reserved = TRUE))))
ids <- xml_text(xml_find_all(s, "//IdList/Id"))
cat("检索日期：", format(Sys.Date()), "；命中", length(ids), "条\n")

summ <- map_dfr(split(ids, ceiling(seq_along(ids) / 100)), function(chunk) {
  x <- read_xml(get(paste0(eu, "/esummary.fcgi?db=gds&id=", paste(chunk, collapse = ","))))
  map_dfr(xml_find_all(x, "//DocSum"), function(d) {
    it <- function(n) xml_text(xml_find_first(d, paste0("./Item[@Name='", n, "']")))
    tibble(accession = it("Accession"), title = it("title"), type = it("gdsType"),
           n_samples = as.integer(it("n_samples")), platform = it("GPL"), pubdate = it("PDAT"),
           summary = str_trunc(it("summary"), 400))
  })
})
write_csv(summ, file.path(tab, "rev_geo_search_raw.csv"))
cat("其中 GSE 系列：", sum(grepl("^GSE", summ$accession)), "\n")
