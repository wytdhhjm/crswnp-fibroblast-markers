# 07_fetch_references.R
# 根据已人工核对的 PubMed ID 列表，从 PubMed 获取完整书目信息并生成 BibTeX（manuscript/references.bib）
#   数据集原始文献的 PMID 取自 GEO 页面的 Series_pubmed_id 字段（最可靠来源）
#   其余 PMID 由检索得到，并逐条核对过题目，避免引用编造或错配的文献
# 下载 XML 用 curl（R 自带的网络连接访问 NCBI 会 SSL 报错），解析用 xml2

suppressPackageStartupMessages({ library(xml2); library(tidyverse) })
setwd(file.path(dirname(normalizePath(sub("--file=", "", grep("--file=", commandArgs(FALSE), value = TRUE)))), ".."))
dir.create("manuscript", showWarnings = FALSE)

refs <- tribble(
  ~key,               ~pmid,
  "peng2019",         "31439685",   # GSE136825
  "nakayama2022",     "34863854",   # GSE179265
  "plager2010",       "20625511",   # GSE23552
  "stevens2015",      "26067893",   # GSE36830
  "gse194282",        "35712705",   # GSE194282
  "liao2025",         "40945518",   # GSE235711（单细胞）
  "fokkens2020",      "32077450",
  "bachert2019",      "31543428",
  "il4r2026",         "42375348",
  "hypoxia2026",      "41855680",
  "erstress2026",     "42299589",
  "m2mac2022",        "36466903",
  "emt2026",          "41924281",
  "scg2_2025",        "40946800",
  "fibro2026",        "42809484",
  "korsunsky2019",    "31740819",
  "hao2024",          "37231261",
  "ritchie2015",      "25605792",
  "law2014",          "24485249",
  "robinson2010",     "19910308",
  "wu2021",           "34557778",
  "friedman2010",     "20808728",
  "robin2011",        "21414208",
  "liberzon2015",     "26771021",
  "dersimonian1986",  "3802833",
  "delong1988",       "3203132",
  "hu2025col10a1",    "40826474",
  "szczypka2020",     "32854348",
  "porras2023",       "36745102",
  "hastan2011",       "21605125",
  "haynes2017",       "27896970",
  "squair2021",       "34584091",
  # 2026-11-02 修订新增（PMID 均在 PubMed 核对题目）
  "intHout2014",      "24548571",   # Hartung-Knapp 校正
  "chen2016edger",    "27508061",   # edgeR 准似然流程
  "xie2024",          "39438439",   # GSE276503（单细胞复现数据集）；IDO1+ 成纤维细胞
  "guo2024",          "39614076",   # CXCL12+ 成纤维细胞
  "gse338426",        "42718792",   # GSE338426
  "gse72713",         "27216292"    # GSE72713
)

tmp <- tempfile(fileext = ".xml")
url <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=pubmed&retmode=xml&id=",
              paste(refs$pmid, collapse = ","))
stopifnot(system2("curl", c("-s", "--ssl-no-revoke", "--max-time", "120", "--retry", "3", "-o", tmp, shQuote(url))) == 0)
doc <- read_xml(tmp)

tx <- function(node, xp) { v <- xml_text(xml_find_first(node, xp)); ifelse(is.na(v), "", v) }
# 作者名中的重音字母转为 LaTeX 写法（plos2025.bst + pdfLaTeX 直接输出 UTF-8 重音字母会错位）
tex_accents <- function(s) {
  map <- c("á" = "\\'a", "é" = "\\'e", "í" = "\\'i", "ó" = "\\'o", "ú" = "\\'u",
           "à" = "\\`a", "è" = "\\`e", "ò" = "\\`o", "ä" = '\\"a', "ö" = '\\"o', "ü" = '\\"u',
           "ñ" = "\\~n", "ç" = "\\c c", "ø" = "\\o", "å" = "\\aa", "Ö" = '\\"O', "Ü" = '\\"U')
  for (k in names(map)) s <- gsub(k, paste0("{", map[[k]], "}"), s, fixed = TRUE)
  s
}
esc <- function(s) str_replace_all(s, c("&" = "\\\\&", "%" = "\\\\%", "_" = "\\\\_", "#" = "\\\\#"))

arts <- xml_find_all(doc, "//PubmedArticle")
info <- map_dfr(arts, function(a) {
  au <- xml_find_all(a, ".//AuthorList/Author")
  authors <- map_chr(au, function(x) {
    ln <- tx(x, "./LastName"); fn <- tx(x, "./ForeName"); cn <- tx(x, "./CollectiveName")
    if (nzchar(ln)) paste0(ln, ", ", fn) else paste0("{", cn, "}")
  }) |> tex_accents()
  year <- tx(a, ".//JournalIssue/PubDate/Year")
  if (!nzchar(year)) year <- str_extract(tx(a, ".//JournalIssue/PubDate/MedlineDate"), "\\d{4}")
  pages <- tx(a, ".//Pagination/MedlinePgn")
  if (!nzchar(pages)) pages <- tx(a, ".//ELocationID[@EIdType='pii']")
  tibble(pmid = tx(a, ".//MedlineCitation/PMID"),
         author = paste(authors, collapse = " and "),
         title = str_remove(xml_text(xml_find_first(a, ".//ArticleTitle")), "\\.$") |>
           str_replace_all(c("²" = "$^{2}$", "β" = "$\\\\beta$", "α" = "$\\\\alpha$", "γ" = "$\\\\gamma$")) |>   # 上标²、希腊字母 pdfLaTeX 无法直接处理
           (\(t) if (t == toupper(t)) str_to_sentence(t) else t)(),            # 全大写标题改为句首大写
         journal = tx(a, ".//Journal/ISOAbbreviation"),
         year = year, volume = tx(a, ".//JournalIssue/Volume"), number = tx(a, ".//JournalIssue/Issue"),
         pages = str_replace(pages, "-", "--"),
         # 只取本文自身的 DOI（PubmedData/ArticleIdList）；".//" 会误取到参考文献列表中别的文章的 DOI
         doi = tx(a, "./PubmedData/ArticleIdList/ArticleId[@IdType='doi']"))
})
refs <- left_join(refs, info, by = "pmid")
# PubMed 未收录 DOI 的条目：DOI 经出版商页面/Crossref 人工核对后补入
doi_manual <- c(friedman2010 = "10.18637/jss.v033.i01", delong1988 = "10.2307/2531595")
refs <- refs %>% mutate(doi = ifelse(!nzchar(doi) & key %in% names(doi_manual), doi_manual[key], doi))
# 在线优先发表、尚无卷期页码的条目，页码处标注
refs <- refs %>% mutate(pages = ifelse(key == "fibro2026", "Online ahead of print", pages))
missing <- refs %>% filter(is.na(title) | !nzchar(title))
if (nrow(missing)) stop("未取到：", paste(missing$key, collapse = ", "))

bib <- pmap_chr(refs, function(key, pmid, author, title, journal, year, volume, number, pages, doi) {
  f <- c(author = author, title = paste0("{", title, "}"), journal = journal, year = year,
         volume = volume, number = number, pages = pages, doi = doi, pmid = pmid)
  f <- f[nzchar(f) & f != "{}"]
  f[names(f) %in% c("title", "journal", "doi")] <- esc(f[names(f) %in% c("title", "journal", "doi")])   # plos2025.bst 把 DOI 当普通文本输出，下划线须转义
  paste0("@article{", key, ",\n", paste0("  ", format(names(f)), " = {", f, "}", collapse = ",\n"), "\n}\n")
})
writeLines(bib, "manuscript/references.bib")
file.copy("manuscript/references.bib", "manuscript_plos/references.bib", overwrite = TRUE)
print(refs %>% transmute(key, year, journal, title = str_trunc(title, 70)), n = 50)
cat("\n共", nrow(refs), "条，已写入 manuscript/references.bib\n")
