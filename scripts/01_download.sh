#!/usr/bin/env bash
# 01_download.sh
# 从 GEO 下载本研究用到的数据
# 网络到 NCBI 很慢且连接约 75 秒就会被断开，因此：
#   用 HTTPS（TLS 保证数据不被篡改）+ 把文件切成 256KB 小块并行下载，
#   每块单独校验大小、失败重试，最后拼接并做 gzip 完整性检查。
# 用法：bash scripts/01_download.sh

cd "$(dirname "$0")/../data/raw" || exit 1
B=https://ftp.ncbi.nlm.nih.gov/geo/series
CHUNK=262144   # 256 KB
PAR=8          # 同时下载的块数
C="curl -s --ssl-no-revoke --max-time 70"

get_chunk() {  # url start end outfile
  local need=$(( $3 - $2 + 1 ))
  for i in $(seq 1 30); do
    $C -r "$2-$3" -o "$4" "$1"
    [ -f "$4" ] && [ "$(stat -c %s "$4")" -eq "$need" ] && return 0
  done
  return 1
}

dl() {
  local url=$1 out=$2
  if [ -f "$out" ] && gzip -t "$out" 2>/dev/null; then echo "SKIP $out (已完整)"; return; fi
  local size
  size=$($C -I "$url" | tr -d '\r' | awk 'tolower($1)=="content-length:"{print $2}')
  [ -z "$size" ] && { echo "FAIL $out (取不到文件大小)"; return; }
  local tmp="$out.chunks"; mkdir -p "$tmp"
  local n=$(( (size + CHUNK - 1) / CHUNK ))
  for k in $(seq 0 $((n - 1))); do
    local s=$((k * CHUNK)) e=$(( (k + 1) * CHUNK - 1 ))
    [ $e -ge $size ] && e=$((size - 1))
    local f; f=$(printf "%s/%05d" "$tmp" "$k")
    if [ -f "$f" ] && [ "$(stat -c %s "$f")" -eq $((e - s + 1)) ]; then continue; fi
    get_chunk "$url" $s $e "$f" &
    while [ "$(jobs -rp | wc -l)" -ge $PAR ]; do sleep 1; done
  done
  wait
  cat "$tmp"/* > "$out"
  if [ "$(stat -c %s "$out")" -eq "$size" ] && gzip -t "$out" 2>/dev/null; then
    rm -rf "$tmp"; echo "OK   $out $(du -h "$out" | cut -f1)"
  else
    rm -f "$out"; echo "FAIL $out (可重新运行本脚本续传)"
  fi
}

# 训练集：RNA-seq 计数矩阵
dl $B/GSE136nnn/GSE136825/suppl/GSE136825_genecounts_20190903.txt.gz GSE136825_counts.txt.gz
dl $B/GSE179nnn/GSE179265/suppl/GSE179265_readcount_genename.txt.gz  GSE179265_counts.txt.gz
# 验证集：芯片 series matrix
dl $B/GSE36nnn/GSE36830/matrix/GSE36830_series_matrix.txt.gz   GSE36830_series_matrix.txt.gz
dl $B/GSE23nnn/GSE23552/matrix/GSE23552_series_matrix.txt.gz   GSE23552_series_matrix.txt.gz
dl $B/GSE194nnn/GSE194282/matrix/GSE194282_series_matrix.txt.gz GSE194282_series_matrix.txt.gz

# 基因集（富集分析用；KEGG 为当日快照，日期记录在 kegg_download_date.txt）
mkdir -p ../genesets && cd ../genesets || exit 1
[ -s h.all.v2024.1.Hs.symbols.gmt ] || curl -s --max-time 120 -O https://data.broadinstitute.org/gsea-msigdb/msigdb/release/2024.1.Hs/h.all.v2024.1.Hs.symbols.gmt
[ -s kegg_hsa_link.txt ] || { curl -s --max-time 120 -o kegg_hsa_link.txt https://rest.kegg.jp/link/hsa/pathway; date +%F > kegg_download_date.txt; }
[ -s kegg_hsa_list.txt ] || curl -s --max-time 120 -o kegg_hsa_list.txt https://rest.kegg.jp/list/pathway/hsa
