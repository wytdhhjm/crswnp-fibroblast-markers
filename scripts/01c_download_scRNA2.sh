#!/usr/bin/env bash
# 01c_download_scRNA2.sh
# 第二个单细胞数据集：GSE276503（Xie 2024 Nat Commun, PMID 39438439；IL1B 与鼻黏膜重塑）
# 33 个样本：健康下鼻甲 2、CRSwNP 下鼻甲 9、CRSwNP 中鼻甲 7、CRSwNP 鼻息肉 15；约 1.3 GB
# NCBI 下载慢且经常返回 503，所以：
#   - 同一时间只下一个文件，用 curl -C - 断点续传（中断后从断点继续，不会重复下载）
#   - 503 等错误时等待后重试；每个文件下完都用 gzip -t 校验，坏文件删除重下
#   - 用 .lock 防止同时运行多个下载进程（多个进程同时写同一文件会损坏数据）
# 用法：bash scripts/01c_download_scRNA2.sh（可随时重新运行续传）

cd "$(dirname "$0")" || exit 1
META=$(pwd)/../data/raw/meta/GSE276503.gsm.txt
mkdir -p ../data/raw/GSE276503 && cd ../data/raw/GSE276503 || exit 1
if [ -f .lock ] && kill -0 "$(cat .lock)" 2>/dev/null; then echo "已有下载进程在运行（PID $(cat .lock)），退出"; exit 0; fi
echo $$ > .lock; trap 'rm -f .lock' EXIT

# 下载顺序：健康下鼻甲 → 鼻息肉 → CRSwNP 下鼻甲 → 中鼻甲（复现分析最需要的先下）
grep -E "supplementary_file_[123]" "$META" | sed 's/.*= ftp:/https:/' > urls_all.txt
for pat in HCIT CRSNP GSM8499604 CRSIT CRSMT; do grep -- "$pat" urls_all.txt; done | awk '!s[$0]++' > urls.txt

fetch() {
  local u=$1 f code; f=$(basename "$1")
  gzip -t "$f" 2>/dev/null && return 0
  code=$(curl -s --ssl-no-revoke --max-time 900 --speed-limit 1000 --speed-time 60 -C - -w "%{http_code}" -o "$f.part" "$u")
  # 服务器拒绝续传（416）或返回错误页（503 等）时，丢弃半成品重新开始
  if [ "$code" = "416" ] || head -c 100 "$f.part" 2>/dev/null | grep -qi "<html\|<!DOCTYPE"; then rm -f "$f.part"; return 1; fi
  if gzip -t "$f.part" 2>/dev/null; then mv "$f.part" "$f"; echo "OK   $f $(du -h "$f" | cut -f1)  $(date +%H:%M)"; return 0; fi
  return 1                                         # 未下完（超时/断线），保留 .part 下次续传
}

while read -r u; do
  f=$(basename "$u")
  until fetch "$u"; do
    if [ -f "$f.part" ]; then echo "...  $f 已下载 $(du -h "$f.part" | cut -f1)，继续  $(date +%H:%M)"; else sleep 20; fi
  done
done < urls.txt
echo "完成：$(for f in *_matrix.mtx.gz; do gzip -t "$f" 2>/dev/null && echo x; done | wc -l) / 33 个矩阵文件"
# 内容校验：gzip -t 查不出文件中间被写入的零字节块（并发写入造成），逐个核对矩阵条目数与文件头是否一致
for m in *_matrix.mtx.gz; do
  h=$(zcat "$m" | sed -n 3p | awk '{print $3}'); n=$(zcat "$m" | tail -n +4 | tr -d '\000' | grep -c .)
  [ "$h" = "$n" ] || echo "BAD  $m 条目数 $n ≠ 文件头 $h，请删除后重新运行本脚本"
done
