#!/usr/bin/env bash
# 01b_download_scRNA.sh
# 单细胞数据：GSE235711（Immunity 2025；10x 过滤后的计数矩阵，每个样本一个 tar.gz）
# 全部 20 个样本：健康对照 5、CRSsNP 筛窦 5、CRSwNP 筛窦 6、CRSwNP 鼻息肉 6
# 复用 01_download.sh 的分块并行下载函数
# 用法：bash scripts/01b_download_scRNA.sh

cd "$(dirname "$0")" || exit 1
source <(sed -n '/^CHUNK=/,/^}$/p; /^dl() {/,/^}$/p' 01_download.sh)   # 只取函数和参数，不执行下载
mkdir -p ../data/raw/GSE235711 && cd ../data/raw/GSE235711 || exit 1
B=https://ftp.ncbi.nlm.nih.gov/geo/samples/GSM7508nnn

# 先下载对照和鼻息肉（主分析），再下载 CRSsNP、CRSwNP 筛窦（补充分析）
for s in GSM7508289_Control_1 GSM7508290_Control_2 GSM7508291_Control_3 GSM7508292_Control_4 GSM7508293_Control_5 \
         GSM7508303_CRSwNP_NP_1 GSM7508304_CRSwNP_NP_2 GSM7508305_CRSwNP_NP_3 GSM7508306_CRSwNP_NP_4 \
         GSM7508307_CRSwNP_NP_5 GSM7508308_CRSwNP_NP_6 \
         GSM7508294_CRSsNP_Eth_1 GSM7508295_CRSsNP_Eth_2 GSM7508296_CRSsNP_Eth_3 GSM7508297_CRSsNP_Eth_4 \
         GSM7508298_CRSsNP_Eth_5 GSM7508299_CRSwNP_Eth_1 GSM7508300_CRSwNP_Eth_2 GSM7508301_CRSwNP_Eth_3 \
         GSM7508302_CRSwNP_Eth_6; do
  gsm=${s%%_*}
  dl "$B/$gsm/suppl/${s}_filtered.tar.gz" "${s}_filtered.tar.gz"
done
