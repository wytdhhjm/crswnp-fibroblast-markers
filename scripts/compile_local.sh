#!/usr/bin/env bash
# compile_local.sh —— 在本地用 TeX Live 编译 PLOS 稿件（pdflatex + bibtex，共四遍），输出错误/警告摘要
# 用法：bash scripts/compile_local.sh [manuscript_plos/plos_main.tex]

TEX=${1:-manuscript_plos/plos_main.tex}
cd "$(dirname "$0")/.." || exit 1
DIR=$(dirname "$TEX"); BASE=$(basename "$TEX" .tex)
export PATH="/d/texlive/2026/bin/windows:$PATH"
cd "$DIR" || exit 1

run() { pdflatex -interaction=nonstopmode -halt-on-error -file-line-error "$BASE.tex" > /dev/null 2>&1; }
run; bibtex "$BASE" > "$BASE.bibtex.log" 2>&1; run; run

echo "== PDF: $DIR/$BASE.pdf ($(ls -l "$BASE.pdf" 2>/dev/null | awk '{print $5}') bytes)"
echo "== 错误:";       grep -E "^.*:[0-9]+: |^! " "$BASE.log" | head -20 || true
echo "== 警告:";       grep -E "LaTeX Warning|Package .* Warning|Overfull|Undefined|undefined" "$BASE.log" | sort | uniq -c | head -30 || true
echo "== BibTeX:";     grep -iE "warning|error" "$BASE.bibtex.log" | head -10 || true
grep -o "Output written on .*" "$BASE.log"
