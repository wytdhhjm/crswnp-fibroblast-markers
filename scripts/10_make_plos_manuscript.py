# 10_make_plos_manuscript.py
# 由 manuscript/main.tex（通用版）生成 PLOS ONE 投稿版 manuscript_plos/plos_main.tex
# 只改格式，不改科学内容：
#   - 使用 PLOS 官方模板（v3.8, Apr 2026）的导言区与标题页
#   - 章节不编号；"Figure" -> "Fig"
#   - 正文不嵌入图片：图注保留在首次引用段落之后，图片另行上传（Fig1.tif ...）
#   - 声明类信息（资助、利益冲突、作者贡献、数据可用性）不放正文，在投稿系统中填写；
#     AI 使用声明按 PLOS 政策放在 Methods 专门小节
#   - 补充材料：S1 Fig、S1–S4 Table 作为单独文件上传，正文末尾只列图注
#   - 参考文献：plos2025.bst（Vancouver）
# 用法：python scripts/10_make_plos_manuscript.py

import io, re, os, shutil

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
src = io.open(os.path.join(ROOT, 'manuscript', 'main.tex'), encoding='utf-8').read()
OUT = os.path.join(ROOT, 'manuscript_plos')
os.makedirs(OUT, exist_ok=True)

def between(s, a, b):
    i = s.index(a); j = s.index(b, i)
    return s[i + len(a):j]

title = re.search(r'\\title\{\\textbf\{(.*?)\}\}\n', src, re.S).group(1)
abstract = between(src, '\\begin{abstract}\n', '\\end{abstract}').strip()
body = between(src, '\\section{Introduction}', '\\section*{Declarations}')
body = '\\section{Introduction}' + body

# ---------- 声明部分：取出 AI 声明正文，其余交给投稿系统 ----------
decl = between(src, '\\section*{Declarations}', '\\clearpage\n\\bibliographystyle')
ai_text = re.search(r'\\textbf\{Use of AI tools:\}\s*(.*?)\n\n', decl, re.S).group(1).strip()
ack_text = re.search(r'\\textbf\{Acknowledgements:\}\s*(.*?)\n', decl, re.S).group(1).strip()

# ---------- 图：去掉嵌入的图片，只保留图注 ----------
fig_caps = {
 'fig:meta': ('Random-effects meta-analysis of five bulk transcriptomic cohorts.', None),
}
def strip_figure(m):
    blk = m.group(0)
    if '\\caption*' in blk:                      # Fig 4 续页：并入 Fig 4 图注
        return ''
    cap = re.search(r'\\caption\{(.*)\}\s*\n\\label\{([^}]*)\}', blk, re.S)
    text, label = cap.group(1).strip(), cap.group(2)
    first, rest = text.split('. ', 1) if '. ' in text else (text, '')
    if label == 'fig:sc':
        rest = rest + ' (C) Sample-level pseudobulk expression of each hub gene by cell type in controls and nasal polyps; each point is one sample, and only sample--cell type combinations with at least 20 cells are shown.'
    return ('\\begin{figure}[!h]\n\\caption{{\\bf ' + first.rstrip('.') + '.}\n' + rest.strip() +
            '}\n\\label{' + label + '}\n\\end{figure}')
body = re.sub(r'\\begin\{figure\}\[p\].*?\\end\{figure\}', strip_figure, body, flags=re.S)

# ---------- 表：标题加粗首句，宽表用 adjustwidth ----------
def fix_table(m):
    blk = m.group(0)
    blk = blk.replace('\\begin{table}[ht]', '\\begin{table}[!ht]')
    blk = re.sub(r'\\caption\{([^.]*\.)\s*', lambda c: '\\caption{{\\bf ' + c.group(1) + '} ', blk, count=1)
    if 'tabularx' in blk:
        blk = blk.replace('\\begin{tabularx}{\\textwidth}{l l X r r l}', '\\begin{tabular}{l l l r r l}')
        blk = blk.replace('\\end{tabularx}', '\\end{tabular}')
    blk = blk.replace('\\begin{tabular}', '\\begin{adjustwidth}{-2.25in}{0in}\n\\begin{tabular}', 1)
    blk = blk.replace('\\end{tabular}', '\\end{tabular}\n\\end{adjustwidth}', 1)
    return blk
body = re.sub(r'\\begin\{table\}\[ht\].*?\\end\{table\}', fix_table, body, flags=re.S)

# ---------- 章节不编号、引用格式 ----------
body = re.sub(r'\\(sub)?section\{', lambda m: '\\' + (m.group(1) or '') + 'section*{', body)
body = body.replace('\\section*{Conclusions}', '\\section*{Conclusion}')
body = re.sub(r'Figure~\\ref\{([^}]*)\}([A-C,]*)', lambda m: 'Fig~\\ref{' + m.group(1) + '}' + m.group(2), body)
body = body.replace('Supplementary Table~S1', '\\nameref{S1_Table}').replace('Supplementary Table~S2', '\\nameref{S2_Table}')
body = body.replace('Supplementary Table~S3', '\\nameref{S3_Table}').replace('Supplementary Table~S4', '\\nameref{S4_Table}')
body = body.replace('Supplementary Figure~S1', '\\nameref{S1_Fig}')
body = body.replace('(Supplementary data on GitHub)', '(\\nameref{S5_Table})')
body = body.replace('\\citep{', '\\cite{')

# AI 使用声明：PLOS 要求放在 Methods 的专门小节
ai_sec = '\\subsection*{Use of artificial intelligence tools}\n' + ai_text + '\n\n'
body = body.replace('\\section*{Results}', ai_sec + '\\section*{Results}', 1)

# ---------- 补充材料图注 ----------
supp = r'''
\section*{Supporting information}

\paragraph*{S1 Fig.}
\label{S1_Fig}
{\bf Principal component analysis of each cohort.} PCA of the 2,000 most variable genes within each bulk cohort; points are coloured by diagnosis.

\paragraph*{S1 Table.}
\label{S1_Table}
{\bf Differentially expressed genes per cohort.} Number of genes with FDR $<0.05$ and $|\log_2\mathrm{FC}|>1$ in each bulk cohort.

\paragraph*{S2 Table.}
\label{S2_Table}
{\bf Hub genes after excluding possibly duplicated samples.} limma-voom results in GSE136825 after excluding polyp samples TT18 and TT21 (40 polyps vs.\ 28 controls).

\paragraph*{S3 Table.}
\label{S3_Table}
{\bf Threshold sensitivity of hub-gene selection.} Genes selected by all three algorithms across nine combinations of LASSO selection frequency and random forest/SVM-RFE rank thresholds.

\paragraph*{S4 Table.}
\label{S4_Table}
{\bf Adjustment for fibroblast content.} Association of each hub gene with polyp status in each cohort before and after adjustment for a bulk fibroblast marker score.

\paragraph*{S5 Table.}
\label{S5_Table}
{\bf Single-cell pseudobulk tests.} Sample-level comparisons of hub-gene expression between nasal polyps and controls for every cell type.

\paragraph*{S6 Table.}
\label{S6_Table}
{\bf Meta-analysis results for all genes.} Per-gene pooled effect sizes, confidence intervals, heterogeneity statistics and FDR from the random-effects meta-analysis, with robust DEG and leave-one-cohort-out flags.
'''

# ---------- 组装 ----------
tpl_pre = r'''% PLOS ONE submission version, generated by scripts/10_make_plos_manuscript.py
% Template: PLOS LaTeX template v3.8 (Apr 2026)
\documentclass[10pt,letterpaper]{article}
\usepackage[top=0.85in,left=2.75in,footskip=0.75in]{geometry}
\usepackage{amsmath,amssymb}
\usepackage{changepage}
\usepackage{textcomp,marvosym}
\usepackage{cite}
\usepackage{nameref,hyperref}
\usepackage[right]{lineno}
\usepackage[nopatch=eqnum]{microtype}
\DisableLigatures[f]{encoding = *, family = * }
\usepackage[table]{xcolor}
\usepackage{array}
\usepackage{booktabs}
\usepackage{multirow}
\newcolumntype{+}{!{\vrule width 2pt}}
\newlength\savedwidth
\newcommand\thickcline[1]{%
  \noalign{\global\savedwidth\arrayrulewidth\global\arrayrulewidth 2pt}%
  \cline{#1}%
  \noalign{\vskip\arrayrulewidth}%
  \noalign{\global\arrayrulewidth\savedwidth}%
}
\newcommand\thickhline{\noalign{\global\savedwidth\arrayrulewidth\global\arrayrulewidth 2pt}%
\hline
\noalign{\global\arrayrulewidth\savedwidth}}
\raggedright
\setlength{\parindent}{0.5cm}
\textwidth 5.25in
\textheight 8.75in
\usepackage[aboveskip=1pt,labelfont=bf,labelsep=period,justification=raggedright,singlelinecheck=off]{caption}
\renewcommand{\figurename}{Fig}
\bibliographystyle{plos2025}
\makeatletter
\renewcommand{\@biblabel}[1]{\quad#1.}
\makeatother
\usepackage{lastpage,fancyhdr,graphicx}
\usepackage{epstopdf}
\pagestyle{fancy}
\fancyhf{}
\rfoot{\thepage/\pageref{LastPage}}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrule}{\hrule height 2pt \vspace{2mm}}
\fancyheadoffset[L]{2.25in}
\fancyfootoffset[L]{2.25in}
\lfoot{\today}
\newcommand{\gene}[1]{\textit{#1}}
\def\UrlBreaks{\do\/\do-}

\begin{document}
\vspace*{0.2in}

\begin{flushleft}
{\Large
\textbf\newline{''' + title + r'''}
}
\newline
\\
Yitao Wu\textsuperscript{1*}
\\
\bigskip
\textbf{1} School of Basic Medical Sciences, Xinjiang Medical University, Urumqi, Xinjiang, China
\\
\bigskip
* wytdhhjm@stu.xjmu.edu.cn
\end{flushleft}

\section*{Abstract}
''' + abstract + r'''

\linenumbers
'''

tpl_post = r'''
\section*{Acknowledgments}
''' + ack_text + r'''

\nolinenumbers
\bibliography{references}

\end{document}
'''

# 补充材料（图、表）从正文中移除 —— 已在 body 截止于 Declarations 之前，故无需额外处理
tex = tpl_pre + body.strip() + '\n' + supp + tpl_post
io.open(os.path.join(OUT, 'plos_main.tex'), 'w', encoding='utf-8').write(tex)

# 附带文件
tpl_dir = os.path.join(ROOT, 'overleaf_check', 'plos_tpl')
shutil.copy(os.path.join(tpl_dir, 'plos2025.bst'), OUT)
shutil.copy(os.path.join(ROOT, 'manuscript', 'references.bib'), OUT)
print('written', os.path.join(OUT, 'plos_main.tex'), len(tex), 'chars')
