# 11_fix_plos_template.py
# 把 manuscript_plos/plos_main.tex 严格改成 PLOS 官方 LaTeX 模板（v3.8, Apr 2026）的格式
# 正文文字逐字保留，只改：导言区/标题页、摘要（单段、少用缩写）、表格样式、图注格式、
# 图表位置（紧跟首次引用段落）、在 PDF 中显示图片（预览开关 \showfigs）
# 此后 manuscript_plos/plos_main.tex 为稿件主文件（manuscript/main.tex 不再同步）

import io, os, re

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..')
P = os.path.join(ROOT, 'manuscript_plos', 'plos_main.tex')
s = io.open(P, encoding='utf-8').read()

def sub1(old, new):
    global s
    assert s.count(old) == 1, 'not unique / not found: ' + old[:80]
    s = s.replace(old, new)

# ---------- 1. 导言区 + 标题页 + 摘要（模板原样） ----------
head_end = s.index('\\linenumbers\n\\section*{Introduction}')
header = r'''% Template for PLoS
% Version 3.8 Apr 2026
% PLOS One submission. Built from the official PLOS LaTeX template.
\documentclass[10pt,letterpaper]{article}
\usepackage[top=0.85in,left=2.75in,footskip=0.75in]{geometry}

% amsmath and amssymb packages, useful for mathematical formulas and symbols
\usepackage{amsmath,amssymb}

% Use adjustwidth environment to exceed column width (see example table in text)
\usepackage{changepage}

% textcomp package and marvosym package for additional characters
\usepackage{textcomp,marvosym}

% cite package, to clean up citations in the main text. Do not remove.
\usepackage{cite}

% Use nameref to cite supporting information files (see Supporting Information section for more info)
\usepackage{nameref,hyperref}

% line numbers
\usepackage[right]{lineno}

% ligatures disabled
\usepackage[nopatch=eqnum]{microtype}
\DisableLigatures[f]{encoding = *, family = * }

% color can be used to apply background shading to table cells only
\usepackage[table]{xcolor}

% array package and thick rules for tables
\usepackage{array}

% create "+" rule type for thick vertical lines
\newcolumntype{+}{!{\vrule width 2pt}}

% create \thickcline for thick horizontal lines of variable length
\newlength\savedwidth
\newcommand\thickcline[1]{%
  \noalign{\global\savedwidth\arrayrulewidth\global\arrayrulewidth 2pt}%
  \cline{#1}%
  \noalign{\vskip\arrayrulewidth}%
  \noalign{\global\arrayrulewidth\savedwidth}%
}

% \thickhline command for thick horizontal lines that span the table
\newcommand\thickhline{\noalign{\global\savedwidth\arrayrulewidth\global\arrayrulewidth 2pt}%
\hline
\noalign{\global\arrayrulewidth\savedwidth}}

% Text layout
\raggedright
\setlength{\parindent}{0.5cm}
\textwidth 5.25in
\textheight 8.75in

% Bold the 'Figure #' in the caption and separate it from the title/caption with a period
% Captions will be left justified
\usepackage[aboveskip=1pt,labelfont=bf,labelsep=period,justification=raggedright,singlelinecheck=off]{caption}
\renewcommand{\figurename}{Fig}

% Please use the included `plos2025.bst` as your BibTeX style
\bibliographystyle{plos2025}

% Remove brackets from numbering in List of References
\makeatletter
\renewcommand{\@biblabel}[1]{\quad#1.}
\makeatother

% Header and Footer with logo
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

%% Include all user-defined macros below
\newcommand{\gene}[1]{\textit{#1}}
\def\UrlBreaks{\do\/\do-}

% Figures are shown in this PDF for reading convenience only.
% PLOS requires figures to be uploaded as separate files (Fig1.tif ... Fig4.tif) and NOT included in the
% manuscript PDF. Before generating the PDF for submission, change \showfigstrue to \showfigsfalse.
\newif\ifshowfigs
\showfigstrue
\newcommand{\plosfig}[1]{\ifshowfigs\centering\includegraphics[width=\textwidth,height=0.70\textheight,keepaspectratio]{#1}\par\fi}
%% END MACROS SECTION


\begin{document}
\vspace*{0.2in}

% Title must be 250 characters or less.
\begin{flushleft}
{\Large
\textbf{Multi-cohort transcriptomic meta-analysis and single-cell mapping identify fibroblast-enriched \textit{SCG2} and \textit{COL10A1} as robust markers of chronic rhinosinusitis with nasal polyps}
}
\newline
% Insert author names, affiliations and corresponding author email (do not include titles, positions, or degrees).
\\
Yitao Wu\textsuperscript{1*}
\\
\bigskip
\textbf{1} School of Basic Medical Sciences, Xinjiang Medical University, Urumqi, Xinjiang, China
\\
\bigskip

% Use the asterisk to denote corresponding authorship and provide email address in note below.
* wytdhhjm@stu.xjmu.edu.cn

\end{flushleft}
% Please keep the abstract below 300 words
\section*{Abstract}
Chronic rhinosinusitis with nasal polyps is a heterogeneous inflammatory disease of the sinonasal mucosa. Many transcriptomic studies have proposed candidate genes, but most relied on one or two cohorts and seldom resolved which cell type expresses each candidate. We analysed five public bulk transcriptomic cohorts (82 nasal polyps and 61 control mucosa samples; two RNA-sequencing and three microarray platforms) and summarised differential expression with a random-effects meta-analysis of standardised effect sizes. Candidate genes were selected in the two RNA-sequencing cohorts by bootstrap-based LASSO stability selection, random forest and support vector machine recursive feature elimination; the three microarray cohorts were excluded from all selection steps and used for external validation. Cellular sources of the selected genes were mapped in a separate single-cell RNA-sequencing dataset (5 controls, 6 nasal polyps; 32,848 cells) using sample-level statistics. The meta-analysis identified 2,816 genes with large, direction-consistent differences (1,512 up- and 1,304 down-regulated), 96.1\% of which remained significant in leave-one-cohort-out analyses. Up-regulated genes were enriched for adaptive immune activation, T helper cell differentiation and epithelial--mesenchymal transition, and down-regulated genes for amino-acid and fatty-acid metabolism. Four genes were selected by all three algorithms: \gene{PTHLH}, \gene{SCG2}, \gene{COL10A1} and \gene{PDE7A}. A score combining them separated polyps from controls in the three external cohorts (areas under the receiver operating characteristic curve 1.00, 1.00 and 0.96). In single-cell data, \gene{SCG2} and \gene{COL10A1} were predominantly expressed by fibroblasts (78\% and 95\% of positive cells) and were increased in polyp fibroblasts (false discovery rate 0.039 for both) without a detectable change in fibroblast abundance, whereas \gene{PTHLH} localised to epithelial cells and \gene{PDE7A} was broadly expressed. These findings identify four reproducible transcriptomic markers of nasal polyps, indicate that \gene{SCG2} and \gene{COL10A1} mark a remodeling-associated fibroblast state, and refine the cellular attribution of \gene{SCG2}, previously studied in epithelial cells.

\clearpage
\newgeometry{top=0.85in,left=1in,right=1in,footskip=0.75in}
'''
s = header + s[head_end:]

# ---------- 2. 去掉分隔注释 ----------
s = re.sub(r'\n% =+\n', '\n', s)

# ---------- 3. 表 1（模板样式） ----------
t1_old = s[s.index('\\begin{table}[!ht]\n\\centering\n\\caption{{\\bf Bulk'):]
t1_old = t1_old[:t1_old.index('\\end{table}') + len('\\end{table}')]
t1_new = r'''% Place tables after the first paragraph in which they are cited.
\begin{table}[!ht]
\centering
\caption{\textbf{Bulk transcriptomic cohorts included in this study.}}
\small
\begin{tabular}{|l|l|l|r|r|l|}
\hline
\textbf{Dataset} & \textbf{Platform} & \textbf{Polyp / control tissue} & \textbf{Polyps} & \textbf{Controls} & \textbf{Role}\\ \thickhline
GSE136825 & RNA-seq (GPL20301) & Polyp / inferior turbinate & 42 & 28 & Training\\ \hline
GSE179265 & RNA-seq (GPL24676) & Polyp / uncinate process & 17 & 7 & Training\\ \hline
GSE36830 & Affymetrix U133 Plus 2.0 (GPL570) & Polyp / uncinate tissue & 6 & 6 & Validation\\ \hline
GSE23552 & Affymetrix Exon 1.0 ST (GPL5175) & Polyp / nasal mucosa & 10 & 13 & Validation\\ \hline
GSE194282 & Affymetrix Gene 2.1 ST (GPL17692) & Polyp / uncinate tissue & 7 & 7 & Validation\\ \thickhline
Total & & & 82 & 61 & \\ \hline
\end{tabular}
\begin{flushleft} RNA-seq, RNA sequencing. Training cohorts were used for candidate selection and machine learning; validation cohorts were not used in any selection step and contributed only to external validation and to the descriptive meta-analysis.
\end{flushleft}
\label{tab:cohorts}
\end{table}'''
sub1(t1_old, t1_new)

# ---------- 4. 表 2 与 Fig 3：移到首次引用段落之后，并改为模板样式 ----------
t2_old = s[s.index('\\begin{table}[!ht]\n\\centering\n\\caption{{\\bf Hub genes'):]
t2_old = t2_old[:t2_old.index('\\end{table}') + len('\\end{table}')]
s = s.replace(t2_old + '\n', '')
f3_old = s[s.index('\\begin{figure}[!h]\n\\caption{{\\bf Machine learning'):]
f3_old = f3_old[:f3_old.index('\\end{figure}') + len('\\end{figure}')]
s = s.replace(f3_old + '\n', '')

fig3 = r'''% Place figure captions after the first paragraph in which they are cited.
\begin{figure}[!h]
\plosfig{figures/Fig3.png}
\caption{\textbf{Machine learning-based selection and external validation of hub genes.}
A: LASSO cross-validation curve, LASSO selection frequency across 200 stratified bootstraps (dashed line, 0.5), random forest permutation importance (top 20 genes), and overlap of genes selected by the three algorithms; highlighted bars indicate the four hub genes. B: Area under the receiver operating characteristic curve (95\% confidence interval, DeLong method) of each hub gene in the two training cohorts (circles) and three withheld validation cohorts (triangles). C: Receiver operating characteristic curves of the four-gene score.}
\label{fig:ml}
\end{figure}'''
tab2 = r'''\begin{table}[!ht]
\centering
\caption{\textbf{Hub genes selected by all three machine learning algorithms.}}
\small
\begin{tabular}{|l|c|c|c|c|c|}
\hline
 & & & \multicolumn{3}{c|}{\textbf{AUC in validation cohorts}}\\ \cline{4-6}
\textbf{Gene} & \textbf{Pooled \textit{g} (95\% CI)} & \textbf{Meta-analysis FDR} & \textbf{GSE36830} & \textbf{GSE23552} & \textbf{GSE194282}\\ \thickhline
\gene{PTHLH} & 1.94 (1.45--2.44) & $2.5\times10^{-12}$ & 0.94 & 0.81 & 0.92\\ \hline
\gene{PDE7A} & 1.49 (1.02--1.96) & $2.3\times10^{-8}$ & 0.89 & 0.83 & 0.65\\ \hline
\gene{SCG2} & 1.95 (1.21--2.69) & $5.0\times10^{-6}$ & 0.94 & 0.98 & 0.84\\ \hline
\gene{COL10A1} & 1.71 (1.08--2.33) & $2.0\times10^{-6}$ & 1.00 & 0.99 & 0.80\\ \thickhline
Four-gene score & -- & -- & 1.00 & 1.00 & 0.96\\ \hline
\end{tabular}
\begin{flushleft} Pooled effect sizes (Hedges' \textit{g}) and FDR are from the five-cohort random-effects meta-analysis; AUCs are from the three validation cohorts, which were not used for gene selection. The four-gene score is the mean within-cohort \textit{z} score of the four genes. AUC, area under the receiver operating characteristic curve; CI, confidence interval; FDR, false discovery rate.
\end{flushleft}
\label{tab:hub}
\end{table}'''
anchor = '(all FDR $\\leq 5\\times10^{-6}$; Table~\\ref{tab:hub}).\n'
sub1(anchor, anchor + '\n' + fig3 + '\n\n' + tab2 + '\n')

# ---------- 5. 其余图：图注改为 "A: ..." 格式，并加入图片 ----------
def new_fig(label, png, title, panels):
    body = ' '.join(f'{k}: {v}' for k, v in panels)
    return ('\\begin{figure}[!h]\n\\plosfig{figures/' + png + '}\n\\caption{\\textbf{' + title + '}\n' + body +
            '}\n\\label{' + label + '}\n\\end{figure}')

def replace_fig(label, png, title, panels):
    global s
    m = re.search(r'\\begin\{figure\}\[!h\]\n\\caption\{\{\\bf [^\n]*\n[^\n]*\}\n\\label\{' + re.escape(label) + r'\}\n\\end\{figure\}', s)
    assert m, label
    s = s[:m.start()] + new_fig(label, png, title, panels) + s[m.end():]

replace_fig('fig:meta', 'Fig1.png', 'Random-effects meta-analysis of five bulk transcriptomic cohorts.', [
 ('A', "Volcano plot of pooled Hedges' $g$ against meta-analytic false discovery rate (FDR). Coloured points are robust differentially expressed genes (FDR $<0.05$, $|g|\\geq 0.8$, consistent direction in all cohorts); the ten genes with the largest effects in each direction are labelled."),
 ('B', 'Forest plots for the six up-regulated and six down-regulated genes with the largest pooled effects. Circles show cohort-level estimates and diamonds show pooled estimates, with 95\\% confidence intervals.')])
replace_fig('fig:enrich', 'Fig2.png', 'Functional enrichment of the nasal polyp meta-signature.', [
 ('A', 'Top Gene Ontology biological processes and KEGG pathways over-represented among up-regulated (red) and down-regulated (blue) robust differentially expressed genes; numbers indicate gene counts.'),
 ('B', 'Hallmark gene sets significantly enriched by gene set enrichment analysis (adjusted $P<0.05$), ranked by normalised enrichment score.')])
replace_fig('fig:sc', 'Fig4.png', 'Cellular sources of hub genes in single-cell RNA-sequencing data (GSE235711; 5 controls, 6 nasal polyps).', [
 ('A', 'Uniform manifold approximation and projection (UMAP) of 32,848 cells coloured by major cell type, expression of the four hub genes across cell types, and canonical marker genes used for annotation.'),
 ('B', 'Expression of hub genes projected on the UMAP.'),
 ('C', 'Sample-level pseudobulk expression of each hub gene by cell type in controls and nasal polyps; each point is one sample, and only sample--cell type combinations with at least 20 cells are shown.')])

# ---------- 6. 补充材料图注：模板写法 \textbf{...} ----------
s = re.sub(r'\{\\bf ([^}]*)\}', r'\\textbf{\1}', s)

# ---------- 7. 检查 ----------
for bad in ['\\textbf\\newline', '\\toprule', '\\midrule', '\\bottomrule', '\\cmidrule', '\\multirow', 'adjustwidth}{-2.25in', '{\\bf ']:
    assert bad not in s, 'left over: ' + bad
io.open(P, 'w', encoding='utf-8').write(s)
print('ok', len(s))
