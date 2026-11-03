# check_manuscript.py — 稿件一致性检查：摘要字数、引用与 bib 是否对应、补充文件是否都在正文中被引用
import re, sys
tex = open(sys.argv[1] if len(sys.argv) > 1 else 'manuscript_plos/plos_main.tex', encoding='utf-8').read()
bib = open('manuscript_plos/references.bib', encoding='utf-8').read()

a = tex[tex.index(r'\section*{Abstract}'):tex.index(r'\clearpage')]
a = re.sub(r'\\gene\{([^}]*)\}', r'\1', a)
a = re.sub(r'\\textbf\{([^}]*)\}', r'\1', a)
a = re.sub(r'\\[a-zA-Z]+\*?(\{[^}]*\})?', '', a)
a = re.sub(r'[$\{\}]', '', a)
print('Abstract words:', len(a.split()) - 1, '(PLOS limit 300)')

cites = set(k.strip() for c in re.findall(r'\\cite\{([^}]*)\}', tex) for k in c.split(','))
keys = set(re.findall(r'@article\{([^,]+),', bib))
print('Cited but missing from bib:', sorted(cites - keys) or 'none')
print('In bib but never cited:', sorted(keys - cites) or 'none')

refs = set(re.findall(r'\\nameref\{(S\d+_\w+)\}', tex))
labs = set(re.findall(r'\\label\{(S\d+_\w+)\}', tex))
print('S files never cited in text:', sorted(labs - refs) or 'none')
print('S files cited but undefined:', sorted(refs - labs) or 'none')
figs = re.findall(r'\\label\{(fig:[^}]*)\}', tex)
print('Figures:', len(figs), '| unreferenced:', [f for f in figs if tex.count('\\ref{' + f + '}') == 0] or 'none')

# 正文（摘要到补充材料之前、去掉图表环境）中的圆括号数量：括号过多读起来生硬
body = tex[tex.index(r'\section*{Abstract}'):tex.index(r'\section*{Supporting information}')]
body = re.sub(r'\\begin\{(figure|table)\}.*?\\end\{\1\}', '', body, flags=re.S)
words = len(re.sub(r'\\[a-zA-Z]+', ' ', body).split())
print('Parentheses in main text:', body.count('('), f'(~{1000 * body.count("(") / words:.1f} per 1,000 words)')
