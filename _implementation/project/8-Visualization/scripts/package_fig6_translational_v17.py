"""Collect independent R-rendered panels, preserve sources, inspect PDFs, write index."""
from pathlib import Path
import csv, hashlib, html, json, os, re, shutil, subprocess
import xml.etree.ElementTree as ET

root=Path(__file__).resolve().parents[2]
base=root/'8-Visualization/refinements'
out=base/'fig6_translational_20261101_v17'
for d in ['figures','supplement','review','tables','qa','provenance']:(out/d).mkdir(exist_ok=True)
def rd(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def wr(p,rows):
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(p).open('w') as f:
        w=csv.DictWriter(f,keys);w.writeheader();w.writerows(rows)
def md5(p):return hashlib.md5(Path(p).read_bytes()).hexdigest()
def rel(p):return str(Path(p).relative_to(root))
rows=rd(out/'new_figure_manifest.csv');carrychecks=[]
def carry(folder,old,new,track,tables=None):
    owner=base/folder
    m=next(r for r in rd(owner/'figure_manifest.csv') if r.get('Panel',r.get('ID'))==old)
    row=dict(Panel=new,Track=track,WidthMM=m['WidthMM'],HeightMM=m['HeightMM'],Reuse='Byte-identical reuse',Origin=rel(owner))
    for key,destdir,suffix in [('PDF','figures' if track=='Main' else 'supplement','.pdf'),('PNG','review','.png')]:
        src=root/m[key];dest=out/destdir/(new+suffix);shutil.copy2(src,dest)
        assert md5(src)==md5(dest)
        if key=='PDF' and 'PDFMD5' in m: assert md5(src)==m['PDFMD5']
        row[key]=rel(dest);carrychecks.append(dict(Panel=new,Kind=key,Origin=rel(src),Output=rel(dest),MD5=md5(dest),Pass=True))
    files=tables or [m.get('SourceTable',m.get('Table'))]
    paths=[]
    for i,p in enumerate(files):
        src=root/p;dest=out/'tables'/(new+('_'+str(i+1) if len(files)>1 else '')+'.csv');shutil.copy2(src,dest);assert md5(src)==md5(dest);paths.append(rel(dest))
        carrychecks.append(dict(Panel=new,Kind='Table',Origin=rel(src),Output=rel(dest),MD5=md5(dest),Pass=True))
    row['SourceTable']=paths[0];row['AdditionalTables']=';'.join(paths[1:]);rows.append(row)

ab='fig6_AB_20261031_v10';abc='fig6_ABC_20261031_v9';bridge='fig6_compact_bridge_20261101_v16'
carry(ab,'Figure6A_clinical_information','Figure6A_clinical_information','Main',[rel(base/ab/'tables/A_displayed_associations.csv'),rel(base/ab/'tables/A_displayed_prediction.csv')])
carry(ab,'Figure6B_measured_duct_increment','Figure6B_measured_duct_increment','Main',[rel(base/ab/'tables/B_displayed_estimates.csv')])
carry(abc,'Figure6C_duct_proportions','Figure6C_duct_proportions','Main',[rel(base/abc/'tables/C_stacked_proportions.csv'),rel(base/abc/'tables/C_patient_source.csv')])
for letter,suffix in [('D','duct_compact_performance'),('E','duct_compact_errors')]:
    carry(bridge,f'Figure6{letter}_{suffix}',f'Figure6{letter}_{suffix}','Main')
carry(bridge,'Figure6F_drug_candidate_uncertainty','Figure6G_drug_candidate_uncertainty','Main')
carry(bridge,'Supplement_duct_reference_sensitivity','Supplement_duct_reference_sensitivity','Supplement')
carry('fig6_duct_anchors_20261031_v11','Duct_anchor_expression_C1_C2_C3','Supplement_duct_anchor_expression','Supplement')
carry('fig6_compact17_expression_20261101_v15','Compact17_expression_C1_C2_C3','Supplement_compact17_expression','Supplement')
rows.sort(key=lambda r:(r['Track']!='Main',r['Panel']))
wr(out/'figure_manifest.csv',rows);wr(out/'qa/carried_file_checks.csv',carrychecks)
for folder in [ab,abc,bridge]:
    for fn in ['README.md','LEGENDS.md','QA.md']:
        src=base/folder/fn
        if src.exists():shutil.copy2(src,out/'provenance'/f'{folder}_{fn}')

# Captions explicitly adapted for the current A/B upper/lower design.
abc_text=(base/abc/'LEGENDS.md').read_text()
abc_text=abc_text[abc_text.index('## A |'):].replace('Left,','Upper,').replace('Right,','Lower,')
abc_text=abc_text.replace('`tables/A_displayed_associations.csv`, `A_displayed_prediction.csv`','`tables/Figure6A_clinical_information_1.csv`, `tables/Figure6A_clinical_information_2.csv`').replace('`tables/B_displayed_estimates.csv`','`tables/Figure6B_measured_duct_increment.csv`').replace('`tables/C_stacked_proportions.csv`, `C_patient_source.csv`, `duct_display_palette.csv`','`tables/Figure6C_duct_proportions_1.csv`, `tables/Figure6C_duct_proportions_2.csv`')
abc_text=re.sub(r'Complete statistical details and sources: \[v8 legend\]\([^\n]+\)\.', 'Inherited statistical definitions remain with the source analysis; no clinical model was refitted.',abc_text)
br=(base/bridge/'LEGENDS.md').read_text()
de=br[br.index('## D.'):br.index('## F.')]
g=br[br.index('## F.'):br.index('## Supplement.')].replace('## F.','## G.').replace(' F ',' G ').replace(' F\n',' G\n').replace('F reuses','G reuses').replace('Figure6F_drug_candidate_uncertainty','Figure6G_drug_candidate_uncertainty')
f='''## F. Patient-level patterns of predicted drug scores

Nine previously selected candidates are shown for all497 Main tumors (C1/C2/C3, n=112/223/162;4,473 unique patient–drug values). The prediction source is GDSC2-trained oncoPredict, not patient PRISM measurements. Each column represents one patient; patient width is equal across subtype blocks. Within each subtype, patients are ordered by Study then Sample, identically for all drugs, without sorting or clustering by drug response. Rows are grouped by the subtype for which each drug was originally nominated; these group labels are not an independent validation result.

The locked score first subtracts the patient's median across198 drug predictions. For display, each drug's score is centered within Study and divided by its pooled centered SD. Blue means a lower relative predicted score for that drug; white denotes zero and rose a higher score. Color saturates at ±2SD (126 values below−2 and92 above+2); all uncapped values are preserved in the source table and distribution supplement. None are deleted. No between-drug potency comparison is supported. No hypothesis testing or significance symbols are added. Candidate selection and these descriptive patterns use the same underlying data; they do not validate treatment response or the17-gene panel as a therapy predictor.

Source: `tables/Figure6F_patient_drug_patterns.csv`; fixed column order in `tables/patient_order.csv`.

'''
supp='''## Supplement: complete candidate distributions

All nine candidates and all497 patients per drug are shown using exactly the same standardized scores as F. Each point is a patient, with deterministic horizontal jitter; boxes show median and interquartile range, whiskers1.5×IQR. All extremes are retained and y axes are not clipped. Panels share a standardized unit and common y limits, not an absolute drug potency scale. No new significance testing or selection of visually stronger representatives was performed. Source: `tables/Supplement_nine_candidate_distributions.csv`.

## Supplement: marker evidence ledger

`tables/marker_evidence_ledger.csv` records the seven duct-related genes and17 compact-development genes, their measured patient/study denominators and roles in classification. Three duct genes (MUC5AC/TFF1/NCAM1) have433 patients; all other listed genes have497. Protein and IHC evidence are explicitly marked as not assessed in this bridge; this does not mean the proteins are absent or antibodies unavailable. The ledger does not claim an RNA/protein concordance or IHC validation experiment.

'''
supp+=br[br.index('## Supplement.'):].replace('tables/Supplement_duct_reference_sensitivity.csv','tables/Supplement_duct_reference_sensitivity.csv')
supp+='\nThe original seven-anchor expression panel uses the same descriptive definitions, with gene-specific denominators. Sources in this package: `tables/Supplement_duct_anchor_expression.csv` and `tables/Supplement_compact17_expression.csv`.\n'
(out/'LEGENDS.md').write_text('# Figure6 | Clinical relevance and translational potential of molecular subtypes\n\nReview candidate, 2026-11-01. Independent panels; no assembled figure or clinical-utility claim. A/B use the v10 upper/lower layout; original estimates and all intervals remain unchanged.\n\n'+abc_text+'\n'+de+f+g+supp)

questions={
'A':('临床预后增量','年龄、性别、分期之外，分型增加多少信息？'),
'B':('实测duct后的增量','已知duct病理后，C3是否仍增加预后信息？'),
'C':('病理与分型的对应','实测、OOF推断与投影duct比例分别如何？'),
'D':('简化识别性能','duct RNA参照与精简流程的识别性能差多少？'),
'E':('识别错误结构','误判发生在哪些亚型、改善发生在哪里？'),
'F':('患者层面治疗假设','候选药物模式和亚型内部异质性如何？'),
'G':('候选效应与不确定性','跨研究效应方向、区间及异质性如何？')}
suppnames={'Supplement_duct_anchor_expression':'七duct-marker表达','Supplement_compact17_expression':'17基因候选表达','Supplement_duct_reference_sensitivity':'duct参照敏感性','Supplement_nine_candidate_distributions':'全部九候选的患者分布'}
cards=[]
for r in rows:
    panel=r['Panel'];main=r['Track']=='Main'
    title,question=questions[panel[7]] if main else (suppnames[panel],'完整数据与分析口径见图注。')
    label=panel[7]+' · ' if main else '补充 · '
    links={k:html.escape(os.path.relpath(root/r[k],out),quote=True) for k in ['PDF','PNG','SourceTable']}
    extras=''.join(f' · <a href="{html.escape(os.path.relpath(root/p,out),quote=True)}">源表{i+2}</a>' for i,p in enumerate(filter(None,r.get('AdditionalTables','').split(';'))))
    cards.append(f'<article><h2>{label}{title}</h2><p>{question}</p><p><a href="{links["PDF"]}">独立PDF</a> · <a href="{links["PNG"]}">PNG</a> · <a href="{links["SourceTable"]}">源表</a>{extras} · {r["WidthMM"]} × {r["HeightMM"]} mm</p><img loading="lazy" src="{links["PNG"]}" alt="{title}" style="width:{r["WidthMM"]}mm"></article>')
atlas=base/'fig6_drug_evidence_20261101_v14';am=rd(atlas/'figure_manifest.csv')
evidence=[('25','候选留一研究稳定性'),('32','组成与增殖代理敏感性'),('33','靶点RNA与蛋白证据'),('34','TIDE跨研究不确定性'),('14','原候选的PRISM细胞系实测'),('35','新版PRISM候选覆盖缺口'),('36','其他共有药物的19线实测')]
links=[];evidence_rows=[]
for prefix,title in evidence:
    r=next(r for r in am if r['ID'].startswith(prefix+'_'));url=html.escape(os.path.relpath(root/r['PDF'],out),quote=True)
    links.append(f'<li><a href="{url}">{title}</a></li>');evidence_rows.append(dict(Topic=title,PDF=r['PDF'],SourceTable=r['Table'],Status='Existing evidence; not newly analyzed'))
wr(out/'supplement_evidence_index.csv',evidence_rows)
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure6 A–G · 临床关联与转化潜力</title><style>body{font:16px/1.6 system-ui,sans-serif;color:#243743;background:#f6f7f8;max-width:1080px;margin:30px auto;padding:0 20px}h1{font-size:27px}h2{font-size:20px}article{background:white;border:1px solid #dce2e7;border-radius:8px;margin:22px 0;padding:20px}img{display:block;max-width:100%;height:auto;margin:18px auto}a{color:#17658b}.note{background:#e8f0f4;padding:16px}</style><h1>Figure6 A–G：临床关联与转化潜力</h1><p>临床增量 → 病理联系 → 简化识别 → 治疗假设</p><p class="note">独立子图审查包。A–E原样沿用；新增F热图，原药物森林图顺延G。D/E评价折内精简流程，不是固定17基因验证；F为GDSC2训练的患者预测，PRISM实测单列于补充。未组装或缩小整张主图。</p><p><a href="README.md">结果与范围</a> · <a href="LEGENDS.md">完整图注</a> · <a href="QA.md">核验</a> · <a href="tables/marker_evidence_ledger.csv">marker证据表</a></p>'''+''.join(cards)+'<article><h2>其余补充证据：按问题索引</h2><ul>'+''.join(links)+'</ul><p>CTRP–PRISM患者预测一致性图暂不纳入：本轮没有可追溯的成对源表。未重做药物模型或补画示意结果。</p></article></html>'
(out/'index.html').write_text(page)

# Real PDF metadata, embedded-font and text-bound checks for all11 deliverables.
checks=[]
for r in rows:
    pdf=root/r['PDF'];name=pdf.stem
    info=subprocess.check_output(['pdfinfo',str(pdf)],text=True)
    fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True)
    bbox=subprocess.check_output(['pdftotext','-bbox',str(pdf),'-'],text=True)
    text=subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True)
    (out/'qa'/f'{name}_pdfinfo.txt').write_text(info);(out/'qa'/f'{name}_fonts.txt').write_text(fonts);(out/'qa'/f'{name}_text.txt').write_text(text)
    pages=list(ET.fromstring(bbox).iter('{http://www.w3.org/1999/xhtml}page'));assert len(pages)==1
    p=pages[0];w,h=float(p.attrib['width']),float(p.attrib['height']);words=list(p.iter('{http://www.w3.org/1999/xhtml}word'))
    outside=[x.text for x in words if float(x.attrib['xMin'])<-.1 or float(x.attrib['yMin'])<-.1 or float(x.attrib['xMax'])>w+.1 or float(x.attrib['yMax'])>h+.1]
    fr=[line.split() for line in fonts.splitlines()[2:] if line.strip()];embedded=bool(fr) and all(f[-5]=='yes' for f in fr)
    sizes=abs(w-float(r['WidthMM'])/25.4*72)<1 and abs(h-float(r['HeightMM'])/25.4*72)<1
    ok=not outside and embedded and sizes and len(words)>20
    checks.append(dict(Panel=name,Pass=ok,FontsEmbedded=embedded,WordsOutside=len(outside),Words=len(words),WidthMM=w/72*25.4,HeightMM=h/72*25.4));assert ok,(name,outside,checks[-1])
wr(out/'qa/pdf_checks.csv',checks)
scan=subprocess.run(['python3','/Users/luye/.agents/skills/nature-figure/scripts/validate_figure.py',str(root/'8-Visualization/scripts/render_fig6_translational_v17.R'),'--backend','r','--json'],capture_output=True,text=True)
(out/'qa/source_scan.json').write_text(scan.stdout)
print(json.dumps(dict(Panels=len(rows),Main=sum(r['Track']=='Main' for r in rows),Supplement=sum(r['Track']=='Supplement' for r in rows),PDFChecksPassed=len(checks),CarriedFiles=len(carrychecks),SourceScan=json.loads(scan.stdout)['summary']),ensure_ascii=False))
