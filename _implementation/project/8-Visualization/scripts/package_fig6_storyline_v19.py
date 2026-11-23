"""Package existing scientific panels and independently verify source-backed additions. No drawing."""
from pathlib import Path
import csv,hashlib,html,json,os,re,shutil,subprocess
import xml.etree.ElementTree as ET
root=Path(__file__).resolve().parents[2];base=root/'8-Visualization/refinements';out=base/'fig6_storyline_20261101_v19';old=base/'fig6_redesign_20261101_v18';atlas=base/'fig6_drug_evidence_20261101_v14'
def rd(p):
 with Path(p).open() as f:return list(csv.DictReader(f))
def wr(p,a):
 keys=list(dict.fromkeys(k for r in a for k in r))
 with Path(p).open('w') as f:w=csv.DictWriter(f,keys);w.writeheader();w.writerows(a)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def rel(p):return str(p.relative_to(root))
rows=rd(out/'new_figure_manifest.csv');copies=[];titles={};questions={}
v18={r['Panel']:r for r in rd(old/'figure_manifest.csv')};v14={r['ID']:r for r in rd(atlas/'figure_manifest.csv')}
def carry(src,new,track,title,question):
 r=dict(src);origin=r.get('Panel',r.get('ID'));r.update(Panel=new,Track=track,Reuse='Byte-identical reuse',OriginPanel=origin)
 if 'Table' in r:r['SourceTable']=r.pop('Table')
 for key,folder in [('PDF','figures' if track=='Main' else 'supplement'),('PNG','review'),('SourceTable','tables')]:
  p=root/r[key];dest=out/folder/(new+p.suffix);shutil.copy2(p,dest);assert sha(p)==sha(dest);copies.append(dict(Panel=new,Type=key,Source=rel(p),Output=rel(dest),SHA256=sha(dest),Pass=True));r[key]=rel(dest)
 extras=[]
 for p in filter(None,r.get('AdditionalTables','').split(';')):
  p=root/p;dest=out/'tables'/p.name;shutil.copy2(p,dest);extras.append(rel(dest))
 r['AdditionalTables']=';'.join(extras);rows.append(r);titles[new]=title;questions[new]=question
for key,title,q in [('Figure6A_clinical_information','A · 临床预后信息','临床变量之外，三亚型及C3分别增加多少信息？'),('Figure6B_measured_duct_increment','B · 实测duct子集','在病理已知的子集中，既有C3专项比较是否仍提示增量？'),('Figure6C_duct_proportions','C · 三亚型与duct表型','三个亚型与实测、OOF及投影duct如何对应？'),('Figure6D_identification_gain','D · 三亚型识别','精简流程相对duct RNA参照改善哪些亚型的识别？'),('Figure6F_drug_heterogeneity','F · 三亚型治疗候选总览','九候选的预测模式与患者异质性如何？')]:carry(v18[key],key,'Main',title,q)
carry(v18['Supplement_confusion_matrices'],'Figure6E_error_structure','Main','E · 完整误判结构','哪些参考亚型被误判，精简流程后误判去向如何变化？')
spec=[('A','Supplement_duct_marker_distributions','duct_expression','七duct marker表达'),('B','Supplement_compact17_distributions','compact17_expression','17基因开发候选表达'),('C','Supplement_reference_sensitivity','reference_sensitivity','四种参照/流程的敏感性'),('D','Figure6E_patient_corrections','patient_transitions','患者纠正与新增错误'),('E','Supplement_drug_distributions','drug_distributions','全部九药物患者分布'),('F','Figure6G_candidate_evidence','drug_uncertainty','效应区间、异质性与排序稳定性')]
for letter,src,suffix,title in spec:carry(v18[src],f'SFig6{letter}_{suffix}','Supplement',f'S{letter} · {title}',title)
for letter,src,suffix,title in [('G','32_candidate_phenotype_sensitivity','phenotype_sensitivity','组成及增殖代理敏感性'),('H','21_TIDE_components','TIDE_components','TIDE分量比较'),('I','34_TIDE_study_meta_uncertainty','TIDE_uncertainty','TIDE跨研究不确定性'),('J','14_PRISM_candidate_observations','PRISM_candidates','旧PRISM候选实测'),('K','35_external_candidate_coverage_gap','PRISM_coverage','外部资源覆盖与缺口'),('L','36_PRISM24_other_shared_drugs','PRISM24_other_drugs','新版PRISM其他共有药物')]:carry(v14[src],f'SFig6{letter}_{suffix}','Supplement',f'S{letter} · {title}',title)
for name,title,q in [('Figure6G_molecular_evidence','G · 九候选的分子证据矩阵','三亚型的相关基因RNA/蛋白证据有哪些一致、矛盾和缺失？'),('Figure6H_C3_hypothesis','H · 最后收束到C3','Ulixertinib预测与MAPK1/MAPK3丰度提供什么可检验假设？')]:titles[name]=title;questions[name]=q
# Copy supporting denominators, summaries and literal legend source references.
for p in (old/'tables').glob('*.csv'):
 if not (out/'tables'/p.name).exists():shutil.copy2(p,out/'tables'/p.name)
for r in v14.values():
 if any(x.get('OriginPanel')==r['ID'] for x in rows):shutil.copy2(root/r['Table'],out/'tables'/Path(r['Table']).name)
for p in [root/'5-ClassifierValidation/review/20261101_drug_evidence_refresh_v1/tables/candidate_target_multiomic_evidence.csv',root/'4-Multiomic/publication/tables/step4_assay_specific_inclusion.csv']:
 shutil.copy2(p,out/'tables'/p.name)
rows.sort(key=lambda r:(r['Track']!='Main',r['Panel']));wr(out/'figure_manifest.csv',rows);wr(out/'qa/reused_files_sha256.csv',copies)
wr(out/'panel_map.csv',[dict(Panel=r['Panel'],Track=r['Track'],Question=questions[r['Panel']],PDF=r['PDF'],SourceTable=r['SourceTable']) for r in rows])
# Independent source checks.
checks=[]
def ck(n,b):checks.append(dict(Check=n,Pass=bool(b)));assert b,n
source=rd(out/'tables/candidate_target_multiomic_evidence.csv');lookup={(r['Drug'],r['Gene'],r['Subtype']):r for r in source}
g=rd(out/'tables/Figure6G_molecular_evidence.csv');h=rd(out/'tables/Figure6H_C3_hypothesis.csv')
def finite(s):
 try:return __import__('math').isfinite(float(s))
 except ValueError:return False
def equal(a,b):return (not finite(a) and not finite(b)) or (finite(a) and finite(b) and abs(float(a)-float(b))<1e-12)
def expected(r):
 s=lookup[r['Drug'],r['Gene'],r['Subtype']];e=s[r['Assay']+'Effect'];f=s[r['Assay']+'FDR']
 state='Unavailable' if not finite(e) or not finite(f) else 'Higher' if float(f)<.05 and float(e)>0 else 'Lower' if float(f)<.05 and float(e)<0 else 'Not significant'
 return equal(r['Effect'],e) and equal(r['FDR'],f) and r['State']==state
ck('G all source effects/FDR and categorical states exact',len(g)==90 and all(expected(r) for r in g))
ck('G complete nine drugs fifteen genes three subtypes two assays',len({r['Drug'] for r in g})==9 and len({r['Gene'] for r in g})==15 and len({(r['Drug'],r['Gene'],r['Subtype'],r['Assay']) for r in g})==90)
ck('H exact twelve molecular observations',len(h)==12 and all(expected(r) for r in h))
pred=rd(out/'tables/Figure6H_patient_predictions.csv');orig={r['Sample']:r for r in rd(old/'tables/Figure6F_drug_heterogeneity.csv') if r['Drug']=='Ulixertinib_1908'}
ck('H all497 predictions and subtype labels exact',len(pred)==497 and all(equal(r['Z'],orig[r['Sample']]['Z']) and r['FinalSubtype']==orig[r['Sample']]['FinalSubtype'] for r in pred))
ck('Eight main and twelve supplementary panels',sum(r['Track']=='Main' for r in rows)==8 and sum(r['Track']=='Supplement' for r in rows)==12)
ck('Every reused artifact byte-identical',all(sha(root/r['Source'])==sha(root/r['Output']) for r in copies))
ck('A B C unchanged',sum(r['Panel'].startswith(('Figure6A','Figure6B','Figure6C')) for r in copies)==9)
ck('Source hashes unchanged',all(sha(root/r['Path'])==r['SHA256'] for r in rd(out/'tables/input_hashes.csv')))
wr(out/'qa/independent_numeric_checks.csv',checks)
# Preserve legacy legends as provenance; authored current legends are supplied separately.
shutil.copy2(old/'LEGENDS.md',out/'LEGENDS_v18_source.md');shutil.copy2(atlas/'LEGENDS.md',out/'LEGENDS_drug_source.md');shutil.copy2(atlas/'LEGENDS_unchanged_v13_source.md',out/'LEGENDS_drug_legacy_source.md')
# PDF text, font and media box inspection.
pdfchecks=[]
for r in rows:
 pdf=root/r['PDF'];fonts=subprocess.check_output(['pdffonts',str(pdf)],text=True);bbox=subprocess.check_output(['pdftotext','-bbox',str(pdf),'-'],text=True);txt=subprocess.check_output(['pdftotext','-layout',str(pdf),'-'],text=True)
 (out/'qa'/f'{r["Panel"]}_text.txt').write_text(txt)
 pages=list(ET.fromstring(bbox).iter('{http://www.w3.org/1999/xhtml}page'));assert len(pages)==1
 p=pages[0];w,h=float(p.attrib['width']),float(p.attrib['height']);words=list(p.iter('{http://www.w3.org/1999/xhtml}word'));outside=[v.text for v in words if float(v.attrib['xMin'])<-.1 or float(v.attrib['yMin'])<-.1 or float(v.attrib['xMax'])>w+.1 or float(v.attrib['yMax'])>h+.1]
 fr=[x.split() for x in fonts.splitlines()[2:] if x.strip()];emb=bool(fr) and all(f[-5]=='yes' for f in fr);sz=abs(w-float(r['WidthMM'])/25.4*72)<1 and abs(h-float(r['HeightMM'])/25.4*72)<1;ok=not outside and emb and sz and len(words)>20
 pdfchecks.append(dict(Panel=r['Panel'],Pass=ok,TextOutside=len(outside),FontsEmbedded=emb,SizeWithinOnePoint=sz));assert ok,(r['Panel'],outside)
wr(out/'qa/pdf_checks.csv',pdfchecks)
scan=subprocess.run(['python3','/Users/luye/.agents/skills/nature-figure/scripts/validate_figure.py',str(root/'8-Visualization/scripts/render_fig6_storyline_v19.R'),'--backend','r','--json'],capture_output=True,text=True);(out/'qa/source_scan.json').write_text(scan.stdout)
# Local review gallery; dimensions shown at intended size, no figure assembly.
cards=[]
for r in rows:
 links={k:html.escape(os.path.relpath(root/r[k],out),quote=True) for k in ['PDF','PNG','SourceTable']}
 cards.append(f'<article><h2>{titles[r["Panel"]]}</h2><p>{questions[r["Panel"]]}</p><p><a href="{links["PDF"]}">PDF</a> · <a href="{links["PNG"]}">PNG</a> · <a href="{links["SourceTable"]}">Source data</a> · {r["WidthMM"]} × {r["HeightMM"]} mm</p><img loading="lazy" src="{links["PNG"]}" style="width:{r["WidthMM"]}mm" alt="{titles[r["Panel"]]}"></article>')
page='''<!doctype html><html lang="zh-CN"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>Figure6 / SFig6 · 三亚型到C3</title><style>body{font:16px/1.65 system-ui;max-width:1080px;margin:28px auto;padding:0 20px;color:#263944;background:#f4f6f7}article{background:white;padding:20px;margin:24px 0;border:1px solid #d9e0e4;border-radius:8px}h1{font-size:28px}h2{font-size:21px}img{max-width:100%;height:auto;display:block;margin:18px auto}a{color:#17658b}.note{background:#e4eef3;padding:16px}</style><h1>Figure6 / SFig6：从C1/C2/C3出发，最后收束C3</h1><p class="note">8张主图子图 + 12张补充子图。临床/病理信息 → 三亚型识别 → 三亚型治疗候选 → 分子证据及缺口 → C3具体假设。独立子图内容版，供后续全篇整合；未进行最终拼版或科学锁定。</p><p><a href="README.md">内容说明</a> · <a href="LEGENDS.md">Figure6图注</a> · <a href="SFIG6_LEGENDS.md">SFig6图注</a> · <a href="QA.md">核验</a> · <a href="../fig6_drug_evidence_20261101_v14/index.html">完整药物证据图集</a></p>'''+''.join(cards)+'</html>'
(out/'index.html').write_text(page)
missing=[]
for v in re.findall(r'(?:href|src)="([^"]+)"',page):
 if not (out/html.unescape(v)).exists():missing.append(v)
(out/'qa/package_links.json').write_text(json.dumps(dict(Missing=missing),ensure_ascii=False));print(json.dumps(dict(Main=8,Supplement=12,PDFChecks=len(pdfchecks),NumericChecks=len(checks),MissingLinks=missing,SourceScan=json.loads(scan.stdout)['summary']),ensure_ascii=False))
