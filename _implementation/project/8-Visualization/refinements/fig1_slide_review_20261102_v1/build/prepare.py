"""Prepare audited inputs and native slide geometry; do not alter source artifacts."""
from pathlib import Path
import csv, json, zipfile, hashlib, subprocess, re
from collections import Counter
import xml.etree.ElementTree as ET

ROOT=Path.cwd()
OUT=ROOT/'8-Visualization/refinements/fig1_slide_review_20261102_v1'
SRC=ROOT/'manuscript/figure_drafts/iCCC_subtypes_figures.pptx'
TAB=ROOT/'8-Visualization/refinements/fig1_sfig1_typography_20261030/source_data'
def read(p):
    with open(p,newline='') as f:return list(csv.DictReader(f))
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
palette={r['FinalSubtype']:r['Color'] for r in read(ROOT/'governance/data/subtype_label_color_dictionary.csv')}
study_palette=read(ROOT/'0-Data/qc/study_color_dictionary.csv')
cols=read(TAB/'Figure1B_columns.csv')
counts=Counter(r['FinalSubtype'] for r in cols)
assert counts=={'C1':112,'C2':223,'C3':162,'C1_C2':41,'C1_C3':39,'C2_C3':9}
assets={}; original_geometry={}; source_hashes={str(SRC.relative_to(ROOT)):sha(SRC)}
ns={'p':'http://schemas.openxmlformats.org/presentationml/2006/main','a':'http://schemas.openxmlformats.org/drawingml/2006/main'}
with zipfile.ZipFile(SRC) as z:
    size=ET.fromstring(z.read('ppt/presentation.xml')).find('p:sldSz',ns).attrib
    for label,index in {'B':8,'C':2,'D':1,'G':3,'H':10}.items():
        content=z.read(f'ppt/media/image{index}.emf'); start=content.index(b'%PDF');end=content.index(b'%%EOF',start)+5
        pdf=OUT/'assets'/f'Figure1{label}_embedded.pdf';pdf.write_bytes(content[start:end])
        svg=pdf.with_suffix('.svg');subprocess.run(['pdftocairo','-svg',str(pdf),str(svg)],check=True)
        subprocess.run(['pdftoppm','-singlefile','-scale-to','2400','-png',str(pdf),str(pdf.with_suffix(''))],check=True)
        # Match the current PowerPoint crop; it removes only the same existing C footer area.
        if label=='C':
            s=svg.read_text();m=re.search(r'viewBox="([^"]+)"',s);v=list(map(float,m.group(1).split()));v[3]*=.84867
            s=s[:m.start(1)]+' '.join(map(str,v))+s[m.end(1):]
            s=re.sub(r'height="[^"]+"',f'height="{v[3]}pt"',s,count=1);svg.write_text(s)
        assets[label]={'svg':str(svg),'pdf':str(pdf),'sha256_pdf':sha(pdf)}
    for pic in ET.fromstring(z.read('ppt/slides/slide1.xml')).findall('.//p:pic',ns):
        rid=pic.find('p:blipFill/a:blip',ns).get('{http://schemas.openxmlformats.org/officeDocument/2006/relationships}embed')
        x=pic.find('p:spPr/a:xfrm',ns)
        original_geometry[rid]={'off':x.find('a:off',ns).attrib,'ext':x.find('a:ext',ns).attrib}

groups=[]; combined=[]
definitions=[(1,'2019Ahn','Original A/B'),(3,'2022Dong','Transcriptomic'),(2,'2022Dong','Proteomic'),(4,'2023Deng','Proteomic')]
for idx,study,modality in definitions:
    path=TAB/f'Figure1E{idx}.csv';d=read(path);samples=read(TAB/f'Figure1E{idx}_samples.csv')
    assert len(samples)==int(d[0]['TotalN'])
    counter=Counter((r['OriginalLabel'],r['FinalSubtype']) for r in samples)
    rows=[]
    for lab in sorted({r['OriginalLabel'] for r in d},key=lambda l:int(next(r['RowOrder'] for r in d if r['OriginalLabel']==l))):
        rr=[next(r for r in d if r['OriginalLabel']==lab and r['FinalSubtype']==c) for c in ('C1','C2','C3')]
        n=[int(r['N']) for r in rr];den=int(rr[0]['RowN']);assert sum(n)==den
        for r in rr:
            assert counter[(lab,r['FinalSubtype'])]==int(r['N'])
            assert abs(int(r['N'])/den-float(r['RowFraction']))<1e-12
            combined.append({'PanelSource':f'Figure1E{idx}',**r})
        rows.append({'label':lab.replace('class: ',''),'n':n,'N':den,'p':[v/den for v in n]})
    groups.append({'study':study,'modality':modality,'rows':rows,'N':len(samples),'available':int(d[0]['AvailableLabelsN']),'columns':[sum(r['n'][j] for r in rows) for j in range(3)]})
    source_hashes[str(path.relative_to(ROOT))]=sha(path)
    source_hashes[str((TAB/f'Figure1E{idx}_samples.csv').relative_to(ROOT))]=sha(TAB/f'Figure1E{idx}_samples.csv')
with open(OUT/'source_data/Figure1E_combined.csv','w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=combined[0]);w.writeheader();w.writerows(combined)
clinical=read(TAB/'Figure1F_compact.csv')
assert len(clinical)==12
source_hashes[str((TAB/'Figure1F_compact.csv').relative_to(ROOT))]=sha(TAB/'Figure1F_compact.csv')

# Geometry is expressed on a 1200 x 1800 review grid and scaled to the original slide.
elements=[]
def rect(panel,name,x,y,w,h,fill,line=None,lw=1):
    elements.append(dict(panel=panel,name=name,kind='rect',x=x,y=y,w=w,h=h,fill=fill,line=line,lw=lw))
def text(panel,name,t,x,y,w,h,size=17,bold=False,align='left',color='#222222'):
    elements.append(dict(panel=panel,name=name,kind='text',text=t,x=x,y=y,w=w,h=h,size=size,bold=bold,align=align,color=color))
def arrow(name,x,y,w,h):elements.append(dict(panel='A',name=name,kind='arrow',x=x,y=y,w=w,h=h,fill='#7B8794'))
def label(c,x,y):text(c,'panel-'+c,c,x,y,32,36,28,True)
def interp(a,b,t):
    aa=[int(a[i:i+2],16) for i in (1,3,5)];bb=[int(b[i:i+2],16) for i in (1,3,5)]
    return '#'+''.join(f'{round(x+(y-x)*t):02X}' for x,y in zip(aa,bb))
def heat(p):return interp('#FFFFFF','#FFFFE0',p*2) if p<=.5 else interp('#FFFFE0',palette['C1'],(p-.5)*2)
def clinical_heat(p):return interp('#F8FAFC','#D5E7F0',p/50) if p<=50 else interp('#D5E7F0','#246085',(p-50)/50)

text('Title','figure-title','Figure 1. Molecular subtypes and clinical associations in iCCA',26,12,1150,48,34,True)
label('A',5,70)
text('A','resource','1,092 iCCA expression samples across 12 studies',40,72,1120,28,21,True)
for nm,x,w in [('discovery',40,210),('clustering',290,210),('principal',540,410),('intermediate',990,170)]:
    rect('A',nm+'-frame',x,113,w,78,'#F5F7FA','#CBD2DA')
text('A','discovery-label','RNA-seq discovery',50,123,190,24,20,True,'center')
text('A','discovery-n','586 samples / 7 studies',50,154,190,22,17,False,'center')
arrow('discovery-to-cluster',260,144,23,15)
text('A','clustering-label','NMF + k-means',300,123,190,24,20,True,'center')
text('A','clustering-note','State discovery',300,154,190,22,17,False,'center')
arrow('cluster-to-main',508,144,23,15)
text('A','principal-title','Principal states: n = 497',550,118,390,24,20,True,'center')
for j,c in enumerate(['C1','C2','C3']):
    rect('A',c+'-count-box',550+j*131,150,123,32,palette[c])
    text('A',c+'-count',f'{c}   n = {counts[c]}',550+j*131,152,123,26,17,True,'center','#FFFFFF' if c=='C1' else '#151515')
text('A','intermediate-title','Intermediate states',995,121,160,23,16.5,True,'center')
text('A','intermediate-n','n = 89',995,149,160,25,19,False,'center')
# Separate branch from clustering to intermediate states, routed above the principal box.
rect('A','intermediate-branch-h',395,104,680,1.2,'#7B8794')
rect('A','intermediate-branch-left',395,104,1.2,9,'#7B8794')
rect('A','intermediate-branch-right',1075,104,1.2,9,'#7B8794')
for nm,x,w in [('projection-input',40,210),('frozen-model',290,210),('projection-result',540,620)]:
    rect('A',nm+'-frame',x,215,w,52,'#FFFFFF','#CBD2DA')
text('A','projection-name','Microarray projection',50,217,190,22,18,True,'center')
text('A','projection-n','506 samples / 5 studies',50,240,190,21,16,False,'center')
arrow('projection-to-model',260,233,23,15)
text('A','model-name','Frozen classifier',300,220,190,23,19,True,'center')
text('A','model-training','Trained on principal states',300,243,190,18,13.5,False,'center')
arrow('model-to-projection',508,233,23,15)
text('A','projection-result','External subtype projection and program preservation',550,227,600,26,19,False,'center')
text('A','evidence-streams','Complementary evidence: multi-omics, clinical outcomes and single-cell localization',40,279,1120,22,16,False,'left','#555555')

placements={'B':[40,337,742,533],'C':[830,320,335,284.4],'D':[830,622,335,248.9],
            'G':[40,1318,540,465.3],'H':[645,1318,520,355.8]}
label('B',5,308);text('B','heatmap-population','586 discovery tumors: 497 principal + 89 intermediate states',40,308,742,24,17,True)
label('C',795,308);label('D',795,610)
label('E',5,891);label('F',682,891);label('G',5,1312);label('H',610,1312)

# E: one common row-normalized matrix, all original observations and missing labels retained.
text('E','E-title','Correspondence with existing classifications',40,892,620,27,20,True)
text('E','E-study-header','Study / classification',40,930,175,23,15.5,True)
text('E','E-original-header','Original',216,930,82,23,15.5,True)
for j,c in enumerate(['C1','C2','C3']):
    rect('E','E-header-'+c,306+94*j,928,90,27,palette[c]);text('E','E-header-text-'+c,c,306+94*j,930,90,23,17,True,'center','#FFFFFF' if c=='C1' else '#151515')
text('E','E-N-header','n',590,930,50,23,16,True,'center')
y=960
for gi,g in enumerate(groups):
    nr=len(g['rows']);mid=y+(nr*21)/2
    text('E',f'E-study-{gi}',g['study'],40,mid-24,170,22,16.5,True)
    text('E',f'E-modality-{gi}',g['modality'],40,mid-3,170,20,14.5)
    text('E',f'E-availability-{gi}',f"Labels {g['available']}/{g['N']}",40,mid+17,170,18,13.5,color='#555555')
    for ri,r in enumerate(g['rows']):
        text('E',f'E-label-{gi}-{ri}',r['label'],216,y,85,21,14.5)
        for j,c in enumerate(['C1','C2','C3']):
            p=r['p'][j];fill='#EFEFEF' if r['label']=='Not reported' else heat(p)
            rect('E',f'E-cell-{gi}-{ri}-{c}',306+94*j,y,90,20,fill)
            text('E',f'E-value-{gi}-{ri}-{c}',f"{r['n'][j]} ({100*p:.0f}%)",306+94*j,y,90,20,14.5,False,'center','#FFFFFF' if p>.76 and r['label']!='Not reported' else '#222222')
        text('E',f'E-N-{gi}-{ri}',str(r['N']),590,y,50,21,14.5,False,'center')
        y+=21
    if gi<3:
        rect('E',f'E-separator-{gi}',216,y+5,424,1,'#D5DADF');y+=14
text('E','E-scale-title','Within original class (%)',40,1280,200,18,13)
for i in range(40):rect('E',f'E-scale-{i}',252+i*5,1285,5.2,7,heat(i/39))
for x,t in [(247,'0'),(345,'50'),(441,'100')]:text('E','E-tick-'+t,t,x,1293,26,17,12,False,'center')
text('E','E-note','Cells: n (row %)',468,1280,180,18,13)

# F: reflow the same descriptive table to a narrower but taller frame.
text('F','F-title','Clinical characteristics',716,892,450,27,20,True)
for j,c in enumerate(['C1','C2','C3']):
    rect('F','F-header-'+c,900+j*88,936,84,27,palette[c]);text('F','F-header-text-'+c,c,900+j*88,938,84,23,17,True,'center','#FFFFFF' if c=='C1' else '#151515')
titles=['Tumor size','Stage III/IV','Vascular invasion','Large-duct type']
notes=['Median [IQR], cm','Observed proportion','Observed proportion','Observed pathology']
for i in range(4):
    rr=[r for r in clinical if int(r['TraitOrder'])==i+1];y=976+i*73
    text('F',f'F-trait-{i}',titles[i],716,y+5,180,22,17,True)
    text('F',f'F-note-{i}',notes[i],716,y+29,180,18,13.5)
    text('F',f'F-n-{i}',f"Available n = {sum(int(r['AvailableN']) for r in rr)}",716,y+48,180,18,13.5,color='#555555')
    for j,c in enumerate(['C1','C2','C3']):
        r=next(r for r in rr if r['FinalSubtype']==c);x=900+j*88
        rect('F',f'F-cell-{i}-{c}',x,y,84,67,'#F2F3F4' if i==0 else clinical_heat(float(r['Value'])))
        text('F',f'F-value-{i}-{c}',r['Primary'],x,y+6,84,25,20,True,'center')
        text('F',f'F-secondary-{i}-{c}',r['Secondary'],x,y+32,84,19,14,False,'center')
        if i==0:text('F',f'F-denominator-{i}-{c}',r['Tertiary'],x,y+50,84,16,12.5,False,'center','#555555')
text('F','F-scale-title','Proportion (%)',716,1280,125,18,13)
for i in range(40):rect('F',f'F-scale-{i}',850+i*3,1285,3.2,7,clinical_heat(i/39*100))
for x,t in [(840,'0'),(898,'50'),(956,'100')]:text('F','F-tick-'+t,t,x,1293,26,17,12,False,'center')
text('F','F-denom-note','n/N: positive / available',985,1280,190,18,12.5)

spec={'slide_size_emu':size,'source':str(SRC),'source_sha256':sha(SRC),'assets':assets,'placements':placements,'elements':elements,'groups':groups,'clinical':clinical,'palette':palette,'study_palette':study_palette,'source_hashes':source_hashes,'original_geometry':original_geometry}
(OUT/'build/spec.json').write_text(json.dumps(spec,indent=2))
(OUT/'review/source_hashes.json').write_text(json.dumps(source_hashes,indent=2))
(OUT/'source_data/Figure1F_compact.csv').write_bytes((TAB/'Figure1F_compact.csv').read_bytes())
(OUT/'review/data_checks.json').write_text(json.dumps({'E_rows':len(combined),'E_class_rows':sum(len(g['rows']) for g in groups),'E_all_counts_match_sample_rows':True,'E_group_populations':[g['N'] for g in groups],'E_missing_label_patients':39,'F_available_denominators':[sum(int(r['AvailableN']) for r in clinical if int(r['TraitOrder'])==i) for i in range(1,5)],'B_state_counts':counts},indent=2))
print(json.dumps({'elements':len(elements),'assets':list(assets),'E_classes':13,'source_sha256':sha(SRC),'size':size}))
