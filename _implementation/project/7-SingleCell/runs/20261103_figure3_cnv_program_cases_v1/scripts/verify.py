from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,subprocess,re,xml.etree.ElementTree as ET
from scipy import sparse
O=Path(__file__).resolve().parents[1];ROOT=O.parents[2];T=O/'tables';C=ROOT/'7-SingleCell/runs/20261102_patient_native_cnv_v3';V=ROOT/'7-SingleCell/runs/20261103_figure3_ecosystem_v3';out=[]
def check(n,ok,e=''):out.append(dict(Check=n,Pass=bool(ok),Evidence=str(e)))
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb')as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
 return h.hexdigest()
for path,digest in json.loads((O/'input_hashes.json').read_text()).items():check('unchanged:'+path,sha(ROOT/path)==digest)
z=np.load(V/'tables/cell_ucell_scores.npz');sc=pd.DataFrame(z['UCell'],index=z['CellKey'],columns=z['programs']);cv=pd.read_csv(V/'tables/ucell_coverage.csv').set_index(['sample','ProgramID']);ident=pd.read_csv(ROOT/'7-SingleCell/runs/20261102_identity_adjudication_v1/tables/adjudicated_cell_annotations.csv.gz',low_memory=False).set_index('CellKey')
for p in ['ICC20','ICC18','ICC24']:
 d=pd.read_csv(T/f'{p}_rows.csv');m=pd.read_csv(T/f'{p}_cnv.csv.gz',header=None).to_numpy();original=sparse.load_npz(C/'matrices'/f'GSE138709__{p}__primary.npz').toarray();ids=pd.read_csv(C/'matrices'/f'GSE138709__{p}__primary_cells.csv').CellKey;order=pd.Index(ids).get_indexer(d.CellKey)
 check(p+':CNV_exact_roundtrip',np.allclose(m,original[order],atol=5e-9,rtol=1e-8),np.max(abs(m-original[order])));check(p+':no_rows_omitted',set(d.CellKey)==set(ids) and len(d)==len(ids) and d.CellKey.is_unique,len(d));check(p+':independent_candidates',np.array_equal(d.StrictCandidate,ident.loc[d.CellKey].StrictCNVCandidateReviewed))
 ss=pd.read_csv(T/f'{p}_programs.csv');expected=sc.loc[d.CellKey,ss.columns].to_numpy();mask=np.array([[cv.loc[(s,pr),'Coverage']>=.7 and cv.loc[(s,pr),'NAvailable']>=10 for pr in ss.columns]for s in d['sample']]);check(p+':scores_by_cellkey',np.allclose(ss.to_numpy()[mask],expected[mask]));check(p+':coverage_mask',np.array_equal(ss.isna(),~mask))
 for v in ['primary','lymphoid','nonlymphoid','window200']:
  mm=sparse.load_npz(C/'matrices'/f'GSE138709__{p}__{v}.npz').toarray()[order];gg=d.Group.copy();ix=d.CleanTumorEpithelial;gg.loc[ix]=d.loc[ix,'SampleLabel']+np.where(d.loc[ix,'StrictCandidate'],' candidates',' other epithelium');meta=pd.read_csv(T/'profile_groups.csv');meta=meta[(meta.Patient==p)&(meta.Variant==v)];ex=np.array([np.median(mm[gg==g],axis=0)for g in meta.Group]);got=pd.read_csv(T/f'{p}_{v}_profiles.csv.gz',header=None).to_numpy();check(p+':'+v+':medians',np.allclose(got,ex,atol=5e-9,rtol=1e-8))
# Independently validate native marker panels and source-availability semantics.
import anndata as ad
atlas=ad.read_h5ad(ROOT/'7-SingleCell/runs/20261102_cluster_guided_annotation_v1/cluster_guided_reviewed_atlas.h5ad',backed='r')
for p in ['ICC20','ICC18','ICC24']:
 d=pd.read_csv(T/f'{p}_rows.csv');mk=pd.read_csv(T/f'{p}_markers.csv');native=atlas[d.CellKey.tolist(),mk.columns.tolist()].X;native=native.toarray() if sparse.issparse(native) else np.asarray(native)
 check(p+':markers_exact',np.allclose(native,mk.to_numpy()));check(p+':markers_available',atlas.varm['native_availability'].loc[mk.columns,d['sample'].unique()].to_numpy().all())
atlas.file.close()
support=pd.read_csv(C/'tables/patient_support.csv');support=support[support.patient.isin(['ICC20','ICC18','ICC24'])];support.to_csv(T/'original_patient_CNV_support.csv',index=False)
for r in pd.read_csv(O/'figure_manifest.csv').itertuples():
 name=r.Stem;pdf=O/'figures'/f'{name}.pdf';svg=(O/'figures'/f'{name}.svg').read_text();check(name+':four_exports',all((O/'figures'/f'{name}.{e}').stat().st_size>1000 for e in ['pdf','svg','png','tiff']))
 tree=ET.fromstring(subprocess.check_output(['/opt/homebrew/bin/pdftotext','-bbox',str(pdf),'-'],text=True));ns={'x':'http://www.w3.org/1999/xhtml'};pages=tree.findall('.//x:page',ns);outside=[];words=[]
 for page in pages:
  w,h=float(page.attrib['width']),float(page.attrib['height'])
  for v in page.findall('.//x:word',ns):
   words.append(v.text);q={k:float(x)for k,x in v.attrib.items()}
   if q['xMin']<0 or q['yMin']<0 or q['xMax']>w or q['yMax']>h:outside.append(v.text)
 check(name+':dimensions',len(pages)==1 and abs(w-r.WidthMM*72/25.4)<1 and abs(h-r.HeightMM*72/25.4)<1);check(name+':text_bounds',len(words)>30 and not outside,outside)
 fonts=subprocess.check_output(['/opt/homebrew/bin/pdffonts',str(pdf)],text=True);(O/'review'/f'{name}_fonts.txt').write_text(fonts);check(name+':embedded_fonts',all('Arial' in f and re.search(r'yes\s+yes\s+yes',f)for f in fonts.strip().splitlines()[2:])and 'Type 3'not in fonts,fonts)
 sizes=[float(q)for q in re.findall(r'font-size:\s*([\d.]+)px',svg)];check(name+':editable_text',bool(re.search(r'<text\b',svg))and len(sizes)>0 and min(sizes)>=5,(min(sizes),max(sizes)))
pd.DataFrame(out).to_csv(O/'review/verification.csv',index=False);print('Checks',len(out),'PASS',sum(x['Pass']for x in out));print(pd.DataFrame([x for x in out if not x['Pass']]).to_string(index=False));assert all(x['Pass']for x in out)
