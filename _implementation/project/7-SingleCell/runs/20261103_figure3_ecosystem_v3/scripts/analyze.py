from pathlib import Path
import json,hashlib,itertools,shutil
import numpy as np,pandas as pd,anndata as ad
from scipy import sparse,stats
from scipy.io import mmwrite
O=Path(__file__).resolve().parents[1];ROOT=O.parents[2];T=O/'tables';V=ROOT/'7-SingleCell/runs/20261103_figure3_cellsource_v2';rng=np.random.default_rng(20261103)
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
 return h.hexdigest()
cfg=json.loads((ROOT/'7-SingleCell/config/current_run.json').read_text());atlas=ROOT/cfg['ActiveAnalysisInputs']['Atlas']
a=ad.read_h5ad(atlas,backed='r');obs=pd.read_csv(V/'tables/cell_analysis_eligibility.csv.gz');assert obs.CellKey.tolist()==a.obs_names.tolist()
for col in ['n_genes_native','pct_mt_native','Figure3CleanEpithelialEligible']:obs[col]=a.obs[col].to_numpy()
obs['UMAP1']=a.obsm['X_umap_harmony'][:,0];obs['UMAP2']=a.obsm['X_umap_harmony'][:,1]
for f in ['programs.csv','membership.csv','coverage.csv']:shutil.copy2(V/'tables'/f,T/f)
defs=pd.read_csv(T/'programs.csv');mem=pd.read_csv(T/'membership.csv');programs=defs.ProgramID.tolist(); genes=np.array(a.var_names);avail=a.varm['native_availability'];sets={p:set(mem.loc[mem.ProgramID==p,'Gene']) for p in programs};common_genes=set(avail.index[avail.all(axis=1)]);common={p:s&common_genes for p,s in sets.items()}
score={v:np.full((len(obs),17),np.nan,dtype=np.float32) for v in ['UCell','UCell1000','UCell3000','UCellCommon']};cov=[];validations=[]
for sample,g in obs.groupby('sample',sort=False):
 ix=g.index.to_numpy();mask=avail[sample].to_numpy(bool);native=genes[mask];x=sparse.csr_matrix(a[ix,mask].X);x.eliminate_zeros();gi={s:i for i,s in enumerate(native)}
 members=[np.array([gi[s] for s in sets[p] if s in gi],dtype=int) for p in programs];cmembers=[np.array([gi[s] for s in common[p] if s in gi],dtype=int) for p in programs]
 for j,p in enumerate(programs):cov.append(dict(sample=sample,ProgramID=p,NDefined=len(sets[p]),NAvailable=len(members[j]),Coverage=len(members[j])/len(sets[p]),NCommon=len(cmembers[j]),CommonCoverage=len(cmembers[j])/len(sets[p]),NUniverse=len(native)))
 for k,oi in enumerate(ix):
  row=x.getrow(k);r=stats.rankdata(-row.data,method='average');zero_rank=len(row.data)+(len(native)-len(row.data)+1)/2
  for variant,mx,mi in [('UCell',1500,members),('UCell1000',1000,members),('UCell3000',3000,members),('UCellCommon',1500,cmembers)]:
   rr=np.full(len(native),min(zero_rank,mx),dtype=float);rr[row.indices]=np.minimum(r,mx)
   for j,inds in enumerate(mi):
    n=len(inds)
    if n:score[variant][oi,j]=1-(rr[inds].sum()-n*(n+1)/2)/(n*mx-n*(n+1)/2)
 # Independent official R reference on 10 cells/library, all programs and rank cutoffs.
 ids=np.unique(np.linspace(0,len(ix)-1,min(10,len(ix)),dtype=int));d=O/'review'/f'ucell_{sample}';d.mkdir(exist_ok=True)
 mmwrite(d/'expression.mtx',x[ids,:].T);pd.Series(native).to_csv(d/'genes.csv',index=False,header=False)
 for v in score:
  pd.DataFrame(score[v][ix[ids]],columns=programs).to_csv(d/f'expected_{v}.csv',index=False)
 validations.append(dict(sample=sample,Directory=str(d.relative_to(O)),NCells=len(ids)))
 print('scored',sample,len(ix),flush=True)
a.file.close();np.savez_compressed(T/'cell_ucell_scores.npz',**score,programs=np.array(programs),CellKey=obs.CellKey.to_numpy(str));obs.to_csv(T/'cells.csv.gz',index=False)
pd.DataFrame(cov).to_csv(T/'ucell_coverage.csv',index=False);pd.DataFrame(validations).to_csv(O/'review/reference_validation_inputs.csv',index=False)
pd.DataFrame([dict(ProgramID=p,Gene=g) for p in programs for g in sorted(common[p])]).to_csv(T/'common_membership.csv',index=False)
cov=pd.DataFrame(cov).set_index(['sample','ProgramID']);oldz=np.load(ROOT/'7-SingleCell/runs/20261102_figure3_program_localization_v1/tables/cell_scores.npz');assert oldz['programs'].tolist()==programs
score['NativeRank']=oldz['full'];score['LogControl']=oldz['log_control'];oldz.close()
# Broad, moderate-resolution and clean epithelial-aware compartments.
rows=[]
for level,flag,typ in [('Broad','BroadIncluded','AnnotationBroad'),('Fine','FineDescriptiveIncluded','FineType'),('Clean','BroadCleanIncluded','AnnotationBroad')]:
 for key,g in obs[obs[flag]].groupby(['cohort','patient','PatientKey','sample','tissue',typ],observed=True):
  meta=dict(zip(['cohort','patient','PatientKey','sample','tissue','CellType'],key));ix=g.index.to_numpy();low=bool(g.library_low_support.any())
  for j,p in enumerate(programs):
   c=cov.loc[(key[3],p)];base=(not low and len(ix)>=20 and c.Coverage>=.7 and c.NAvailable>=10)
   for method,ss in score.items():
    eligible=base and (method!='UCellCommon' or c.NCommon>=10)
    rows.append(dict(**meta,Level=level,ProgramID=p,Method=method,Score=float(np.median(ss[ix,j])),NCells=len(ix),Eligible=eligible,LowSupport=low,Coverage=c.Coverage,NAvailable=c.NAvailable,CommonCoverage=c.CommonCoverage))
lib=pd.DataFrame(rows);lib.to_csv(T/'library_program_scores.csv',index=False)
keys=['Level','Method','cohort','patient','PatientKey','tissue','CellType','ProgramID'];pat=lib[lib.Eligible].groupby(keys,observed=True).agg(Score=('Score','mean'),NCells=('NCells','sum'),NLibraries=('sample','size'),Libraries=('sample',lambda x:';'.join(x)),MinCoverage=('Coverage','min')).reset_index();pat.to_csv(T/'patient_program_scores.csv',index=False)
# Descriptive display includes pooled patient means plus cohort-specific support.
sumkeys=['Level','Method','tissue','CellType','ProgramID']
summary=pat.groupby(sumkeys).agg(Score=('Score','mean'),NPatients=('PatientKey','nunique'),NCells=('NCells','sum')).reset_index()
for cohort in ['GSE151530','GSE138709']:
 q=pat[pat.cohort==cohort].groupby(sumkeys).agg(**{cohort:('PatientKey','nunique')}).reset_index();summary=summary.merge(q,on=sumkeys,how='left')
summary.to_csv(T/'dotplot_summary.csv',index=False)
pat.groupby(['Level','Method','cohort','tissue','CellType','ProgramID']).agg(Score=('Score','mean'),NPatients=('PatientKey','nunique')).reset_index().to_csv(T/'cohort_summary.csv',index=False)
# All same-patient broad contrasts preserve the original testing family.
cr=[];crp=[]
for key,g in pat[pat.Level=='Broad'].groupby(['Method','cohort','tissue','ProgramID']):
 info=dict(zip(['Method','cohort','tissue','ProgramID'],key));ep=g[g.CellType=='Epithelial'].set_index('PatientKey')
 for ct,c in g[g.CellType!='Epithelial'].groupby('CellType'):
  c=c.set_index('PatientKey');ids=ep.index.intersection(c.index);d=ep.loc[ids,'Score'].to_numpy()-c.loc[ids,'Score'].to_numpy();n=len(d)
  if not n:continue
  p=np.nan;lo=hi=np.nan
  if n>=4:
   signs=np.array(list(itertools.product([-1,1],repeat=n)));p=float((np.abs(signs@d/n)>=abs(d.mean())-1e-12).mean());lo,hi=np.quantile(rng.choice(d,(2000,n),replace=True).mean(axis=1),[.025,.975])
  cr.append(dict(**info,Other=ct,NPatients=n,Effect=d.mean(),NPositive=int((d>0).sum()),NNegative=int((d<0).sum()),P=p,Lower95=lo,Upper95=hi))
  for pk,delta in zip(ids,d):crp.append(dict(**info,Other=ct,PatientKey=pk,patient=ep.loc[pk,'patient'],EpithelialScore=ep.loc[pk,'Score'],OtherScore=c.loc[pk,'Score'],Difference=delta))
cr=pd.DataFrame(cr);cr['FDR']=np.nan
for _,g in cr.dropna(subset=['P']).groupby(['Method','cohort','tissue']):cr.loc[g.index,'FDR']=stats.false_discovery_control(g.P)
cr.to_csv(T/'source_contrasts.csv',index=False);pd.DataFrame(crp).to_csv(T/'source_patient_pairs.csv',index=False)
# Method concordance: retain changes rather than cherry-pick source directions.
base=cr[cr.Method=='UCell'];compar=[]
for method in ['NativeRank','LogControl','UCell1000','UCell3000','UCellCommon']:
 j=base.merge(cr[cr.Method==method],on=['cohort','tissue','ProgramID','Other'],suffixes=('_UCell','_other'));j['SameDirection']=np.sign(j.Effect_UCell)==np.sign(j.Effect_other);j['Comparator']=method;compar.append(j)
pd.concat(compar).to_csv(T/'method_concordance.csv',index=False)
# Captured proportions: all QC cells including Uncertain in denominator, non-low libraries only.
props=[]
for key,g in obs[~obs.library_low_support].groupby(['cohort','patient','PatientKey','sample','tissue']):
 for ct in ['Fibroblast','Myeloid','T_NK']:
  props.append(dict(zip(['cohort','patient','PatientKey','sample','tissue'],key),CellType=ct,Fraction=float(g.AnnotationBroad.eq(ct).mean()),Numerator=int(g.AnnotationBroad.eq(ct).sum()),Denominator=len(g)))
props=pd.DataFrame(props);props.to_csv(T/'library_captured_fractions.csv',index=False);prop=props.groupby(['cohort','patient','PatientKey','tissue','CellType']).Fraction.mean().reset_index();prop.to_csv(T/'patient_captured_fractions.csv',index=False)
assocspec=[('C1_metabolic_IFN','Pattern_P2','Myeloid','WGCNA_MEdarkturquoise'),('C2_epithelial_IFN','Pattern_P3','Myeloid','WGCNA_MEdarkturquoise'),('C3_epithelial_IFN','Pattern_P4','Myeloid','WGCNA_MEdarkturquoise'),('C3_epithelial_Tcycle','Pattern_P4','T_NK','WGCNA_MEmagenta'),('C3_epithelial_fibroblast','Pattern_P4','Fibroblast','Fraction'),('C2_epithelial_myeloid','Pattern_P3','Myeloid','Fraction')]
pd.DataFrame(assocspec,columns=['Association','EpithelialProgram','TMECellType','TMEProgram']).to_csv(T/'association_specification.csv',index=False)
def adjusted_r(x,y,cohort):
 z=pd.get_dummies(pd.Series(cohort),dtype=float).to_numpy();rx=stats.rankdata(x);ry=stats.rankdata(y);rx=rx-z@np.linalg.lstsq(z,rx,rcond=None)[0];ry=ry-z@np.linalg.lstsq(z,ry,rcond=None)[0]
 if np.std(rx)<1e-10 or np.std(ry)<1e-10:return np.nan,rx,ry
 return float(np.corrcoef(rx,ry)[0,1]),rx,ry
ar=[];points=[]
for method in ['UCell','UCellCommon','NativeRank','LogControl']:
 for name,xp,ct,yp in assocspec:
  ep=pat[(pat.Level=='Clean')&(pat.Method==method)&(pat.tissue=='Tumor')&(pat.CellType=='Epithelial')&(pat.ProgramID==xp)][['cohort','patient','PatientKey','Score']].rename(columns={'Score':'X'})
  if yp=='Fraction':yy=prop[(prop.tissue=='Tumor')&(prop.CellType==ct)][['PatientKey','Fraction']].rename(columns={'Fraction':'Y'})
  else:yy=pat[(pat.Level=='Clean')&(pat.Method==method)&(pat.tissue=='Tumor')&(pat.CellType==ct)&(pat.ProgramID==yp)][['PatientKey','Score']].rename(columns={'Score':'Y'})
  d=ep.merge(yy,on='PatientKey');n=len(d);r=pp=np.nan;loo=[];lo=hi=np.nan
  if n>=6:
   r,rx,ry=adjusted_r(d.X,d.Y,d.cohort);den=np.linalg.norm(rx)*np.linalg.norm(ry)
   if np.isfinite(r):
    permutations=np.tile(ry,(9999,1));
    for cohort,gg in d.reset_index(drop=True).groupby('cohort'):
     idx=gg.index.to_numpy();permutations[:,idx]=np.stack([rng.permutation(ry[idx]) for _ in range(9999)])
    pp=(1+int((np.abs(permutations@rx/den)>=abs(r)-1e-12).sum()))/10000
    for drop in range(n):
     q=d.drop(d.index[drop]);loo.append(adjusted_r(q.X,q.Y,q.cohort)[0])
  ar.append(dict(Association=name,Method=method,NPatients=n,N151530=int((d.cohort=='GSE151530').sum()),N138709=int((d.cohort=='GSE138709').sum()),AdjustedR=r,P=pp,LOOMin=np.nanmin(loo) if loo else np.nan,LOOMax=np.nanmax(loo) if loo else np.nan,XProgram=xp,YCellType=ct,YProgram=yp))
  for _,row in d.iterrows():points.append(dict(Association=name,Method=method,**row.to_dict()))
ar=pd.DataFrame(ar);ar['FDR']=np.nan
for _,g in ar.dropna(subset=['P']).groupby('Method'):ar.loc[g.index,'FDR']=stats.false_discovery_control(g.P)
ar.to_csv(T/'patient_associations.csv',index=False);pd.DataFrame(points).to_csv(T/'patient_association_points.csv',index=False)
# Within clean epithelial cells, library + depth-adjusted program correlation (descriptive only).
cellcorr=[]
for pk,g in obs[(obs.tissue=='Tumor')&obs.BroadCleanIncluded&(obs.AnnotationBroad=='Epithelial')&~obs.library_low_support].groupby('PatientKey'):
 for p1,p2 in itertools.combinations(['Pattern_P2','Pattern_P3','Pattern_P4'],2):
  valid=g['sample'].map(lambda s:cov.loc[(s,p1),'Coverage']>=.7 and cov.loc[(s,p2),'Coverage']>=.7);q=g[valid];ix=q.index.to_numpy()
  if len(q)<50:continue
  design=np.column_stack([pd.get_dummies(q['sample'],dtype=float).to_numpy(),stats.rankdata(q.n_genes_native),stats.rankdata(q.pct_mt_native)])
  x=score['UCell'][ix,programs.index(p1)];y=score['UCell'][ix,programs.index(p2)];rx=stats.rankdata(x);ry=stats.rankdata(y);rx-=design@np.linalg.lstsq(design,rx,rcond=None)[0];ry-=design@np.linalg.lstsq(design,ry,rcond=None)[0]
  cellcorr.append(dict(PatientKey=pk,patient=q.patient.iloc[0],cohort=q.cohort.iloc[0],Program1=p1,Program2=p2,NCells=len(q),RawR=stats.spearmanr(x,y).statistic,AdjustedR=np.corrcoef(rx,ry)[0,1],SharedGenes=len(sets[p1]&sets[p2])))
pd.DataFrame(cellcorr).to_csv(T/'within_epithelial_program_correlations.csv',index=False)
# Score-depth diagnostics and coverage-aware feature data (all atlas rows retained).
dep=[]
for key,g in obs[obs.BroadIncluded&~obs.library_low_support].groupby(['sample','AnnotationBroad']):
 if len(g)<20:continue
 for method in ['UCell','NativeRank','LogControl']:
  for j,p in enumerate(programs):
   if cov.loc[(key[0],p),'Coverage']<.7:continue
   dep.append(dict(sample=key[0],CellType=key[1],Method=method,ProgramID=p,NCells=len(g),DepthR=stats.spearmanr(score[method][g.index,j],g.n_genes_native).statistic))
pd.DataFrame(dep).to_csv(T/'score_depth_diagnostics.csv',index=False)
selected=['Pattern_P2','Pattern_P3','Pattern_P4','WGCNA_MEcyan','WGCNA_MEmagenta','WGCNA_MEdarkturquoise'];f=obs[['CellKey','cohort','patient','sample','tissue','UMAP1','UMAP2','library_low_support','BroadIncluded']].copy()
for p in selected:
 c=obs['sample'].map(lambda s:cov.loc[(s,p),'Coverage']>=.7 and cov.loc[(s,p),'NAvailable']>=10);f[p]=np.where(c&~obs.library_low_support,score['UCell'][:,programs.index(p)],np.nan)
f.to_csv(T/'featureplot_cells.csv.gz',index=False)
# Actual joint score availability; not a positive/active classification.
coverage=[]
for method in ['UCell']:
 for level in ['Broad','Clean']:
  for typ in ['Epithelial','Myeloid','T_NK']:
   q=pat[(pat.Method==method)&(pat.Level==level)&(pat.CellType==typ)&(pat.tissue=='Tumor')]
   for name,ps in [('Patterns',['Pattern_P2','Pattern_P3','Pattern_P4']),('Modules',['WGCNA_MElightgreen','WGCNA_MEcyan','WGCNA_MEmagenta'])]:
    z=q[q.ProgramID.isin(ps)].groupby('PatientKey').ProgramID.nunique();ids=z[z==3].index
    coverage.append(dict(Level=level,CellType=typ,ProgramGroup=name,NPatients=len(ids),PatientKeys=';'.join(ids)))
pd.DataFrame(coverage).to_csv(T/'coavailability.csv',index=False)
inputs={'atlas':{'path':str(atlas.relative_to(ROOT)),'sha256':sha(atlas)},'membership':{'path':str((T/'membership.csv').relative_to(ROOT)),'sha256':sha(T/'membership.csv')},'prior_run':V.name,'UCell_reference':json.loads((O/'review/UCell_reference.json').read_text())}
(O/'review/input_manifest.json').write_text(json.dumps(inputs,indent=2)+'\n')
(O/'analysis_complete.json').write_text(json.dumps({'cells':len(obs),'patients':obs.PatientKey.nunique(),'libraries':obs['sample'].nunique(),'programs':17,'patient_score_rows':len(pat),'library_score_rows':len(lib),'score_min':float(score['UCell'].min()),'score_max':float(score['UCell'].max())},indent=2)+'\n')
print(ar[ar.Method=='UCell'].to_string(index=False));print('analysis complete',flush=True)
