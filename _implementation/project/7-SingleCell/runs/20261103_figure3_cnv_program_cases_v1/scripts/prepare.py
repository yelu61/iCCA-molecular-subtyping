from pathlib import Path
import numpy as np,pandas as pd,anndata as ad,json,hashlib
from scipy import sparse
from scipy.sparse.linalg import svds
O=Path(__file__).resolve().parents[1];ROOT=O.parents[2];B=ROOT/'7-SingleCell/runs';T=O/'tables';C=B/'20261102_patient_native_cnv_v3';V=B/'20261103_figure3_ecosystem_v3';A=B/'20261102_identity_adjudication_v1'
def sha(p):
 h=hashlib.sha256()
 with open(p,'rb')as f:
  for x in iter(lambda:f.read(8*1024*1024),b''):h.update(x)
 return h.hexdigest()
inputs={}
def track(p):inputs[str(p.relative_to(ROOT))]=sha(p);return p
obs=pd.read_csv(track(V/'tables/cells.csv.gz')).set_index('CellKey');ident=pd.read_csv(track(A/'tables/adjudicated_cell_annotations.csv.gz'),low_memory=False).set_index('CellKey')
z=np.load(track(V/'tables/cell_ucell_scores.npz'));assert z['CellKey'].tolist()==obs.index.tolist()
programs=['Pattern_P2','Pattern_P3','WGCNA_MEcyan','Pattern_P4','WGCNA_MEmagenta','WGCNA_MEdarkturquoise'];markers=['KRT8','KRT18','EPCAM','PTPRC','COL1A1'];score=pd.DataFrame(z['UCell'],index=obs.index,columns=z['programs']);cov=pd.read_csv(track(V/'tables/ucell_coverage.csv')).set_index(['sample','ProgramID'])
atlas=track(B/'20261102_cluster_guided_annotation_v1/cluster_guided_reviewed_atlas.h5ad');a=ad.read_h5ad(atlas,backed='r')
assert a.obs_names.tolist()==obs.index.tolist()
summary=[];dist=[];checks=[];saturation=[];robmeta=[];exclusion=[]
for patient in ['ICC20','ICC18','ICC24']:
 stem=f'GSE138709__{patient}';ids=pd.read_csv(track(C/'matrices'/f'{stem}__primary_cells.csv')).CellKey
 mat=sparse.load_npz(track(C/'matrices'/f'{stem}__primary.npz')).toarray();win=pd.read_csv(track(C/'matrices'/f'{stem}__primary_windows.csv'))
 d=obs.loc[ids].copy();d['StrictCandidate']=ident.loc[ids,'StrictCNVCandidateReviewed'];d['OriginalTarget']=ident.loc[ids,'AnnotationBroad'].isin(['Epithelial','Hepatocyte']);d['CleanTumorEpithelial']=(d.tissue=='Tumor')&(d.AnnotationBroad=='Epithelial')&d.BroadCleanIncluded&~d.library_low_support
 d['SampleLabel']=d['sample'].str.replace(r'ICC_\d+_','',regex=True)
 d['Group']=np.where(~d.OriginalTarget,'Reference controls',np.where(d.AnnotationBroad.eq('Hepatocyte'),d.SampleLabel+' hepatocytes',np.where(d.BroadCleanIncluded,d.SampleLabel+' epithelium',d.SampleLabel+' mixed / review')))
 groups=['Reference controls']+sorted(set(d.Group)-{'Reference controls'});d['GroupIndex']=d.Group.map({x:i for i,x in enumerate(groups)})
 # Deterministic CNV-only ordering; no program values enter sorting.
 centered=mat-mat.mean(axis=0);u,s,vt=svds(centered,k=1,v0=np.ones(min(centered.shape)),random_state=0);pc=u[:,0]*s[0];pc*=np.sign(vt[0,np.argmax(abs(vt[0]))]);d['CNV_PC1']=pc;d['OriginalRow']=range(len(d));d=d.sort_values(['GroupIndex','sample','AnnotationBroad','CNV_PC1'],kind='stable');order=d.OriginalRow.to_numpy();mat=mat[order]
 d['Row']=range(1,len(d)+1);ss=score.loc[d.index,programs].copy();eligible=np.ones(ss.shape,bool)
 for i,(cell,row) in enumerate(d.iterrows()):
  for j,p in enumerate(programs):
   cc=cov.loc[(row['sample'],p)];eligible[i,j]=(cc.Coverage>=.7 and cc.NAvailable>=10)
 ss=ss.where(eligible)
 xx=a[d.index.tolist(),markers].X;xx=xx.toarray()if sparse.issparse(xx)else np.asarray(xx)
 pd.DataFrame(mat).to_csv(T/f'{patient}_cnv.csv.gz',index=False,header=False,float_format='%.9g');ss.to_csv(T/f'{patient}_programs.csv',index=False);pd.DataFrame(xx,columns=markers).to_csv(T/f'{patient}_markers.csv',index=False)
 d.to_csv(T/f'{patient}_rows.csv');win.to_csv(T/f'{patient}_windows.csv',index=False)
 for group,g in d.groupby('Group',sort=False):summary.append(dict(Patient=patient,Group=group,NCells=len(g),NCandidates=int(g.StrictCandidate.sum()),NCleanTumorEpithelial=int(g.CleanTumorEpithelial.sum())))
 for sample,g in d[d.CleanTumorEpithelial].groupby('sample'):
  for p in programs:
   cc=cov.loc[(sample,p)];ok=cc.Coverage>=.7 and cc.NAvailable>=10
   for status,q in g.groupby('StrictCandidate'):
    vals=ss.loc[q.index,p];exclusion.append(dict(Patient=patient,Sample=sample,ProgramID=p,StrictCandidate=bool(status),NCells=len(q),Coverage=cc.Coverage,AvailableGenes=cc.NAvailable,Eligible=bool(ok)))
    for cell,val in vals.items():dist.append(dict(Patient=patient,Sample=sample,ProgramID=p,CellKey=cell,StrictCandidate=bool(status),Score=val,Coverage=cc.Coverage,Eligible=bool(ok)))
 for variant in ['primary','lymphoid','nonlymphoid','window200']:
  vv=sparse.load_npz(track(C/'matrices'/f'{stem}__{variant}.npz')).toarray();vi=pd.read_csv(track(C/'matrices'/f'{stem}__{variant}_cells.csv')).CellKey;ww=pd.read_csv(track(C/'matrices'/f'{stem}__{variant}_windows.csv'));assert vi.tolist()==ids.tolist();vv=vv[order]
  group2=d.Group.copy();group2.loc[d.CleanTumorEpithelial]=d.loc[d.CleanTumorEpithelial,'SampleLabel']+np.where(d.loc[d.CleanTumorEpithelial,'StrictCandidate'],' candidates',' other epithelium')
  gg=['Reference controls']+sorted(set(group2)-{'Reference controls'});out=np.array([np.median(vv[group2.to_numpy()==g],axis=0)for g in gg]);pd.DataFrame(out).to_csv(T/f'{patient}_{variant}_profiles.csv.gz',index=False,header=False,float_format='%.9g');ww.to_csv(T/f'{patient}_{variant}_windows.csv',index=False)
  for k,g in enumerate(gg):robmeta.append(dict(Patient=patient,Variant=variant,Row=k+1,Group=g,NCells=int((group2==g).sum())))
  saturation.append(dict(Patient=patient,Variant=variant,NValues=vv.size,AbsAbove03=int((abs(vv)>.3).sum()),FractionAbove03=float((abs(vv)>.3).mean())))
 checks.extend([dict(Check=f'{patient}:all_saved_rows',Pass=len(d)==len(ids)and set(d.index)==set(ids)),dict(Check=f'{patient}:candidate_count',Pass=int(d.StrictCandidate.sum())=={'ICC20':27,'ICC18':6,'ICC24':9}[patient]),dict(Check=f'{patient}:candidate_clean',Pass=bool(d.loc[d.StrictCandidate,'CleanTumorEpithelial'].all())),dict(Check=f'{patient}:program_join',Pass=np.allclose(ss.to_numpy()[eligible],score.loc[d.index,programs].to_numpy()[eligible])),dict(Check=f'{patient}:all_current_clean_epithelial_present',Pass=set(obs[(obs.patient==patient)&obs.BroadCleanIncluded&(obs.AnnotationBroad=='Epithelial')&~obs.library_low_support].index).issubset(set(ids)))])
 print(patient,len(d),'candidate',d.StrictCandidate.sum(),flush=True)
a.file.close();pd.DataFrame(summary).to_csv(T/'cell_support.csv',index=False);pd.DataFrame(dist).to_csv(T/'program_distributions.csv.gz',index=False);pd.DataFrame(exclusion).to_csv(T/'program_support_coverage.csv',index=False);pd.DataFrame(robmeta).to_csv(T/'profile_groups.csv',index=False);pd.DataFrame(saturation).to_csv(T/'display_saturation.csv',index=False)
q=pd.DataFrame(dist);s=q[q.Eligible].groupby(['Patient','Sample','ProgramID','StrictCandidate']).Score.agg(NCells='size',Median='median',Q25=lambda x:x.quantile(.25),Q75=lambda x:x.quantile(.75)).reset_index();s.to_csv(T/'program_distribution_summary.csv',index=False)
pd.DataFrame(checks).to_csv(O/'review/numerical_checks.csv',index=False);assert all(x['Pass']for x in checks)
(O/'input_hashes.json').write_text(json.dumps(inputs,indent=2));(O/'selection.json').write_text(json.dumps(dict(Patients=['ICC20','ICC18','ICC24'],Programs=programs,Markers=markers,Score='existing UCell; maxRank1500',CNVDisplayLimit=.3,Claims='case description only'),indent=2))
