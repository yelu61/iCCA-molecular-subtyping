from pathlib import Path
import ast,json,hashlib,sys
import anndata as ad,numpy as np,pandas as pd
from scipy.stats import rankdata
ROOT=Path(__file__).resolve().parents[4];OUT=Path(__file__).resolve().parents[1];T=OUT/'tables'
if (OUT/'scoring_complete.json').exists():raise SystemExit('Complete output exists')
source=ROOT/'7-SingleCell/runs/20261102_two_cohort_repair_v1/r3_annotation/broad_annotated_atlas.h5ad'
cand=ROOT/'7-SingleCell/review/20261102_figure3_score_candidates';m=pd.read_csv(cand/'candidate_gene_membership.csv');defs=pd.read_csv(cand/'candidate_programs.csv');defs.to_csv(T/'programs.csv',index=False);m.to_csv(T/'membership.csv',index=False)
programs=defs.ProgramID.tolist();genes={p:set(m.loc[m.ProgramID==p,'Gene']) for p in programs};core={p:set(m.loc[(m.ProgramID==p)&((m.Affinity>=.6) if p.startswith('Pattern') else (m.kME>=.5)),'Gene']) for p in programs}
astree=ast.parse((ROOT/'7-SingleCell/scripts/repair_20261102/r3a_annotation_evidence.py').read_text());panels=next(ast.literal_eval(n.value) for n in astree.body if isinstance(n,ast.Assign) and any(isinstance(t,ast.Name) and t.id=='PANELS' for t in n.targets));markers=set(sum(panels.values(),[]));nomarkers={p:g-markers for p,g in genes.items()};alltarget=set.union(*genes.values())
a=ad.read_h5ad(source,backed='r');obs=a.obs.copy();v=a.varm['native_availability'];obs['PatientKey']=obs.cohort.astype(str)+'|'+obs.patient.astype(str);obs['cell_index']=np.arange(a.n_obs)
for axis,key in [('harmony','X_umap_harmony'),('native','X_umap_unintegrated')]:
 for j in range(2):obs[f'{axis}_{j+1}']=a.obsm[key][:,j]
obs.to_csv(T/'cells.csv.gz',index=False)
variants=['full','core','marker_removed','common','log_control'];arrays={q:np.full((a.n_obs,len(programs)),np.nan,np.float32) for q in variants};detect=np.full_like(arrays['full'],np.nan);coverage=[];control_rows=[]
common={}
for (cohort,tissue),g in obs.groupby(['cohort','tissue'],observed=True):common[(str(cohort),str(tissue))]=set(v.index[v[g['sample'].unique()].all(axis=1)])
for sample,g in obs.groupby('sample',observed=True,sort=False):
 idx=g.cell_index.to_numpy();available=v[str(sample)].to_numpy(bool);native=np.array(a.var_names)[available];native_set=set(native);G=len(native);X=a[idx,available].X.tocsr();nativeindex={s:i for i,s in enumerate(native)};cohort=str(g.cohort.iloc[0]);tissue=str(g.tissue.iloc[0]);cg=common[(cohort,tissue)];ci=np.array([j for j,s in enumerate(native) if s in cg]);CG=len(ci)
 means=np.asarray(X.mean(axis=0)).ravel();bins=np.minimum(19,((rankdata(means,method='average')-1)/G*20).astype(int));rng=np.random.default_rng(1711);bin_ctrl={}
 for b in range(20):
  choices=np.array([j for j,s in enumerate(native) if bins[j]==b and s not in alltarget],dtype=int);bin_ctrl[b]=np.sort(rng.choice(choices,min(50,len(choices)),replace=False))
  for j in bin_ctrl[b]:control_rows.append(dict(sample=sample,bin=b,Gene=native[j]))
 sels={};cms={};weights={}
 for k,p in enumerate(programs):
  for variant,ss in [('full',genes[p]),('core',core[p]),('marker_removed',nomarkers[p]),('common',genes[p]&cg)]:
   used=ss&native_set;sels[(variant,k)]=np.array([nativeindex[s] for s in sorted(used)],dtype=int)
   defined=len(genes[p]) if variant=='common' else len(ss)
   coverage.append(dict(sample=sample,cohort=cohort,tissue=tissue,ProgramID=p,Variant=variant,NDefined=defined,NAvailable=len(used),Coverage=len(used)/defined if defined else 0,NUniverse=CG if variant=='common' else G))
  cms[k]=np.flatnonzero(np.isin(native[ci],list(genes[p])));fb=bins[sels[('full',k)]];weights[k]=np.bincount(fb,minlength=20)/max(1,len(fb))
 for start in range(0,len(idx),256):
  stop=min(start+256,len(idx));xx=X[start:stop].toarray();rr=(rankdata(xx,axis=1,method='average')-(G+1)/2)/(G-1);cr=(rankdata(xx[:,ci],axis=1,method='average')-(CG+1)/2)/(CG-1)
  bc=np.column_stack([xx[:,bin_ctrl[b]].mean(axis=1) if len(bin_ctrl[b]) else np.full(stop-start,np.nan) for b in range(20)])
  for k,p in enumerate(programs):
   for variant in ['full','core','marker_removed']:
    jj=sels[(variant,k)]
    if len(jj):arrays[variant][idx[start:stop],k]=rr[:,jj].mean(axis=1)
   if len(cms[k]):arrays['common'][idx[start:stop],k]=cr[:,cms[k]].mean(axis=1)
   jj=sels[('full',k)]
   if len(jj):
    usedbins=np.flatnonzero(weights[k]);arrays['log_control'][idx[start:stop],k]=xx[:,jj].mean(axis=1)-(bc[:,usedbins]*weights[k][usedbins]).sum(axis=1);detect[idx[start:stop],k]=(xx[:,jj]>0).mean(axis=1)
 print('scored',sample,len(idx),'native genes',G,flush=True)
pd.DataFrame(coverage).to_csv(T/'coverage.csv',index=False);pd.DataFrame(control_rows).to_csv(T/'control_genes.csv',index=False);pd.DataFrame([{'ProgramID':p,'NMarkerOverlap':len(genes[p]&markers),'MarkerGenes':';'.join(sorted(genes[p]&markers))} for p in programs]).to_csv(T/'annotation_overlap.csv',index=False)
np.savez_compressed(T/'cell_scores.npz',**arrays,detection=detect,programs=np.array(programs));a.file.close()
(OUT/'scoring_complete.json').write_text(json.dumps({'cells':len(obs),'programs':len(programs),'variants':variants,'membership_sha256':hashlib.sha256((T/'membership.csv').read_bytes()).hexdigest(),'contract_sha256':hashlib.sha256((OUT/'analysis_contract.md').read_bytes()).hexdigest()},indent=2))
