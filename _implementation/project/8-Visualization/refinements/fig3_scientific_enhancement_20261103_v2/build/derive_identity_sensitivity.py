"""Broad-versus-Clean eligibility sensitivity from saved patient scores only.
Preserve full context x fixed-program x comparator grid, and use identical
patients on each side of every plotted comparison. No new p values or scoring.
"""
from pathlib import Path
import hashlib, itertools, json
import numpy as np
import pandas as pd
O=Path(__file__).resolve().parents[1]
T=O/'source_data'
p=pd.read_csv(T/'patient_program_scores.csv')
l=pd.read_csv(T/'library_program_scores.csv')
cells=pd.read_csv(T/'cells.csv.gz')
reg=pd.read_csv(T/'shared_program_registry.csv')
keys=['Level','Method','cohort','patient','PatientKey','tissue','CellType','ProgramID']
calc=l[l.Eligible].groupby(keys).Score.mean(); saved=p.set_index(keys).Score
assert calc.index.equals(saved.index)
err=float(np.max(np.abs(calc-saved)));assert err<1e-12
p=p[p.Method.eq('UCell') & p.Level.isin(['Broad','Clean'])].copy()
assert not p.duplicated(keys).any()
contexts=cells[['cohort','tissue']].drop_duplicates().sort_values(['cohort','tissue'])
others=['Hepatocyte','T_NK','B_Plasma','Myeloid','Fibroblast','Endothelial']
records=[];pts=[]
for ctx in contexts.itertuples():
 for pr,other in itertools.product(reg.ProgramID,others):
  q=p[(p.cohort==ctx.cohort)&(p.tissue==ctx.tissue)&(p.ProgramID==pr)]
  meta=dict(cohort=ctx.cohort,tissue=ctx.tissue,ProgramID=pr,Other=other)
  maps={}
  for level,ct in itertools.product(['Broad','Clean'],['Epithelial',other]):
   z=q[(q.Level==level)&(q.CellType==ct)].set_index('PatientKey')
   maps[level,ct]=z
  broad=set(maps['Broad','Epithelial'].index)&set(maps['Broad',other].index)
  clean=set(maps['Clean','Epithelial'].index)&set(maps['Clean',other].index)
  matched=sorted(broad&clean);union=sorted(broad|clean)
  db=[];dc=[];nlibraries_same=0
  for patient in union:
   row=dict(**meta,PatientKey=patient,Matched=patient in matched)
   for (level,ct),z in maps.items():
    role='Epithelial' if ct=='Epithelial' else 'Comparator'
    present=patient in z.index
    row[level+role+'Score']=float(z.loc[patient,'Score']) if present else np.nan
    row[level+role+'Libraries']=str(z.loc[patient,'Libraries']) if present else ''
   row['BroadDifference']=row['BroadEpithelialScore']-row['BroadComparatorScore']
   row['CleanDifference']=row['CleanEpithelialScore']-row['CleanComparatorScore']
   row['SameEligibleLibraries']=all(set(row['Broad'+role+'Libraries'].split(';'))==set(row['Clean'+role+'Libraries'].split(';')) for role in ['Epithelial','Comparator']) if row['Matched'] else False
   row['PatientStatus']='matched' if row['Matched'] else ('broad_only_pair' if patient in broad else 'clean_only_pair')
   if row['Matched']:
    db.append(row['BroadDifference']);dc.append(row['CleanDifference']);nlibraries_same+=row['SameEligibleLibraries']
   pts.append(row)
  n=len(matched);eb=float(np.mean(db)) if n else np.nan;ec=float(np.mean(dc)) if n else np.nan
  status='qualified_tumor_n4' if ctx.tissue=='Tumor' and n>=4 else ('descriptive_adjacent' if ctx.tissue=='Adjacent' and n else ('insufficient_matched_patients' if n else ('no_pairs_in_either' if not union else 'no_common_patients')))
  records.append(dict(**meta,NBroadPairs=len(broad),NCleanPairs=len(clean),NMatched=n,NBroadOnly=len(broad-clean),NCleanOnly=len(clean-broad),NMatchedSameEligibleLibraries=nlibraries_same,BroadEffect=eb,CleanEffect=ec,EffectChange=ec-eb,SameDirection=(bool(np.sign(eb)==np.sign(ec)) if n else None),Status=status,MatchedPatients=';'.join(matched),BroadOnlyPatients=';'.join(sorted(broad-clean)),CleanOnlyPatients=';'.join(sorted(clean-broad))))
d=pd.DataFrame(records);pt=pd.DataFrame(pts)
assert len(d)==len(contexts)*17*6
assert not d.duplicated(['cohort','tissue','ProgramID','Other']).any()
qual=d[d.Status.eq('qualified_tumor_n4')]
assert np.isfinite(qual[['BroadEffect','CleanEffect']]).all().all()
assert qual.NMatched.ge(4).all()
# Independently check each aggregate against the saved point records.
for r in d[d.NMatched.gt(0)].itertuples():
 g=pt[(pt.cohort==r.cohort)&(pt.tissue==r.tissue)&(pt.ProgramID==r.ProgramID)&(pt.Other==r.Other)&pt.Matched]
 assert len(g)==r.NMatched and g.PatientKey.nunique()==r.NMatched
 assert abs(g.BroadDifference.mean()-r.BroadEffect)<1e-12 and abs(g.CleanDifference.mean()-r.CleanEffect)<1e-12
sp=pd.DataFrame(dict(ProgramID=['Pattern_P2','WGCNA_MElightgreen','Pattern_P3','WGCNA_MEcyan','Pattern_P4','WGCNA_MEmagenta','WGCNA_MEdarkturquoise'],Other=['Hepatocyte','Hepatocyte','Myeloid','Myeloid','Myeloid','T_NK','Myeloid'],tissue=['Adjacent']*2+['Tumor']*5,QuestionOrder=range(7)))
sel=sp.merge(d,on=['ProgramID','Other','tissue'],how='left')
# Keep only cohort/tissue contexts which actually exist; source full grid is unchanged.
d.to_csv(T/'identity_sensitivity_contrasts.csv',index=False);pt.to_csv(T/'identity_sensitivity_patient_pairs.csv',index=False);sel.to_csv(T/'identity_sensitivity_selected.csv',index=False)
summary=dict(Status='COMPUTED_DESCRIPTIVE_SENSITIVITY',NoRescoring=True,NoNewInference=True,AggregationMaxError=err,FullContrastGrid=len(d),ContextCount=len(contexts),TumorContrastGrid=int(d.tissue.eq('Tumor').sum()),QualifiedTumorContrasts=len(qual),QualifiedByCohort=qual.groupby('cohort').size().to_dict(),Statuses=d.Status.value_counts().to_dict(),SameDirection=int(qual.SameDirection.sum()),DirectionChanges=int((~qual.SameDirection.astype(bool)).sum()),MaxAbsoluteEffectChange=float(qual.EffectChange.abs().max()),MedianAbsoluteEffectChange=float(qual.EffectChange.abs().median()),QualifiedPatientRange=[int(qual.NMatched.min()),int(qual.NMatched.max())],SelectedContrasts=len(sel),SelectedSupported=int(sel.NMatched.gt(0).sum()),QualifiedWithChangedEligibleLibraries=int((qual.NMatchedSameEligibleLibraries<qual.NMatched).sum()),Inputs={f:hashlib.sha256((T/f).read_bytes()).hexdigest() for f in ['patient_program_scores.csv','library_program_scores.csv','cells.csv.gz','shared_program_registry.csv']},Interpretation='Same patients across Broad/Clean; native scores reused. Cell eligibility and eligible-library availability can differ. Descriptive identity-eligibility sensitivity, not independent validation or malignancy confirmation.')
(O/'review/identity_sensitivity_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
print('DIRECTION CHANGES');print(qual[~qual.SameDirection.astype(bool)].to_string(index=False))
