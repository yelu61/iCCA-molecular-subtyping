from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil
from scipy.stats import false_discovery_control
R=Path(__file__).resolve().parents[2];O=R/'8-Visualization/refinements/fig2D_selection_20261031_v33';P=R/'1-Clustering/publication/tables';V=R/'8-Visualization/refinements/fig2_transcriptomic_20261026_v9/tables'
sources=[P/'step1_dorothea_tf_activity_main_only_study_adjusted.csv',P/'step1_dorothea_regulon_coverage.csv',P/'step1_dorothea_tf_leave_one_study_out.csv',V/'tf_reconstructed_measured_regulon.csv',V/'tf_regulon_annotation_provenance.csv',P/'step1_dorothea_tf_method_contract.csv']
hashes={str(p.relative_to(R)):hashlib.md5(p.read_bytes()).hexdigest() for p in sources}
d=pd.read_csv(sources[0]);cv=pd.read_csv(sources[1]).rename(columns={'source':'TF'});lo=pd.read_csv(sources[2]);net=pd.read_csv(sources[3]);assert len(d)==789 and d.TF.nunique()==263 and d.NTotal.eq(497).all() and d.StudyAdjusted.all()
assert not net.duplicated(['TF','Gene']).any() and np.isfinite(net.mor).all()
counts=net.groupby('TF').Gene.nunique();assert len(counts)==271 and cv.set_index('TF').TargetsAvailable.eq(counts).all()
net['Sign']=np.sign(net.mor).astype(int);nz=net[net.Sign!=0];signed={tf:set(zip(z.Gene,z.Sign)) for tf,z in nz.groupby('TF')};unsigned={tf:set(z.Gene) for tf,z in nz.groupby('TF')};signedcounts=nz.groupby('TF').Gene.nunique().rename('NonzeroSignedTargets')
# The static method-contract sentence says combined contrasts; original numbers actually use per-subtype BH.
fdr=[]
for st,z in d.groupby('Subtype'):
 q=false_discovery_control(z.PValue.to_numpy(),method='bh');fdr.append(dict(Subtype=st,N=len(z),MaximumAbsoluteDifference=float(np.max(np.abs(q-z.FDR.to_numpy())))));assert np.max(np.abs(q-z.FDR.to_numpy()))<1e-12
pd.DataFrame(fdr).to_csv(O/'tables/FDR_family_verification.csv',index=False)
assert not lo.duplicated(['TF','Subtype','HeldOutStudy']).any(); assert set(lo.HeldOutStudy)=={'2019Ahn','2022Dong','2023Deng','2023Jeon','2023Lee','2024Fan','TCGA-CHOL'}
merged=lo.merge(d[['TF','Subtype','ActivityDifferenceSubtypeMinusOthers','FDR']],on=['TF','Subtype'],validate='many_to_one');assert np.allclose(merged.FullDifference,merged.ActivityDifferenceSubtypeMinusOthers) and np.allclose(merged.FullFDR,merged.FDR)
assert ((np.sign(merged.FullDifference)==np.sign(merged.ReducedDifference))==merged.SignAgreement).all()
ls=lo.groupby(['TF','Subtype']).agg(LOSOStudies=('HeldOutStudy','nunique'),DirectionRetained=('SignAgreement','sum'),LOSOmin=('ReducedDifference','min'),LOSOmax=('ReducedDifference','max'),LOSOSignificant=('ReducedFDR',lambda v:(v<.05).sum())).reset_index()
a=d.merge(cv,on='TF',validate='many_to_one').merge(signedcounts,on='TF',validate='many_to_one').merge(ls,on=['TF','Subtype'],how='left',validate='one_to_one');a['Direction']=np.where(a.ActivityDifferenceSubtypeMinusOthers>0,'Positive','Negative');a['AbsEffect']=a.ActivityDifferenceSubtypeMinusOthers.abs();a['CoveragePass']=(a.TargetsAvailable>=10)&(a.NonzeroSignedTargets>=10);a['OriginalSignificant']=a.FDR<.05;a['LOSOStable']=a.LOSOStudies.eq(7)&a.DirectionRetained.eq(7);a['Eligible']=a.CoveragePass&a.OriginalSignificant&a.LOSOStable&np.isfinite(a.AbsEffect)&a.AbsEffect.gt(0)
a['ExclusionReason']=a.apply(lambda r:'Eligible' if r.Eligible else '; '.join([s for ok,s in [(r.CoveragePass,'Insufficient target coverage'),(r.OriginalSignificant,'Original FDR >= 0.05'),(r.LOSOStudies==7,'Incomplete existing LOSO'),(r.DirectionRetained==7,'LOSO sign not retained in all seven')] if not ok]),axis=1)
a.to_csv(O/'tables/D_all263_eligibility.csv',index=False);a.groupby(['Subtype','Direction']).agg(Tested=('TF','size'),Eligible=('Eligible','sum')).to_csv(O/'tables/D_candidate_counts.csv')
strata=[(s,di) for di in ['Positive','Negative'] for s in ['C1','C2','C3']];queues={(s,di):a[a.Eligible&a.Subtype.eq(s)&a.Direction.eq(di)].sort_values(['FDR','AbsEffect','TF'],ascending=[True,False,True]).reset_index(drop=True) for s,di in strata}
ranked=pd.concat([z.assign(QueueRank=np.arange(1,len(z)+1)) for z in queues.values()],ignore_index=True);ranked.to_csv(O/'tables/D_all_eligible_queue_ranks.csv',index=False)
def jac(x,y,mode):
 sets=signed if mode=='signed' else unsigned
 return len(sets[x]&sets[y])/len(sets[x]|sets[y])
def select(n,cut,mode='signed'):
 selected=[];ids=set();ptr={s:0 for s in strata};audit=[];roundno=0
 while len(selected)<n:
  roundno+=1;before=len(ids)
  for s,di in strata:
   q=queues[s,di]
   while ptr[s,di]<len(q):
    row=q.iloc[ptr[s,di]];ptr[s,di]+=1;tf=row.TF
    rec=dict(TF=tf,AnchorSubtype=s,AnchorDirection=di,QueueRank=ptr[s,di],Round=roundno,AnchorFDR=row.FDR,AnchorEffect=row.ActivityDifferenceSubtypeMinusOthers,TargetsAvailable=row.TargetsAvailable,NonzeroSignedTargets=row.NonzeroSignedTargets,LOSOStudies=row.LOSOStudies,DirectionRetained=row.DirectionRetained)
    if tf in ids:rec.update(Decision='Already selected',RedundantWith=tf,Jaccard=1);audit.append(rec);continue
    ov,rep=max([(jac(tf,x,mode),x) for x in ids],default=(0,''))
    if ov>=cut:rec.update(Decision='Redundant measured regulon',RedundantWith=rep,Jaccard=ov);audit.append(rec);continue
    rec.update(Decision='Selected',RedundantWith='',Jaccard=ov,SelectionOrder=len(selected)+1);selected.append(rec);ids.add(tf);audit.append(rec);break
   if len(selected)==n:break
  if len(ids)==before:break
 return pd.DataFrame(selected),pd.DataFrame(audit)
main,trace=select(12,.5);assert len(main)==12;main.to_csv(O/'tables/D_selection12.csv',index=False);trace.to_csv(O/'tables/D_selection_decisions.csv',index=False)
sel=main[['TF','AnchorSubtype','AnchorDirection','SelectionOrder']].merge(a,on='TF',validate='one_to_many');assert len(sel)==36;sel.to_csv(O/'tables/D_selected36.csv',index=False)
rows=[];summaries=[];base=set(main.TF)
for mode,cut,n in [('signed',c,n) for c in [.35,.5,.65] for n in [6,12,18]]+[('unsigned',.5,12)]:
 z,t=select(n,cut,mode);z['Mode']=mode;z['Threshold']=cut;z['Budget']=n;rows.append(z);ids=set(z.TF);summaries.append(dict(Mode=mode,Threshold=cut,Budget=n,Selected=len(z),SharedWithMain12=len(ids&base),JaccardWithMain12=len(ids&base)/len(ids|base),MaximumPairwiseOverlap=max(jac(x,y,mode) for i,x in enumerate(sorted(ids)) for y in sorted(ids)[i+1:])))
sens=pd.concat(rows,ignore_index=True);sens.to_csv(O/'tables/D_sensitivity_membership.csv',index=False);pd.DataFrame(summaries).to_csv(O/'tables/D_sensitivity_summary.csv',index=False)
p18=sens[(sens.Mode=='signed')&(sens.Threshold==.5)&(sens.Budget==18)][['TF','AnchorSubtype','AnchorDirection','SelectionOrder']];ext=p18.merge(a,on='TF',validate='one_to_many');ext.to_csv(O/'tables/SFigure2D_extended54.csv',index=False)
old=pd.read_csv(R/'8-Visualization/refinements/fig2_transcriptomic_20261026_v8/tables/tf_selection_plan.csv');dis=a[['TF']].drop_duplicates();dis['Old9']=dis.TF.isin(old.TF);dis['Main12']=dis.TF.isin(base);dis['EligibleAny']=dis.TF.isin(a[a.Eligible].TF);dis['Extended18']=dis.TF.isin(ext.TF);dis.to_csv(O/'tables/D_old_new_disposition.csv',index=False)
pairs=[]
for i,x in enumerate(sorted(base)):
 for y in sorted(base)[i+1:]:pairs.append(dict(TF1=x,TF2=y,SignedJaccard=jac(x,y,'signed'),UnsignedJaccard=jac(x,y,'unsigned')))
pd.DataFrame(pairs).to_csv(O/'tables/D_selected_regulon_overlap.csv',index=False)
net[net.TF.isin(sens.TF)].to_csv(O/'tables/D_display_union_regulon_edges.csv',index=False)
(O/'review/input_hashes.json').write_text(json.dumps(hashes,indent=2));shutil.copy2(R/'8-Visualization/tables/figure_manifest.csv',O/'review/canonical_manifest_before.csv');assert all(hashlib.md5((R/p).read_bytes()).hexdigest()==v for p,v in hashes.items())
print(main[['TF','AnchorSubtype','AnchorDirection','QueueRank','AnchorEffect','TargetsAvailable']].to_string(index=False));print(pd.DataFrame(summaries).to_string(index=False));print(a.groupby(['Subtype','Direction']).Eligible.sum());print('zero edges',int(net.Sign.eq(0).sum()),'old retained',sorted(set(old.TF)&base),'LOSO incomplete',a[a.LOSOStudies.isna()].TF.nunique())
