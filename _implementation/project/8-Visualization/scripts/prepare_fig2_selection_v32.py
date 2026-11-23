from pathlib import Path
import pandas as pd,numpy as np,json,hashlib,shutil
R=Path(__file__).resolve().parents[2];O=R/'8-Visualization/refinements/fig2_selection_20261031_v32'
contract=json.loads((O/'SELECTION_CONTRACT.json').read_text());p=R/'1-Clustering/publication/tables'
sources=[p/'step1_ranked_gsea_main_only_study_adjusted.csv',p/'step1_main_only_study_adjusted_deg.csv',R/'2-TranscriptomicCharacterization/runs/fig2_deg_patterns_20261026_v7/reference_gene_sets.csv.gz',R/'2-TranscriptomicCharacterization/publication/tables/step2_wgcna_module_subtype_classification.csv',R/'2-TranscriptomicCharacterization/publication/tables/step2_wgcna_subtype_random_effects_meta.csv',p/'step1_dorothea_tf_activity_main_only_study_adjusted.csv',p/'step1_dorothea_regulon_coverage.csv',p/'step1_dorothea_tf_leave_one_study_out.csv']
hashes={str(x.relative_to(R)):hashlib.md5(x.read_bytes()).hexdigest() for x in sources}
g=pd.read_csv(sources[0]);deg=pd.read_csv(sources[1]);ref=pd.read_csv(sources[2]);universe=set(deg.Gene);assert len(universe)==10531
sets={k:set(v.Gene)&universe for k,v in ref.groupby('Pathway')};assert all(len(sets[r.pathway])==r.size for r in g.itertuples())
g['Direction']=np.where(g.NES>0,'Positive','Negative');eligible=g[(g.padj<.05)&np.isfinite(g.NES)&g['size'].between(10,500)].copy()
strata=[(s,d) for d in ['Positive','Negative'] for s in ['C1','C2','C3']]
queues={(s,d):eligible[(eligible.Subtype==s)&(eligible.Direction==d)].sort_values(['padj','AbsNES','pathway'],ascending=[True,False,True]).reset_index(drop=True) for s,d in strata}
def jac(a,b):return len(sets[a]&sets[b])/len(sets[a]|sets[b])
def select(n,threshold):
 selected=[];ids=set();ptr={s:0 for s in strata};audit=[];roundno=0
 while len(selected)<n:
  roundno+=1;before=len(selected)
  for s,d in strata:
   q=queues[s,d]
   while ptr[s,d]<len(q):
    row=q.iloc[ptr[s,d]];ptr[s,d]+=1;pid=row.pathway
    rec=dict(Pathway=pid,AnchorSubtype=s,AnchorDirection=d,QueueRank=ptr[s,d],Round=roundno,AnchorFDR=row.padj,AnchorNES=row.NES,Collection=row.Collection)
    if pid in ids:rec.update(Decision='Already selected',RedundantWith=pid,Jaccard=1);audit.append(rec);continue
    overlaps=[(jac(pid,k),k) for k in ids];ov,rep=max(overlaps,default=(0,''))
    if ov>=threshold:rec.update(Decision='Membership redundant',RedundantWith=rep,Jaccard=ov);audit.append(rec);continue
    rec.update(Decision='Selected',RedundantWith='',Jaccard=ov,SelectionOrder=len(selected)+1);audit.append(rec);selected.append(rec);ids.add(pid);break
   if len(selected)==n:break
  if len(selected)==before:break
 return pd.DataFrame(selected),pd.DataFrame(audit)
main,audit=select(18,.5);assert len(main)==18
main.to_csv(O/'tables/B_selection18.csv',index=False);audit.to_csv(O/'tables/B_selection_decisions.csv',index=False)
selected=main[['Pathway','AnchorSubtype','AnchorDirection','SelectionOrder']].merge(g,left_on='Pathway',right_on='pathway',validate='one_to_many');assert len(selected)==54
selected.to_csv(O/'tables/B_selected54.csv',index=False)
eligible.to_csv(O/'tables/B_all_significant_candidates.csv',index=False)
summary=[];ss=[]
for threshold in contract['SensitivityB']['JaccardThreshold']:
 for n in contract['SensitivityB']['UniqueSets']:
  s,a=select(n,threshold);s['Threshold']=threshold;s['Budget']=n;ss.append(s)
  ids=set(s.Pathway);base=set(main.Pathway);z=g[g.pathway.isin(ids)]
  summary.append(dict(Threshold=threshold,Budget=n,Selected=len(ids),SharedWithMain18=len(ids&base),JaccardWithMain18=len(ids&base)/len(ids|base),Main18Retained=len(ids&base)/18,GO_BP=(s.Collection=='GO_BP').sum(),Reactome=(s.Collection=='Reactome').sum(),Hallmark=(s.Collection=='Hallmark').sum(),KEGG=(s.Collection=='KEGG').sum(),MaximumPairwiseJaccard=max(jac(a,b) for i,a in enumerate(sorted(ids)) for b in sorted(ids)[i+1:])))
pd.DataFrame(summary).to_csv(O/'tables/B_sensitivity_summary.csv',index=False);sens=pd.concat(ss,ignore_index=True);sens.to_csv(O/'tables/B_sensitivity_membership.csv',index=False)
stability=sens.groupby('Pathway').agg(SelectedInScenarios=('Budget','size')).reset_index();stability['ShownInMain18']=stability.Pathway.isin(main.Pathway);stability.to_csv(O/'tables/B_selection_frequency.csv',index=False)
old=pd.read_csv(R/'8-Visualization/refinements/fig2_redesign_20261031_v30/tables/B_GSEA_selected48.csv');disp=g[['pathway','Collection']].drop_duplicates();disp['PreviouslyShown']=disp.pathway.isin(old.pathway);disp['SignificantAnySubtype']=disp.pathway.isin(eligible.pathway);disp['ShownV32']=disp.pathway.isin(main.Pathway);disp.to_csv(O/'tables/B_all_set_disposition.csv',index=False)
# All eleven, no broad suffix matching that would accidentally include Non_specific.
cls=pd.read_csv(sources[3]);meta=pd.read_csv(sources[4]);c=cls[cls.SpecificityClass.isin(['C1_specific','C2_specific','C3_specific'])].copy();assert len(c)==11
oldc=pd.read_csv(R/'8-Visualization/refinements/fig2_layout_audit_20261031_v31/C/tables/C_all24_effects.csv');labels=oldc.drop_duplicates('Module').set_index('Module').Program.to_dict();labels.update(MEdarkgreen='Cell-cycle checkpoints',MEdarkturquoise='Antiviral / interferon response',MEtan='Muscle cell proliferation')
order=['MElightgreen','MEpink','MEblack','MEcyan','MEbrown','MEdarkgreen','MEmagenta','MEgreenyellow','MEred','MEdarkturquoise','MEtan'];c['DisplayOrder']=c.Module.map({m:i+1 for i,m in enumerate(order)});c['Program']=c.Module.map(labels);c=c.sort_values('DisplayOrder');c['DisplayClass']=c.AssignedSubtypeOrPair;c['PreviousV31']=c.Module.isin(oldc.Module)
c.to_csv(O/'tables/C_all11_modules.csv',index=False);cc=c.merge(meta,on='Module',validate='one_to_many').sort_values(['DisplayOrder','TargetSubtype']);assert len(cc)==33;cc.to_csv(O/'tables/C_all33_effects.csv',index=False)
# TF audit. LOSO is only available for a subset of the full 263 tested TFs.
d=pd.read_csv(sources[5]);cov=pd.read_csv(sources[6]).rename(columns={'source':'TF'});lo=pd.read_csv(sources[7]);assert len(d)==789 and d.TF.nunique()==263
ls=lo.groupby(['TF','Subtype']).agg(LOSOStudies=('HeldOutStudy','nunique'),DirectionRetained=('SignAgreement','sum'),ReducedMinimum=('ReducedDifference','min'),ReducedMaximum=('ReducedDifference','max'),ReducedSignificant=('ReducedFDR',lambda x:(x<.05).sum())).reset_index()
tf=d.merge(cov,on='TF',validate='many_to_one').merge(ls,on=['TF','Subtype'],how='left',validate='one_to_one');tf['CoveragePass']=tf.TargetsAvailable>=10;tf['PositiveSignificant']=(tf.FDR<.05)&(tf.ActivityDifferenceSubtypeMinusOthers>0)&tf.CoveragePass;tf['LOSO7Positive']=tf.LOSOStudies.eq(7)&tf.DirectionRetained.eq(7)&tf.ReducedMinimum.gt(0);tf['PassAllCurrentDisplayChecks']=tf.PositiveSignificant&tf.LOSO7Positive
plan=pd.read_csv(R/'8-Visualization/refinements/fig2_transcriptomic_20261026_v8/tables/tf_selection_plan.csv');tf['CurrentlyShownTF']=tf.TF.isin(plan.TF);ctx=plan.set_index('TF').ContextSubtype.to_dict();tf['IsCurrentAnchor']=tf.apply(lambda r:ctx.get(r.TF)==r.Subtype,axis=1)
rank=tf[tf.PassAllCurrentDisplayChecks].sort_values(['Subtype','FDR','ActivityDifferenceSubtypeMinusOthers','TF'],ascending=[True,True,False,True]).copy();rank['EvidenceRankWithinSubtype']=rank.groupby('Subtype').cumcount()+1
tf=tf.merge(rank[['TF','Subtype','EvidenceRankWithinSubtype']],on=['TF','Subtype'],how='left',validate='one_to_one');tf.to_csv(O/'tables/D_all263_TF_audit.csv',index=False);tf[tf.IsCurrentAnchor].to_csv(O/'tables/D_current9_anchor_audit.csv',index=False)
rank.groupby('Subtype').head(3).to_csv(O/'tables/D_top3_per_subtype_diagnostic_not_replacement.csv',index=False)
tfsum=tf.groupby('Subtype').agg(TestedTFs=('TF','nunique'),PositiveSignificant=('PositiveSignificant','sum'),CompleteLOSO=('LOSOStudies',lambda x:x.eq(7).sum()),PassAllChecks=('PassAllCurrentDisplayChecks','sum'),CurrentAnchors=('IsCurrentAnchor','sum'));tfsum.to_csv(O/'tables/D_selection_summary.csv')
with (O/'review/input_hashes.json').open('w') as f:json.dump(hashes,f,indent=2)
assert all(hashlib.md5((R/k).read_bytes()).hexdigest()==v for k,v in hashes.items())
shutil.copy2(R/'8-Visualization/tables/figure_manifest.csv',O/'review/canonical_manifest_before.csv')
print('B main selection:');print(main[['AnchorSubtype','AnchorDirection','Pathway','AnchorNES']].to_string(index=False));print('\nSensitivity:');print(pd.DataFrame(summary).to_string(index=False));print('\nD audit:');print(tfsum.to_string());print(tf[tf.IsCurrentAnchor][['TF','Subtype','TargetsAvailable','EvidenceRankWithinSubtype']].to_string(index=False))
