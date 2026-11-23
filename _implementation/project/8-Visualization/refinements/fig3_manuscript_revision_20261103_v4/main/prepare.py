from pathlib import Path
import pandas as pd,numpy as np,json,hashlib
O=Path(__file__).resolve().parent;ROOT=O.parents[3];S=ROOT/'8-Visualization/refinements/fig3_scientific_enhancement_20261103_v2/source_data';T=O/'source_data';RUN=ROOT/'7-SingleCell/runs/20261103_figure3_ecosystem_v3/tables'
bridge=ROOT/'8-Visualization/refinements/fig2_global_revision_20261103_v2/pattern_bridge/tables/pattern_module_counts_rowpercent.csv';b=pd.read_csv(bridge)
sp=[('Pattern_P1','Hepatocyte','Adjacent'),('WGCNA_MEpink','Hepatocyte','Adjacent'),('Pattern_P3','Myeloid','Tumor'),('WGCNA_MEbrown','Myeloid','Tumor'),('Pattern_P6','Myeloid','Tumor'),('WGCNA_MEred','Myeloid','Tumor')];prs=[p for p,_,_ in sp]
selected=pd.concat([b[(b.Pattern==p)&(b.Module==m)]for p,m in [('P1','MEpink'),('P3','MEbrown'),('P6','MEred')]])
selected.to_csv(T/'selected_overlap_pairs.csv',index=False);b.to_csv(T/'all_pattern_module_overlap.csv',index=False)
pair=pd.read_csv(S/'source_patient_pairs.csv');stat=pd.read_csv(S/'source_contrasts.csv');pp=[];ss=[]
for j,(pr,other,ti)in enumerate(sp):
 d=pair[(pair.Method=='UCell')&(pair.ProgramID==pr)&(pair.Other==other)&(pair.tissue==ti)].copy();d['PanelOrder']=j+1;pp.append(d)
 d=stat[(stat.Method=='UCell')&(stat.ProgramID==pr)&(stat.Other==other)&(stat.tissue==ti)].copy();d['PanelOrder']=j+1;ss.append(d)
pd.concat(pp).to_csv(T/'D_patient_pairs.csv',index=False);pd.concat(ss).to_csv(T/'D_existing_statistics.csv',index=False);pd.DataFrame(sp,columns=['ProgramID','Other','tissue']).to_csv(T/'D_contrast_specification.csv',index=False)
cells=pd.read_csv(S/'featureplot_cells.csv.gz');z=np.load(RUN/'cell_ucell_scores.npz');assert z['CellKey'].tolist()==cells.CellKey.tolist();cov=pd.read_csv(S/'ucell_coverage.csv').set_index(['sample','ProgramID']);f=cells.iloc[:,:9].copy();masks=[]
for p in prs:
 eligible=f['sample'].map(lambda sm:cov.loc[(sm,p),'Coverage']>=.7 and cov.loc[(sm,p),'NAvailable']>=10)&~f.library_low_support
 f[p]=np.where(eligible,z['UCell'][:,list(z['programs']).index(p)],np.nan).astype(float);masks.append(dict(ProgramID=p,AtlasRows=len(f),EligibleRows=int(eligible.sum()),MaskedRows=int((~eligible).sum())))
f.to_csv(T/'featureplot_selected_cells.csv.gz',index=False);pd.DataFrame(masks).to_csv(T/'featureplot_masks.csv',index=False)
meta=pd.read_csv(S/'cells.csv.gz',usecols=['cohort','patient','PatientKey','tissue']);pats=meta[meta.tissue=='Tumor'][['cohort','patient','PatientKey']].drop_duplicates().sort_values(['cohort','patient'],ascending=[False,True]).reset_index(drop=True);pats['PatientOrder']=range(1,len(pats)+1)
lib=pd.read_csv(S/'library_program_scores.csv');types=['Epithelial','Myeloid','T_NK'];l=lib[(lib.Level=='Clean')&(lib.Method=='UCell')&(lib.tissue=='Tumor')&lib.CellType.isin(types)&lib.ProgramID.isin(prs)].copy();l['Used']=l.Eligible&~l.LowSupport&l.Score.notna();l.to_csv(T/'E_library_inputs.csv',index=False)
u=l[l.Used].groupby(['PatientKey','CellType','ProgramID']).agg(Score=('Score','mean'),NCells=('NCells','sum'),NLibraries=('sample','size')).reset_index();grid=pd.MultiIndex.from_product([pats.PatientKey,types,prs],names=['PatientKey','CellType','ProgramID']).to_frame(index=False).merge(pats,on='PatientKey').merge(u,on=['PatientKey','CellType','ProgramID'],how='left');grid.to_csv(T/'E_patient_compartment_grid.csv',index=False)
for name in ['programs.csv','dotplot_summary.csv','shared_program_registry.csv']:
 import shutil;shutil.copy2(S/name,T/name)
old=pd.read_csv(S/'patient_program_scores.csv');old=old[(old.Level=='Clean')&(old.Method=='UCell')&(old.tissue=='Tumor')&old.CellType.isin(types)&old.ProgramID.isin(prs)];comp=grid.merge(old[['PatientKey','CellType','ProgramID','Score']],on=['PatientKey','CellType','ProgramID'],how='outer',suffixes=('_new','_old'))
checks={'C_all_34227_rows':len(f)==34227,'D_expected_pair_counts':[len(d)for d in pp]==[3,3,9,9,9,9],'E_all_288_grid':len(grid)==288,'E_all16_patients':grid.PatientKey.nunique()==16,'E_patient_means_exact':np.allclose(comp.Score_new,comp.Score_old,equal_nan=True,atol=1e-12),'D_exact_differences':all(np.allclose(d.EpithelialScore-d.OtherScore,d.Difference)for d in pp)}
report=dict(checks=checks,all_pass=all(checks.values()),D_pairs=[len(d)for d in pp],E_estimates=int(grid.Score.notna().sum()),E_missing=int(grid.Score.isna().sum()),feature_masks=masks)
(O/'review/numeric_QA.json').write_text(json.dumps(report,indent=2));assert report['all_pass'],report
print(json.dumps(report,indent=2));print(pd.concat(ss).to_string(index=False))
