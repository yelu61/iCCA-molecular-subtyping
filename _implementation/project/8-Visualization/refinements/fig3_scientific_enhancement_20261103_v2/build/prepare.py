from pathlib import Path
import pandas as pd, numpy as np, hashlib, json, shutil, itertools

ROOT=Path.cwd(); O=ROOT/'8-Visualization/refinements/fig3_scientific_enhancement_20261103_v2'
T=ROOT/'7-SingleCell/runs/20261103_figure3_ecosystem_v3/tables'
for d in ['source_data','figures','review','build']:(O/d).mkdir(exist_ok=True)
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
src=ROOT/'manuscript/figure_drafts/iCCC_subtypes_figures-V2.pptx'
shutil.copy2(src,O/'build/source_snapshot.pptx')
reg=pd.read_csv(ROOT/'8-Visualization/refinements/fig2_global_revision_20261103_v2/pattern_bridge/tables/shared_program_registry.csv')
mem=pd.read_csv(T/'membership.csv'); checks=[]
gene_col=next(c for c in mem if c.lower() in ['gene','symbol','genesymbol'])
for r in reg.itertuples():
    genes=sorted(mem.loc[mem.ProgramID.eq(r.ProgramID),gene_col].tolist())
    h=hashlib.sha256(('\n'.join(genes)+'\n').encode()).hexdigest()
    assert len(genes)==r.NGenes and h==r.MembersSHA256
    checks.append(dict(ProgramID=r.ProgramID,NGenes=len(genes),MembersSHA256=h))
assert len(reg)==17 and len(mem)==3321
reg.to_csv(O/'source_data/shared_program_registry.csv',index=False)
pd.DataFrame(checks).to_csv(O/'review/membership_verification.csv',index=False)
names=['programs.csv','membership.csv','patient_program_scores.csv','library_program_scores.csv','dotplot_summary.csv','cohort_summary.csv','source_contrasts.csv','source_patient_pairs.csv','ucell_coverage.csv','patient_associations.csv','patient_association_points.csv','within_epithelial_program_correlations.csv','method_concordance.csv','coavailability.csv','cells.csv.gz','featureplot_cells.csv.gz','association_specification.csv']
hashes=[]
for n in names:
    shutil.copy2(T/n,O/'source_data'/n);hashes.append(dict(Path=str((T/n).relative_to(ROOT)),SHA256=sha(T/n)))
pd.DataFrame(hashes).to_csv(O/'review/input_hashes.csv',index=False)
cells=pd.read_csv(T/'cells.csv.gz');pat=pd.read_csv(T/'patient_program_scores.csv');lib=pd.read_csv(T/'library_program_scores.csv')
patients=cells[cells.tissue.eq('Tumor')][['cohort','patient','PatientKey']].drop_duplicates().sort_values(['cohort','patient'],ascending=[False,True])
ep=pat.query("Level=='Clean' and Method=='UCell' and tissue=='Tumor' and CellType=='Epithelial'")
grid=pd.DataFrame(itertools.product(patients.PatientKey,reg.ProgramID),columns=['PatientKey','ProgramID']).merge(patients,on='PatientKey').merge(ep.drop(columns=['patient','cohort']),on=['PatientKey','ProgramID'],how='left')
assert len(grid)==16*17 and ep.PatientKey.nunique()==12
grid.to_csv(O/'source_data/main_epithelial_patient_grid.csv',index=False)
# Verify aggregation against the already executed library table, without re-scoring.
eligible=lib.query("Level=='Clean' and Method=='UCell' and tissue=='Tumor' and CellType=='Epithelial' and Eligible")
calc=eligible.groupby(['PatientKey','ProgramID']).Score.mean().rename('Recomputed').reset_index()
cmp=ep.merge(calc,on=['PatientKey','ProgramID']);assert np.allclose(cmp.Score,cmp.Recomputed,rtol=0,atol=1e-12)
sp=pd.DataFrame(dict(ProgramID=['Pattern_P2','WGCNA_MElightgreen','Pattern_P3','WGCNA_MEcyan','Pattern_P4','WGCNA_MEmagenta','WGCNA_MEdarkturquoise'],Other=['Hepatocyte','Hepatocyte','Myeloid','Myeloid','Myeloid','T_NK','Myeloid'],tissue=['Adjacent']*2+['Tumor']*5,Row=range(7)))
c=pd.read_csv(T/'source_contrasts.csv');pairs=pd.read_csv(T/'source_patient_pairs.csv')
sel=sp.merge(c[c.Method.eq('UCell')],on=['ProgramID','Other','tissue']);sel.to_csv(O/'source_data/main_source_contrasts.csv',index=False)
for r in sel.itertuples():
    x=pairs[(pairs.Method=='UCell')&(pairs.ProgramID==r.ProgramID)&(pairs.Other==r.Other)&(pairs.tissue==r.tissue)&(pairs.cohort==r.cohort)]
    assert len(x)==r.NPatients and abs(x.Difference.mean()-r.Effect)<1e-12
    if r.NPatients>=4:
        d=x.Difference.to_numpy();signs=np.array(list(itertools.product([-1,1],repeat=len(d))))
        pp=np.mean(np.abs(signs@d/len(d))>=abs(d.mean())-1e-12)
        assert abs(pp-r.P)<1e-12
assert c.query("Method=='UCell' and tissue=='Tumor'").FDR.min()>=.05
assoc=pd.read_csv(T/'patient_associations.csv').query("Method=='UCell'");assert len(assoc)==6 and assoc.FDR.min()>=.05
# Store the full contrast background, including NE and non-significant results.
sc=c.query("Method=='UCell'");assert len(sc)==272
supp=lib.query("Level=='Broad' and Method=='UCell' and not LowSupport").drop_duplicates(['PatientKey','sample','tissue','CellType']).groupby(['cohort','patient','PatientKey','tissue','CellType']).agg(NCells=('NCells','sum')).reset_index()
# Count unique supported programs across libraries rather than summing library-program rows.
nprog=lib.query("Level=='Broad' and Method=='UCell' and Eligible").groupby(['PatientKey','tissue','CellType']).ProgramID.nunique().rename('NPrograms')
supp=supp.merge(nprog,on=['PatientKey','tissue','CellType'],how='left').fillna({'NPrograms':0})
supp.to_csv(O/'source_data/patient_support_summary.csv',index=False)
summary=dict(Status='REVIEW_NOT_FINALIZED',Source=str(src),SourceSHA256=sha(src),SourceSlideCount=19,FullAtlasCells=len(cells),AllPatients=cells.PatientKey.nunique(),TumorPatients=len(patients),Libraries=cells['sample'].nunique(),TumorLibraries=cells[cells.tissue.eq('Tumor')]['sample'].nunique(),CleanEpithelialPatients=ep.PatientKey.nunique(),PatientProgramGrid=len(grid),PatientProgramSupported=int(grid.Score.notna().sum()),CompletePrograms=len(reg),MembershipRecords=len(mem),MainSourceContrasts=len(sel),CompleteSourceContrasts=len(sc),MinTumorSourceFDR=float(sc[sc.tissue.eq('Tumor')].FDR.min()),MinTMEFDR=float(assoc.FDR.min()),NoRescoring=True)
(O/'review/numeric_verification.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
