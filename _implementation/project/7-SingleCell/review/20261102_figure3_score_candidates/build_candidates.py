from pathlib import Path
import pandas as pd,anndata as ad,json,hashlib
ROOT=Path(__file__).resolve().parents[3];OUT=Path(__file__).resolve().parent
paths={
'pattern':ROOT/'2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/gene_pattern_membership.csv',
'pattern_profile':ROOT/'2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/pattern_state_summary.csv',
'modules':ROOT/'8-Visualization/refinements/fig2_selection_20261031_v32/tables/C_all11_modules.csv',
'membership':ROOT/'2-TranscriptomicCharacterization/publication/tables/step2_wgcna_gene_module_membership.csv'}
p=pd.read_csv(paths['pattern']);mods=pd.read_csv(paths['modules']);w=pd.read_csv(paths['membership']);pr=pd.read_csv(paths['pattern_profile'])
assert len(p)==842 and p.Gene.is_unique and len(mods)==11
rows=[];summary=[]
for prog,g in p.groupby('Pattern',sort=True):
 for r in g.itertuples(index=False):rows.append(dict(ProgramID='Pattern_'+prog,Family='Pattern',Gene=r.Gene,SourceModule=r.WGCNAModule,kME=r.kME,Affinity=r.Affinity,LowAffinity=r.LowAffinity))
 profile=pr[(pr.Pattern==prog)&(pr.SubtypeClass=='Main')].set_index('FinalSubtype')['Mean']
 summary.append(dict(ProgramID='Pattern_'+prog,Family='Pattern',NGenes=len(g),BulkContext=' > '.join(profile.sort_values(ascending=False).index),Description='Bulk expression shape; retain P identifier',LowAffinityGenes=int(g.LowAffinity.sum()),NegativeKMEGenes=None,PreservationZ=None))
for r in mods.sort_values('DisplayOrder').itertuples(index=False):
 g=w[w.Module==r.Module];assert len(g)==r.ModuleSize and g.Gene.is_unique
 for x in g.itertuples(index=False):rows.append(dict(ProgramID='WGCNA_'+r.Module,Family='WGCNA',Gene=x.Gene,SourceModule=x.Module,kME=x.kME,Affinity=None,LowAffinity=None))
 summary.append(dict(ProgramID='WGCNA_'+r.Module,Family='WGCNA',NGenes=len(g),BulkContext=r.DisplayClass,Description=r.Program,LowAffinityGenes=None,NegativeKMEGenes=int((g.kME<0).sum()),PreservationZ=r.ZsummaryPreservation))
member=pd.DataFrame(rows);member.to_csv(OUT/'candidate_gene_membership.csv',index=False);summ=pd.DataFrame(summary);summ.to_csv(OUT/'candidate_programs.csv',index=False)
sets={k:set(g.Gene) for k,g in member.groupby('ProgramID')}
(OUT/'candidate_genesets.gmt').write_text('\n'.join('\t'.join([k,'Figure2-derived_candidate_not_scRNA_validated']+sorted(v)) for k,v in sets.items())+'\n')
ov=[]
for a,sa in sets.items():
 for b,sb in sets.items():
  if a>=b:continue
  ov.append(dict(ProgramA=a,ProgramB=b,NOverlap=len(sa&sb),Jaccard=len(sa&sb)/len(sa|sb)))
pd.DataFrame(ov).to_csv(OUT/'set_overlap.csv',index=False)
ct=pd.crosstab(p.WGCNAModule,p.Pattern).reindex(mods.Module,fill_value=0);ct['PatternDisplayGenes']=ct.sum(axis=1);ct['FullModuleGenes']=mods.set_index('Module').ModuleSize;ct.to_csv(OUT/'pattern_module_membership_counts.csv')
a=ad.read_h5ad(ROOT/'7-SingleCell/runs/20261102_two_cohort_repair_v1/r3_annotation/broad_annotated_atlas.h5ad',backed='r');v=a.varm['native_availability'];obs=a.obs[['cohort','patient','sample','tissue']].drop_duplicates();cov=[]
for r in obs.itertuples(index=False):
 avail=set(v.index[v[r.sample].astype(bool)])
 for prog,genes in sets.items():cov.append(dict(ProgramID=prog,cohort=r.cohort,patient=r.patient,sample=r.sample,tissue=r.tissue,NDefined=len(genes),NAvailable=len(genes&avail),Coverage=len(genes&avail)/len(genes)))
cov=pd.DataFrame(cov);cov.to_csv(OUT/'native_library_coverage.csv',index=False)
cs=cov.groupby(['ProgramID','cohort','tissue']).agg(NLibraries=('sample','nunique'),MinCoverage=('Coverage','min'),MedianCoverage=('Coverage','median'),MaxCoverage=('Coverage','max')).reset_index();cs.to_csv(OUT/'coverage_summary.csv',index=False)
meta={'status':'CANDIDATE_DEFINITIONS_AND_NATIVE_COVERAGE_ONLY','programs':len(sets),'membership_rows':len(member),'unique_genes':member.Gene.nunique(),'library_program_records':len(cov),'pattern_low_affinity':int(p.LowAffinity.sum()),'main_WGCNA_genes':int(mods.ModuleSize.sum()),'scoring_run':False,'canonical_files_changed':False,'source_hashes':{k:{'path':str(v.relative_to(ROOT)),'sha256':hashlib.sha256(v.read_bytes()).hexdigest()} for k,v in paths.items()}}
(OUT/'verification.json').write_text(json.dumps(meta,indent=2));a.file.close();print(summ.to_string(index=False));print(ct.to_string());print(cs[cs.tissue=='Tumor'][['ProgramID','cohort','MinCoverage','MaxCoverage']].to_string(index=False));print(meta)
