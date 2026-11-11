from pathlib import Path
import csv,hashlib,json
root=Path.cwd();out=root/'8-Visualization/refinements/fig6_evidence_20261103_v21/qa/candidate_audit'
paths=[]
def rd(p):
 paths.append(p)
 with (root/p).open() as f:return list(csv.DictReader(f))
def wr(name,a):
 keys=list(dict.fromkeys(k for r in a for k in r))
 with (out/name).open('w') as f:w=csv.DictWriter(f,keys);w.writeheader();w.writerows(a)
b='5-ClassifierValidation/review/20261101_drug_evidence_refresh_v1/tables/'
ev=rd(b+'candidate_target_multiomic_evidence.csv');ret=rd(b+'original_candidate_retention.csv');pr=rd(b+'PRISM_candidate_model_cells.csv');ann=rd(b+'drug_target_resource_annotations.csv');meta=rd('8-Visualization/refinements/fig6_redesign_20261101_v18/tables/Figure6G_candidate_evidence.csv')
screen=rd('8-Visualization/refinements/fig6_application_suite_20261031_v7/tables/all_subtype_candidate_screen.csv')
gs=rd('1-Clustering/publication/tables/step1_ranked_gsea_main_only_study_adjusted.csv');sel=rd('8-Visualization/refinements/fig2_supp_reorganization_20261103_v3/expression/source_data/GSEA_main_selection_B_selection18.csv');selected={r['Pathway'] for r in sel}
# Retrospective descriptive crosswalk based on annotated mechanisms, not ranked NES/FDR.
spec=[
('AZD1208_1449','C1','PIM kinases','PIM1;PIM2;PIM3','','Unmapped','No uniquely named PIM-specific tested program; do not assign broad survival/proliferation as target activity.'),
('JAK1_8709_1718','C1','JAK1 kinase','JAK1','GOBP_CELL_SURFACE_RECEPTOR_SIGNALING_PATHWAY_VIA_JAK_STAT','Family context','No generic Reactome JAK/STAT entry in tested table; GO family-level signaling is preferred to choosing a ligand-specific Reactome result.'),
('LY2109761_1852','C1','TGF-beta receptors','TGFBR1;TGFBR2','REACTOME_SIGNALING_BY_TGF_BETA_RECEPTOR_COMPLEX','Receptor context','Receptor-complex identity; correct erroneous GDSC TGFB1 ligand annotation using original study PMID16452224.'),
('Ibrutinib_1799','C2','BTK kinase','BTK','REACTOME_SIGNALING_BY_THE_B_CELL_RECEPTOR_BCR','Family context','BCR pathway is broad BTK-related context, not malignant-cell BTK activity or tumor dependence.'),
('PRIMA-1MET_1131','C2','p53-associated','TP53','REACTOME_TRANSCRIPTIONAL_REGULATION_BY_TP53','Regulatory context','Mechanism-associated TP53; RNA amount does not establish mutation state, reactivation potential or compound dependence.'),
('Sepantronium bromide_1941','C2','Survivin-associated','BIRC5','REACTOME_CELL_CYCLE_MITOTIC','Broad mitotic context','Broad BIRC5-associated mitotic background only; not an apoptosis or BIRC5-specific activity score. BIRC5 belongs to C2/C3 leading edges.'),
('AZD7762_1022','C3','CHK1/2 kinases','CHEK1;CHEK2','REACTOME_CELL_CYCLE_CHECKPOINTS','Checkpoint context','Named cell-cycle checkpoint program; does not establish checkpoint addiction or inhibitor efficacy.'),
('PD0325901_1060','C3','MEK1/2 kinases','MAP2K1;MAP2K2','REACTOME_MAP2K_AND_MAPK_ACTIVATION','Cascade context','Nearest named MAP2K/MAPK cascade entry; enrichment is transcriptional context, not phosphorylation/activity.'),
('Ulixertinib_1908','C3','ERK1/2 kinases','MAPK1;MAPK3','REACTOME_MAPK1_MAPK3_SIGNALING','Cascade context','Exact named MAPK1/MAPK3 signaling program; retain nonsignificant GSEA despite significant MAPK3 abundance.')]
maprows=[];pg=[];summary=[]
for d,st,cl,genes,pn,scope,reason in spec:
 a=next(r for r in ann if r['Drug']==d)
 mr=dict(Drug=d,CandidateSubtype=st,FunctionalClass=cl,MappedGenes=genes,ProgramID=pn,ProgramScope=scope,MappingReason=reason,SourceTarget=a['SourceTarget'],SourcePathway=a['SourcePathway'],ProgramInMainFigure2_18=pn in selected,MechanismSource='GDSC release8.2; candidate_mechanism_gene_map.csv',FunctionalCaveat='TP53/BIRC5 are mechanism-associated; do not assert direct binding from this dataset' if 'associated' in cl else 'Target annotation does not establish selectivity or efficacy')
 maprows.append(mr)
 pp=[r for r in gs if r['pathway']==pn] if pn else []
 assert len(pp)==3 if pn else len(pp)==0
 for g in pp:pg.append(dict(Drug=d,CandidateSubtype=st,ProgramScope=scope,**g))
 target=next((r for r in pp if r['Subtype']==st),{})
 sc=next(r for r in screen if r['Drug']==d and r['Subtype']==st)
 tt=[r for r in ev if r['Drug']==d and r['Subtype']==st]
 def direction(assay):
  groups={'Up':[],'Down':[],'NS':[],'NA':[]}
  for r in tt:
   x,f=r[assay+'Effect'],r[assay+'FDR'];k='NA' if not x or not f else 'Up' if float(f)<.05 and float(x)>0 else 'Down' if float(f)<.05 and float(x)<0 else 'NS';groups[k].append(r['Gene'])
  return '; '.join(k+':'+','.join(v) for k,v in groups.items() if v)
 mm=[r for r in meta if r['Drug']==d];rr=next(r for r in ret if r['Drug']==d);zz=[r for r in pr if r['Drug']==d];meas=[r for r in zz if r['ObservationStatus']=='Measured']
 summary.append(dict(**mr,OriginalSelected=sc.get('Selected',''),OriginalEligible=sc.get('Eligible',''),TrainingPearson=sc.get('Pearson_r',''),BatchStandardizeSpearman=sc.get('Primary_vs_pooled_standardize_Spearman',''),BatchEBSpearman=sc.get('Primary_vs_pooled_eb_Spearman',''),StudyLowerFraction=sc.get('FractionLower',''),ProgramNES=target.get('NES',''),ProgramFDR=target.get('padj',''),ProgramState='Unmapped' if not target else 'Up' if float(target['padj'])<.05 and float(target['NES'])>0 else 'Down' if float(target['padj'])<.05 and float(target['NES'])<0 else 'NS',RNA=direction('RNA'),Protein=direction('Protein'),PairwiseCIBelowZero=sum(float(r['Upper'])<0 for r in mm),PairwiseContrasts=len(mm),MinI2=min(float(r['I2_percent']) for r in mm),MaxI2=max(float(r['I2_percent']) for r in mm),Top3Retained=rr['TopThreeRetained'],Omissions=rr['Omissions'],PRISMMeasured=len(meas),PRISMGrid=len(zz),PRISMMeasuredC1=sum(r['Prediction']=='C1' for r in meas),PRISMMeasuredC2=sum(r['Prediction']=='C2' for r in meas),PRISMMeasuredC3=sum(r['Prediction']=='C3' for r in meas),PRISMNoAcceptedIdentity=sum(r['ObservationStatus']=='No accepted assay identity' for r in zz),PRISMMatchedMissing=sum(r['ObservationStatus']=='Matched assay; missing measurement' for r in zz),PRISM24OriginalCandidateMeasurements=0))
wr('candidate_program_map.csv',maprows);wr('candidate_program_all_subtypes.csv',pg);wr('candidate_evidence_summary.csv',summary)
ids={r['Drug'] for r in maprows};wr('candidate_PRISM_all_cells.csv',[r for r in pr if r['Drug'] in ids]);wr('candidate_target_all_subtypes.csv',[r for r in ev if r['Drug'] in ids]);wr('candidate_pairwise_meta.csv',[r for r in meta if r['Drug'] in ids])
wr('input_hashes.csv',[dict(Path=p,SHA256=hashlib.sha256((root/p).read_bytes()).hexdigest()) for p in paths])
print(json.dumps(dict(drugs=len(maprows),mapped_programs=sum(bool(r['ProgramID']) for r in maprows),program_rows=len(pg),main18_overlap=sum(r['ProgramInMainFigure2_18'] for r in maprows),prism_measured=sum(r['PRISMMeasured'] for r in summary)),ensure_ascii=False))
