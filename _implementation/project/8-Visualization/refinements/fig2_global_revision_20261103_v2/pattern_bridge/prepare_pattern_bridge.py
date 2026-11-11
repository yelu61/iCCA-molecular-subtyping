"""Display registry and descriptive intersections; no inference or rescoring."""
from pathlib import Path
import csv, collections, hashlib, json, shutil

ROOT=Path.cwd()
OUT=ROOT/'8-Visualization/refinements/fig2_global_revision_20261103_v2/pattern_bridge'
S={
 'fig2_members':'8-Visualization/refinements/fig2_layout_audit_20261031_v31/A/tables/row_order.csv',
 'frozen_members':'2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/gene_pattern_membership.csv',
 'pattern_counts':'2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/pattern_counts.csv',
 'state_summary':'2-TranscriptomicCharacterization/runs/fig2_soft_patterns_20261030_v20/K6/pattern_state_summary.csv',
 'pattern_ORA':'8-Visualization/refinements/fig2_layout_audit_20261031_v31/A/tables/annotation_evidence.csv',
 'module_order':'8-Visualization/refinements/sfig2_wgcna_narrow_20261031_v37/tables/module_order.csv',
 'network_members':'2-TranscriptomicCharacterization/publication/tables/step2_wgcna_gene_module_membership.csv',
 'fig3_programs':'7-SingleCell/runs/20261103_figure3_ecosystem_v3/tables/programs.csv',
 'fig3_members':'7-SingleCell/runs/20261103_figure3_ecosystem_v3/tables/membership.csv',
 'fig3_runbook':'7-SingleCell/runs/20261103_figure3_ecosystem_v3/RUNBOOK.md',
 'fig3_contract':'7-SingleCell/runs/20261103_figure3_ecosystem_v3/CONTRACT.md',
 'fig3_manifest':'7-SingleCell/runs/20261103_figure3_ecosystem_v3/manifest.json',
 'fig3_program_builder':'7-SingleCell/review/20261102_figure3_score_candidates/build_candidates.py',
 'subtype_palette':'governance/data/subtype_label_color_dictionary.csv',
 'study_palette':'0-Data/qc/study_color_dictionary.csv',
}
def read(k):
 with open(ROOT/S[k]) as f:return list(csv.DictReader(f))
def write(name,rows):
 with open(OUT/'tables'/name,'w',newline='') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def ghash(genes):return hashlib.sha256(('\n'.join(sorted(genes))+'\n').encode()).hexdigest()
for d in ['figures','tables','review']:(OUT/d).mkdir(parents=True,exist_ok=True)
hashes=[dict(Key=k,Path=p,SHA256=hashlib.sha256((ROOT/p).read_bytes()).hexdigest()) for k,p in S.items()]
hashby={r['Key']:r['SHA256'] for r in hashes}
for k in ['fig2_members','pattern_ORA','fig3_programs','fig3_members']:
 shutil.copyfile(ROOT/S[k],OUT/'tables'/f'source_{k}.csv')
a=read('fig2_members');frozen=read('frozen_members');pc={r['Pattern']:r for r in read('pattern_counts')}
sm={(r['Pattern'],r['FinalSubtype']):r for r in read('state_summary')}
ora=read('pattern_ORA');mods=read('module_order');nw=read('network_members')
programs=read('fig3_programs');mem=read('fig3_members')
assert len(a)==842 and len({r['Gene'] for r in a})==842
assert {(r['Gene'],r['Pattern'],r['WGCNAModule']) for r in a}=={(r['Gene'],r['Pattern'],r['WGCNAModule']) for r in frozen}
assert len(programs)==17 and len(mem)==3321 and len(mods)==24 and len(nw)==6753
labels={
 'P1':('Metabolism / coagulation','Metabolism / coagulation','C1-high profile','代谢 / 凝血'),
 'P2':('Complement / steroid metabolism','Complement / steroid','C1/C2 gradient relative to C3','补体 / 类固醇代谢'),
 'P3':('WNT / ion transport','WNT / ion transport','C2-high profile','WNT / 离子转运'),
 'P4':('Epithelial development','Epithelial development','C2/C3-high relative to C1','上皮发育相关'),
 'P5':('Epithelial / stress response','Epithelial / stress','C3-high profile','上皮 / 应激反应'),
 'P6':('Cytokine / transport','Cytokine / transport','C3 > C1 > C2 descriptive profile','细胞因子 / 转运'),
}
support={
 'P1':{'Retinol metabolism','Coagulation','Xenobiotic metabolism','Monocarboxylic acid metabolism'},
 'P2':{'Complement cascade','Steroid metabolism','Bile acid metabolism'},
 'P3':{'Sodium ion transmembrane transport','Regulation of WNT signaling'},
 'P4':{'Embryonic organ development','Epithelial cell differentiation'},
 'P5':{'Epidermis development','TNF / NF-kB signaling','P53 pathway','mTORC1 signaling'},
 'P6':{'Monocarboxylic acid transport','Interferon-alpha response','Cytokine-mediated signaling'},
}
reg=[];checks=[]
for r in programs:
 pid=r['ProgramID'];genes={x['Gene'] for x in mem if x['ProgramID']==pid}
 if r['Family']=='Pattern':
  p=pid.removeprefix('Pattern_');expected={x['Gene'] for x in a if x['Pattern']==p}
  full,short,context,zh=labels[p];count=pc[p];evidence=[x for x in ora if x['Pattern']==p and x['DisplayLabel'] in support[p]]
  assert len(evidence)==len(support[p])
  assert all(x['SourceModule']==next(z['WGCNAModule'] for z in a if z['Gene']==x['Gene']) for x in mem if x['ProgramID']==pid)
  desc='Post hoc bulk-ORA description; not exclusive, cell identity, mechanism or naming lock'
  row=dict(RegistryVersion='20261103_v2',ProgramID=pid,Family=r['Family'],StableID=p,DisplayDescription=full,CompactDescription=short,DisplayLabel=f'{p} | {full}',CompactLabel=f'{p} | {short}',ChineseDescription=zh,NGenes=len(genes),MembersSHA256=ghash(genes),OriginalBulkContext=r['BulkContext'],SharedDescriptiveContext=context,OriginalDescription=r['Description'],LowAffinityGenes=int(float(r['LowAffinityGenes'])),NamingStatus=desc,NameEvidencePathwayIDs=';'.join(x['Pathway'] for x in evidence),NameEvidenceFDRs=';'.join(x['FDR'] for x in evidence),Figure2MembershipPath=S['fig2_members'],Figure2MembershipSourceSHA256=hashby['fig2_members'],Figure3MembershipPath=S['fig3_members'],Figure3MembershipSourceSHA256=hashby['fig3_members'],Figure3ProgramPath=S['fig3_programs'],Figure3ProgramSourceSHA256=hashby['fig3_programs'],ORAPath=S['pattern_ORA'],ORASourceSHA256=hashby['pattern_ORA'],ShapeC1=count['ShapeC1'],ShapeC2=count['ShapeC2'],ShapeC3=count['ShapeC3'],BulkMeanC1=sm[p,'C1']['Mean'],BulkMeanC2=sm[p,'C2']['Mean'],BulkMeanC3=sm[p,'C3']['Mean'])
 else:
  module=pid.removeprefix('WGCNA_');expected={x['Gene'] for x in nw if x['Module']==module};label=r['Description']
  row=dict(RegistryVersion='20261103_v2',ProgramID=pid,Family=r['Family'],StableID=module,DisplayDescription=label,CompactDescription=label,DisplayLabel=f'{module} | {label}',CompactLabel=f'{module} | {label}',ChineseDescription='',NGenes=len(genes),MembersSHA256=ghash(genes),OriginalBulkContext=r['BulkContext'],SharedDescriptiveContext=r['BulkContext'],OriginalDescription=r['Description'],LowAffinityGenes='',NamingStatus='Existing full-module bulk annotation retained; not cell identity or causal regulation',NameEvidencePathwayIDs='',NameEvidenceFDRs='',Figure2MembershipPath=S['network_members'],Figure2MembershipSourceSHA256=hashby['network_members'],Figure3MembershipPath=S['fig3_members'],Figure3MembershipSourceSHA256=hashby['fig3_members'],Figure3ProgramPath=S['fig3_programs'],Figure3ProgramSourceSHA256=hashby['fig3_programs'],ORAPath='',ORASourceSHA256='',ShapeC1='',ShapeC2='',ShapeC3='',BulkMeanC1='',BulkMeanC2='',BulkMeanC3='')
 assert genes==expected and len(genes)==int(r['NGenes']),pid
 checks.append(dict(ProgramID=pid,NGenes=len(genes),ExactFigure2Figure3Members=True,MissingFromFigure3=0,ExtraInFigure3=0,MembersSHA256=ghash(genes)))
 reg.append(row)
write('shared_program_registry.csv',reg);write('program_membership_equivalence.csv',checks)
write('source_hashes.csv',hashes)
for x in ora:x['SupportsWorkingName']=x['DisplayLabel'] in support[x['Pattern']]
write('pattern_naming_ORA_evidence.csv',ora)

cols=[]
for i,m in enumerate(mods):
 c=dict(Module=m['Module'],ModuleLabel=m['ModuleColor'],ModuleColor=m['ModuleColor'],ColumnOrder=i+1,FullModuleN=int(m['ModuleSize']),InMainFigure2C=m['MainFigure'],InCurrentFigure3=any(p['ProgramID']=='WGCNA_'+m['Module'] for p in programs),MembershipCategory='Named network module')
 assert len([x for x in nw if x['Module']==m['Module']])==c['FullModuleN']
 cols.append(c)
cols += [dict(Module='grey',ModuleLabel='grey (unassigned)',ModuleColor='#B8B8B8',ColumnOrder=25,FullModuleN='',InMainFigure2C='False',InCurrentFigure3=False,MembershipCategory='Unassigned grey'),dict(Module='Not_in_network',ModuleLabel='Outside network',ModuleColor='#EEEEEE',ColumnOrder=26,FullModuleN='',InMainFigure2C='False',InCurrentFigure3=False,MembershipCategory='Outside network')]
counts=collections.Counter((x['Pattern'],x['WGCNAModule']) for x in a)
matrix=[];gene_rows=[]
for p in labels:
 n=sum(x['Pattern']==p for x in a)
 for c in cols:
  ct=counts[p,c['Module']]
  matrix.append(dict(Pattern=p,PatternN=n,PatternDescription=labels[p][0],Module=c['Module'],ModuleLabel=c['ModuleLabel'],ColumnOrder=c['ColumnOrder'],IntersectionGenes=ct,RowPercent=100*ct/n,FullModuleN=c['FullModuleN'],PercentOfFullModule=(100*ct/c['FullModuleN'] if c['FullModuleN'] else ''),InCurrentFigure3=c['InCurrentFigure3']))
 assert sum(x['IntersectionGenes'] for x in matrix if x['Pattern']==p)==n
 for x in a:
  if x['Pattern']==p:gene_rows.append(dict(ProgramID='Pattern_'+p,Pattern=p,Gene=x['Gene'],Module=x['WGCNAModule'],LowAffinity=x['LowAffinity'],Affinity=x['Affinity'],SourceRowOrder=x['RowOrder']))
for c in cols:c['GenesIn842']=sum(r['IntersectionGenes'] for r in matrix if r['Module']==c['Module'])
write('pattern_module_counts_rowpercent.csv',matrix);write('module_column_registry.csv',cols);write('pattern_gene_module_membership.csv',gene_rows)
assert len(matrix)==156 and sum(x['IntersectionGenes'] for x in matrix)==842
assert sum(x['IntersectionGenes'] for x in matrix if x['Module']=='grey')==42
assert sum(x['IntersectionGenes'] for x in matrix if x['Module']=='Not_in_network')==6
for r in hashes:assert hashlib.sha256((ROOT/r['Path']).read_bytes()).hexdigest()==r['SHA256']
summary=dict(Status='Descriptive registry and intersection draft; no upstream modification',CurrentFigure3Run='20261103_figure3_ecosystem_v3',PatternSetsExact=6,CompleteModuleSetsExact=11,Programs=17,ProgramMemberRows=len(mem),PatternGenes=842,NamedModuleGenesIn842=794,UnassignedGreyIn842=42,OutsideNetworkIn842=6,MatrixRows=6,MatrixColumns=26,MatrixCells=156,ZeroNamedModuleColumns=[c['Module'] for c in cols if c['MembershipCategory']=='Named network module' and c['GenesIn842']==0],NamedFullNetworkGenes=6753,SourceHashesUnchanged=True,NoEnrichmentTestsPerformed=True,P4ContextNote='Current Figure3 BulkContext reflects arithmetic bulk means (C3 slightly > C2); fuzzy adjusted profile coordinates put C2 slightly > C3. Both summaries are retained; shared display context is C2/C3-high relative to C1, not a rank-equivalence claim.')
(OUT/'review/data_verification.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
