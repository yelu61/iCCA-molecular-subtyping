from pathlib import Path
import pandas as pd
import numpy as np
import json,hashlib
from scipy.spatial import cKDTree
R=Path(__file__).resolve().parent;P=R.parent/'20261103_figure3_opening_umap_v2'
s=pd.read_csv(P/'source_data/integrated_umap_cells.csv.gz',keep_default_na=False)
s['PlotCellType']=s.CellTypeDisplay.replace({'Monocyte-like':'Mono/macrophage','Macrophage':'Mono/macrophage','Myeloid, other':'Mono/macrophage'})
s['PlotExclusion']=''
s.loc[s.CellTypeDisplay.eq('Unresolved'),'PlotExclusion']='Comment5_9_unresolved_display_only'
s.loc[s.CellTypeDisplay.eq('Endothelial')&s.TargetedConflictFlag.eq('Immune_signal')&s.UMAP2.gt(-6.5),'PlotExclusion']='Comment7_endothelial_tail_display_only'
s.loc[s.CellTypeDisplay.eq('CD8 T')&s.UMAP1.between(-5,-3)&s.UMAP2.between(0,1),'PlotExclusion']='Comment8_isolated_T_display_only'
e=s[s.CellTypeDisplay.eq('Epithelial')]
dist=cKDTree(e[['UMAP1','UMAP2']]).query(e[['UMAP1','UMAP2']],k=6)[0][:,-1]
sel=e.UMAP1.between(-4,1)&e.UMAP2.between(-1.5,3)&(dist>.5)
s.loc[e.index[sel],'PlotExclusion']='Comment8_sparse_epithelial_bridge_display_only'
s['PlotIncluded']=s.PlotExclusion.eq('')
s['PlotMixedOpenSymbol']=s.MixedSignalDisplay&~s.PlotCellType.eq('Mono/macrophage')
assert len(s)==34227 and s.PlotIncluded.sum()==34104 and (~s.PlotIncluded).sum()==123
assert s.PlotCellType.eq('Mono/macrophage').sum()==3601
s.to_csv(R/'source_data/all_cells_display_decisions.csv.gz',index=False)
s[~s.PlotIncluded].to_csv(R/'source_data/display_exclusions.csv',index=False)
p=s[s.PlotIncluded];p.to_csv(R/'source_data/plotted_cells.csv.gz',index=False)
pal=pd.read_csv(P/'source_data/annotation_palette.csv')
pal=pal[~pal.CellTypeDisplay.isin(['Macrophage','Myeloid, other','Unresolved'])].copy()
pal.loc[pal.CellTypeDisplay.eq('Monocyte-like'),'Color']='#B77B48'
pal.CellTypeDisplay=pal.CellTypeDisplay.replace({'Monocyte-like':'Mono/macrophage'})
pal['NCells']=pal.CellTypeDisplay.map(p.PlotCellType.value_counts())
pal.to_csv(R/'source_data/annotation_palette.csv',index=False)
p.groupby(['PlotCellType','tissue']).size().rename('NCells').reset_index().to_csv(R/'source_data/display_counts.csv',index=False)
summary={'source_run':P.name,'source_csv_sha256':hashlib.sha256((P/'source_data/integrated_umap_cells.csv.gz').read_bytes()).hexdigest(),
 'original_cells':len(s),'plotted_cells':len(p),'omitted_display_only':123,'omission_reasons':s.loc[~s.PlotIncluded,'PlotExclusion'].value_counts().to_dict(),
 'plotted_tissue_counts':p.tissue.value_counts().to_dict(),'plotted_tissue_patients':p.groupby('tissue').patient.nunique().to_dict(),
 'patients':p.patient.nunique(),'libraries':p['sample'].nunique(),'merged_myeloid':3601,
 'remaining_open_symbols':int(p.PlotMixedOpenSymbol.sum()),'underlying_annotations_changed':False,'downstream_analysis_subset_changed':False}
(R/'source_data/revision_summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2))
