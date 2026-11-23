from pathlib import Path
import pandas as pd,numpy as np,shutil,json,hashlib
r=Path.cwd();o=r/'8-Visualization/refinements/sfig4_cnv_D_20261104_v1';s=r/'7-SingleCell/runs/20261103_figure3_cnv_program_cases_v1/tables';ins={}
for n in ['ICC20_cnv.csv.gz','ICC20_rows.csv','ICC20_windows.csv','cell_support.csv']:
 f=s/n;shutil.copy2(f,o/'source_data'/n);ins[str(f.relative_to(r))]=hashlib.sha256(f.read_bytes()).hexdigest()
sup=pd.read_csv(s/'cell_support.csv');rows=pd.read_csv(s/'ICC20_rows.csv');m=pd.read_csv(s/'ICC20_cnv.csv.gz',header=None);win=pd.read_csv(s/'ICC20_windows.csv');out=[]
for patient,group,label in [('ICC20','Tumor epithelium','ICC20'),('ICC18','Tumor epithelium','ICC18'),('ICC24','Tumor1 epithelium','ICC24 / T1'),('ICC24','Tumor2 epithelium','ICC24 / T2')]:
 a=sup[(sup.Patient==patient)&(sup.Group==group)].iloc[0];out.append({'Patient':patient,'Region':label,'CleanEpithelial':int(a.NCleanTumorEpithelial),'StrictCandidates':int(a.NCandidates)})
pd.DataFrame(out).to_csv(o/'source_data/patient_support.csv',index=False)
checks={'all2265savedrows':len(rows)==2265 and len(m)==2265,'all22autosomes':win.Chromosome.nunique()==22,'windowaligned':len(win)==m.shape[1],'allfinite':bool(np.isfinite(m).all().all()),'27ICC20candidates':int(rows.StrictCandidate.sum())==27,'originalorderretained':rows.Row.tolist()==list(range(1,2266)),'42candidatesacross3patients':sum(x['StrictCandidates']for x in out)==42,'ICC24T2zero':out[-1]['StrictCandidates']==0}
(o/'review/source_verification.json').write_text(json.dumps({'checks':checks,'all_pass':all(checks.values()),'matrix_shape':m.shape,'source_sha256':ins,'above_display_range_fraction':float((abs(m.to_numpy())>.3).mean())},indent=2));assert all(checks.values())
print('Verified',m.shape,'matrix;42candidates across3patients;all2265ICC20rows retained.')
