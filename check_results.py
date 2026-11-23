#!/usr/bin/env python3
"""Check all selected quantitative outputs against accepted frozen-table redraw content."""
from pathlib import Path
import argparse,json,subprocess,re,csv,hashlib,shutil,xml.etree.ElementTree as E
B=Path(__file__).resolve().parent

def xml_content(p):
 r=E.fromstring(p.read_bytes());return {'text':[''.join(x.itertext())for x in r.iter()if x.tag.endswith('}text')],'images':[hashlib.sha256((x.get('href',x.get('{http://www.w3.org/1999/xlink}href'))or'').encode()).hexdigest()for x in r.iter()if x.tag.endswith('}image')],'viewBox':r.get('viewBox')}
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--run-dir',type=Path,required=True);ap.add_argument('--data-root',type=Path,required=True);a=ap.parse_args();session=a.run_dir/'session';refs=json.loads((B/'config/expected_panel_content.json').read_text());checks=[]
 target=a.run_dir/'selected_panels';target.mkdir(exist_ok=True)
 for r in refs:
  p=session/r['output'];ok=p.exists();detail=''
  if ok and p.suffix=='.pdf':
   text=re.sub(r'\s+',' ',subprocess.check_output(['pdftotext','-layout',str(p),'-'],text=True)).strip();info=subprocess.check_output(['pdfinfo',str(p)],text=True);d=re.search(r'Page size:\s+([\d.]+) x ([\d.]+)',info);dims=[round(float(x),1)for x in d.groups()]if d else[];ok=text==r['text']and dims==r['dimensions'];detail='Exact normalized redraw text and dimensions'
  elif ok:
   # SVG serialization IDs may vary; extract meaningful texts/images/viewport.
   ok=xml_content(p)==r['svg_content'];detail='Exact SVG text/image payload/viewBox'
  if ok:
   dest=target/r['figure']/('Panel_'+r['panel']+p.suffix);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
  checks.append({'asset_id':r['asset_id'],'figure':r['figure'],'panel':r['panel'],'pass':ok,'check':detail,'source':str(p)})
 prefix=Path('8-Visualization/refinements/fig1_sfig1_typography_20261030/source_data');z=a.data_root/'source_data/derived_render_tables'/prefix/'Figure1B_row_zscore.csv';rr=z.with_name('Figure1B_rows.csv');cc=z.with_name('Figure1B_columns.csv')
 with z.open()as f:reader=csv.DictReader(f);fields=reader.fieldnames;rows=list(reader)
 with rr.open()as f:genes=list(csv.DictReader(f))
 with cc.open()as f:samples=list(csv.DictReader(f))
 matrix_ok=len(rows)==2000 and len(fields)==587 and fields[1:]==[x['Sample']for x in samples]and [x['Gene']for x in rows]==[x['Gene']for x in genes]
 status=json.loads((session/'review/production_replay/fig1_legacy_status.json').read_text());matrix_ok=matrix_ok and status['return_code']==0
 checks.append({'asset_id':'S01_P02','figure':'Figure1','panel':'B','pass':matrix_ok,'check':'2000x586 complete matrix and row/column identities; original renderer success'})
 report={'quantitative_assets':len(checks),'passed':sum(x['pass']for x in checks),'failures':[x for x in checks if not x['pass']],'scope':'Current 100 quantitative panel assets;2 workflow schematics and14-program editable synthesis illustration are separate; no scientific acceptance or upstream refit inferred','checks':checks}
 (a.run_dir/'selected_panel_verification.json').write_text(json.dumps(report,indent=2));print(json.dumps({k:v for k,v in report.items()if k!='checks'}));return int(bool(report['failures']))
if __name__=='__main__':raise SystemExit(main())
