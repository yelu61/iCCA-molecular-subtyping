#!/usr/bin/env python3
"""Replay six final display corrections from frozen tables into a separate directory."""
from pathlib import Path
import argparse,json,os,subprocess,tempfile
B=Path(__file__).resolve().parents[1]
def main():
 p=argparse.ArgumentParser();p.add_argument('--output',type=Path);a=p.parse_args()
 out=a.output.resolve() if a.output else Path(tempfile.mkdtemp(prefix='icca_final_panels_'))
 if out==B or B in out.parents:p.error('Choose an output directory outside the frozen package')
 out.mkdir(parents=True,exist_ok=True);env=os.environ.copy()
 cfg=json.loads((B/'code/replay_jobs.json').read_text());rhome=Path(env.get('ICCA_R_HOME',os.environ['R_HOME']))
 env.update(R_HOME=str(rhome),R_LIBS_USER=str(out/'disabled_user_library'),ICCA_BUNDLE_ROOT=str(B),ICCA_FINALIZATION_OUTPUT_ROOT=str(out))
 rows=[]
 for source in sorted((B/'code/final_panel_drivers').glob('*.R')):
  log=out/(source.stem+'.log')
  with log.open('w') as f:rc=subprocess.run([str(rhome/'bin/exec/R'),'--vanilla','--slave','-f',str(source)],env=env,cwd=out,stdout=f,stderr=subprocess.STDOUT).returncode
  rows.append({'driver':str(source.relative_to(B)),'return_code':rc,'log':str(log)})
 result={'output_root':str(out),'drivers':rows,'pdfs':[str(p.relative_to(out)) for p in (out/'figures/panels').glob('*.pdf')]}
 (out/'replay_summary.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
 return int(any(r['return_code']!=0 for r in rows))
if __name__=='__main__':raise SystemExit(main())
