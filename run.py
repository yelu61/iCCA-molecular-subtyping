#!/usr/bin/env python3
"""Canonical public entry: verify inputs, redraw selected figures, or stage upstream analysis."""
from pathlib import Path
import argparse,json,hashlib,os,shutil,subprocess,sys,tempfile,csv
REPO=Path(__file__).resolve().parent

def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def copy(a,b):b.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(a,b)
def verify(data):
 records=json.loads((REPO/'config/input_manifest.json').read_text());bad=[r['path']for r in records if not(data/r['path']).is_file()or digest(data/r['path'])!=r['sha256']]
 if bad:raise SystemExit('Input missing/hash mismatch: '+', '.join(bad))
 return len(records)
def stage_implementation(dest):
 for p in(REPO/'_implementation/project').rglob('*'):
  if not p.is_file():continue
  target=dest/p.relative_to(REPO/'_implementation/project');copy(p,target)
  if target.suffix in{'.R','.py','.ipynb','.json','.yaml','.yml','.txt'}:
   text=target.read_text(errors='replace').replace('@PROJECT_ROOT@',str(dest));target.write_text(text)
def rhome(value):
 if value:return Path(value).resolve()
 if os.environ.get('ICCA_R_HOME'):return Path(os.environ['ICCA_R_HOME']).resolve()
 return Path(subprocess.check_output(['R','RHOME'],text=True).strip()).resolve()
def main():
 p=argparse.ArgumentParser(description=__doc__);sub=p.add_subparsers(dest='command',required=True)
 v=sub.add_parser('verify');v.add_argument('--data-root',type=Path,required=True)
 d=sub.add_parser('figures');d.add_argument('--data-root',type=Path,required=True);d.add_argument('--output',type=Path,required=True);d.add_argument('--r-home');d.add_argument('--job',action='append')
 a=sub.add_parser('analysis');a.add_argument('--list',action='store_true');a.add_argument('--stage');a.add_argument('--project-root',type=Path);a.add_argument('--execute',action='store_true');a.add_argument('--r-home')
 args=p.parse_args()
 if args.command=='verify':print(json.dumps({'input_files_verified':verify(args.data_root.resolve())}));return
 if args.command=='analysis':
  catalog=json.loads((REPO/'config/analysis_stages.json').read_text())
  if args.list or not args.stage:
   print(json.dumps(catalog,indent=2));return
  if args.stage not in catalog:p.error('Unknown stage')
  if not args.project_root:p.error('--project-root must name independently prepared inputs and code')
  root=args.project_root.resolve()
  if root==REPO or REPO in root.parents or any((x/'.git').exists()for x in[root,*root.parents]):p.error('Choose an independent non-Git project directory')
  root.mkdir(parents=True,exist_ok=True);stage_implementation(root);item=catalog[args.stage];entry=root/item['entry'];kind=item['kind']
  if kind=='R':rh=rhome(args.r_home);cmd=[str(rh/'bin/exec/R'),'--vanilla','--slave','-f',str(entry)]+(['--args','run']if args.stage in{'classifier114','compact'}else[])
  elif kind=='python':cmd=[sys.executable,str(entry)]
  else:cmd=['jupyter','nbconvert','--to','notebook','--execute',str(entry),'--output',entry.stem+'_executed.ipynb']
  print(json.dumps({'stage':args.stage,'command':cmd,'executed':args.execute,'full_pipeline_verified':False}))
  if args.execute:
   env=os.environ.copy()
   if kind=='R':env.update(R_HOME=str(rh),R_LIBS_USER=str(root/'user_R_library'))
   raise SystemExit(subprocess.run(cmd,cwd=root,env=env).returncode)
  return
 data=args.data_root.resolve();out=args.output.resolve()
 if out==data or data in out.parents or out==REPO or REPO in out.parents:p.error('Output must be outside data and code trees')
 if out.exists()and any(out.iterdir()):p.error('Output must be new or empty')
 n=verify(data);out.mkdir(parents=True,exist_ok=True);session=out/'session';session.mkdir()
 # Inputs remain immutable. Only session copies are materialized by historical renderers.
 shutil.copytree(data/'source_data',session/'source_data',copy_function=shutil.copyfile)
 stage_implementation(session/'source_data/project')
 for q in(REPO/'_implementation/panel_sources').rglob('*'):
  if q.is_file():copy(q,session/'source_data'/q.relative_to(REPO/'_implementation/panel_sources'))
 for q in(REPO/'_implementation/runner').glob('*'):
  if q.is_file():copy(q,session/'code'/q.name)
 for q in(REPO/'_implementation/final_panels').glob('*.R'):copy(q,session/'code/final_panel_drivers'/q.name)
 copy(REPO/'config/figure_jobs.json',session/'code/replay_jobs.json');(session/'review').mkdir()
 rh=rhome(args.r_home);env=os.environ.copy();env.update(R_HOME=str(rh),ICCA_R_HOME=str(rh),R_LIBS_USER=str(out/'disabled_user_library'),ICCA_IGNORE_RUNTIME_PATHS='1')
 cmd=[sys.executable,str(session/'code/replay_source_tables.py'),'--original-root',str(session/'source_data/project'),'--work-dir',str(out/'working_tree'),'--timeout','300']
 for j in args.job or[]:cmd+=['--job',j]
 rc=subprocess.run(cmd,env=env).returncode
 report={'return_code':rc,'input_files_verified':n,'selected_jobs':args.job or'ALL37_PLUS3_FINAL_DRIVERS','output_root':str(out),'upstream_recomputed':False}
 (out/'run_report.json').write_text(json.dumps(report,indent=2));print(json.dumps(report));raise SystemExit(rc)
if __name__=='__main__':main()
