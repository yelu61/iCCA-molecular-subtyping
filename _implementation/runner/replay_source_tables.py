#!/usr/bin/env python3
"""Execute frozen table renderers in a separate tree; copied figures never count as redraws."""
from pathlib import Path
import argparse,hashlib,json,os,shutil,subprocess,time,re,tempfile
BUNDLE=Path(__file__).resolve().parents[1]
def sha(p):
 try:return hashlib.sha256(p.read_bytes()).hexdigest()
 except OSError:return ''
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--job',action='append');ap.add_argument('--original-root',type=Path,default=BUNDLE/'source_data/project');ap.add_argument('--timeout',type=int,default=180);ap.add_argument('--work-dir',type=Path);a=ap.parse_args()
 cfg=json.loads((BUNDLE/'code/replay_jobs.json').read_text());root=a.original_root.resolve()
 unknown=set(a.job or [])-{j['id'] for j in cfg['jobs']}
 if unknown:ap.error('Unknown replay job(s): '+', '.join(sorted(unknown)))
 out=BUNDLE/'review/production_replay';work=a.work_dir.resolve() if a.work_dir else Path(tempfile.mkdtemp(prefix='icca_table_replay_',dir='/private/tmp' if Path('/private/tmp').exists() else None));logs=out/'logs'
 if work==root or root in work.parents:ap.error('Replay work directory must be separate from the frozen input tree')
 work.mkdir(parents=True,exist_ok=True);logs.mkdir(parents=True,exist_ok=True)
 renv=Path(os.environ.get('ICCA_R_HOME',os.environ['R_HOME']));rex=renv/'bin/exec/R'
 env=os.environ.copy();env.update(R_HOME=str(renv),R_LIBS_USER=str(out/'disabled_user_library'),ICCA_ORIGINAL_ROOT=str(root),ICCA_REPLAY_ROOT=str(work),ICCA_BUNDLE_ROOT=str(BUNDLE),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1')
 statuses=[]
 for job in cfg['jobs']:
  if a.job and job['id'] not in a.job:continue
  env['ICCA_JOB_ID']=job['id']
  # Directory skeletons contain no pre-rendered figures.
  for d in job.get('directories',[]):(work/d).mkdir(parents=True,exist_ok=True)
  queue=[job['script']]+job.get('prefetch',[]);seen=set()
  while queue:
   p=queue.pop()
   if p in seen:continue
   seen.add(p)
   subprocess.run(['python3',str(BUNDLE/'code/materialize_input.py'),str(root),str(work),str(BUNDLE),p],capture_output=True,env=env)
   fp=work/p
   if fp.is_file() and fp.suffix in {'.R','.py'}:
    for s in re.findall(r'[\"\']([^\"\'\n]+)[\"\']',fp.read_text()):
     if s.startswith(('0-Data/','1-Clustering/','2-TranscriptomicCharacterization/','3-DuctType/','4-Multiomic/','5-ClassifierValidation/','6-ClinicalSurvival/','7-SingleCell/','8-Visualization/','governance/')):
      if (root/s).is_file():queue.append(s)
      elif (root/s).is_dir() or not Path(s).suffix:
       (work/s).mkdir(parents=True,exist_ok=True)
       for tail in ['figures','source_data','tables','review','review/pdf_rendered','main','supplement','independent','provenance','pdf','png','svg']:(work/s/tail).mkdir(parents=True,exist_ok=True)
  before={str(p.relative_to(work)):sha(p) for p in work.rglob('*.pdf')}
  runner=logs/(job['id']+'.R');runner.write_text('source('+json.dumps(str(BUNDLE/'code/replay_io.R'))+')\n'+job['command']+'\n')
  log=logs/(job['id']+'.log')
  if log.exists():shutil.copy2(log,logs/(job['id']+'_'+str(int(time.time()))+'.previous.log'))
  start=time.time();rc=None;error=''
  try:
   with log.open('w') as f:rc=subprocess.run([str(rex),'--vanilla','--slave','-f',str(runner)],cwd=work,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=a.timeout).returncode
  except subprocess.TimeoutExpired:error='timeout'
  except Exception as e:error=str(e)
  after={str(p.relative_to(work)):sha(p) for p in work.rglob('*.pdf')}
  copied={}
  inp=BUNDLE/'review/production_inputs.jsonl'
  if inp.exists():
   for line in inp.read_text().splitlines():
    r=json.loads(line)
    if r['kind']=='pre_rendered_asset':copied.setdefault(r['sha256'],[]).append(r['path'])
  changed=[{'path':p,'sha256':h,'status':'COPIED_EXISTING_FIGURE' if h in copied else 'GENERATED_FROM_EXECUTED_RENDERER'} for p,h in after.items() if before.get(p)!=h]
  for item in changed:
   p=work/item['path'];dest=out/'project'/item['path'];dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
   if p.with_suffix('.svg').exists():shutil.copyfile(p.with_suffix('.svg'),dest.with_suffix('.svg'))
  row={'job':job['id'],'script':job['script'],'return_code':rc,'error':error,'runtime_root':str(work),'seconds':round(time.time()-start,2),'log':str(log.relative_to(BUNDLE)),'outputs':changed,'status':'EXECUTED' if rc==0 else 'FAILED'}
  statuses.append(row);(out/(job['id']+'_status.json')).write_text(json.dumps(row,indent=2));print(json.dumps({'job':job['id'],'return_code':rc,'error':error,'new_pdfs':len(changed)}),flush=True)
 (out/'run_summary.json').write_text(json.dumps(statuses,indent=2))
 if not a.job:
  final_dir=work/'final_display_panels'
  rc=subprocess.run(['python3',str(BUNDLE/'code/replay_final_panels.py'),'--output',str(final_dir)],env=env).returncode
  if final_dir.exists():shutil.copytree(final_dir,BUNDLE/'review/production_final_panels',dirs_exist_ok=True,copy_function=shutil.copyfile)
  if rc:raise SystemExit(rc)
 if any(j['status']!='EXECUTED' for j in statuses):raise SystemExit(1)
if __name__=='__main__':main()
