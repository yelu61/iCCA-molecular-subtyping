#!/usr/bin/env python3
"""Copy one requested project input into an isolated replay; never hydrate cloud files."""
import hashlib,json,os,shutil,sys,subprocess
from pathlib import Path
origin,work,bundle,requested=map(Path,sys.argv[1:5])
origin=origin.resolve();work=work.resolve();bundle=bundle.resolve()
if requested.is_absolute():
    try:rel=requested.relative_to(work)
    except ValueError:
        try:rel=requested.relative_to(origin)
        except ValueError:raise SystemExit(3)
else:rel=requested
if '..' in rel.parts:raise SystemExit(3)
src=origin/rel;dst=work/rel
if dst.exists():print(dst);raise SystemExit()
if not src.exists():raise SystemExit(2)
if src.is_dir():dst.mkdir(parents=True,exist_ok=True);print(dst);raise SystemExit()
st=src.stat()
if st.st_size>256_000_000:
    print('INPUT_LARGE_EXTERNAL '+str(rel),file=sys.stderr);raise SystemExit(4)
dst.parent.mkdir(parents=True,exist_ok=True)
try:
    subprocess.run([sys.executable,'-c','import shutil,sys;shutil.copyfile(sys.argv[1],sys.argv[2])',str(src),str(dst)],check=True,timeout=20)
except (subprocess.TimeoutExpired,subprocess.CalledProcessError):
    if dst.exists():dst.unlink()
    print('INPUT_READ_TIMEOUT '+str(rel),file=sys.stderr);raise SystemExit(4)
frozen=bundle/'source_data/project'/rel;frozen.parent.mkdir(parents=True,exist_ok=True)
if not frozen.exists():shutil.copyfile(dst,frozen)
row={'job':os.environ.get('ICCA_JOB_ID','pre_logging_run'),'path':str(rel),'bytes':st.st_size,'sha256':hashlib.sha256(src.read_bytes()).hexdigest(),'kind':'pre_rendered_asset' if src.suffix.lower() in {'.pdf','.png','.svg','.emf'} else 'script' if src.suffix in {'.R','.py','.mjs','.js'} else 'input'}
with (bundle/'review/production_inputs.jsonl').open('a') as f:f.write(json.dumps(row)+'\n')
print(dst)
