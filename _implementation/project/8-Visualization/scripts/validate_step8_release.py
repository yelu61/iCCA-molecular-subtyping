#!/usr/bin/env python3
"""Validate standalone PDF source contracts; a structural PASS is not visual QA."""
from __future__ import annotations
import argparse,csv,json,re,shutil,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[2]
ALLOWED_STATUS={"SourceLocked","AnalysisReady","Blocked","Rendered","Finalized"}

def split_paths(value): return [x.strip() for x in value.split(';') if x.strip()]

def source_fingerprint(root,sources):
    import hashlib
    return ';'.join(f'{s}={hashlib.md5((root/s).read_bytes()).hexdigest()}' for s in sources)

def checked_path(root,rel):
    p=(root/rel).resolve()
    if not p.is_relative_to(root.resolve()): raise ValueError(f'Path escapes project: {rel}')
    return p

def check_step7_release(root):
    """Require explicit current release identity and source hashes, never old gate alone."""
    p=root/'7-SingleCell/publication/tables/step7_current_release_contract.json'
    if not p.is_file(): raise ValueError('Step7 integrated reanalysis is not released; current release contract absent')
    record=json.loads(p.read_text())
    if record.get('DesignID')!='Step7_Integrated_Reanalysis_Design_20260909' or record.get('TechnicalRelease')!='PASS':
        raise ValueError('Step7 current release identity/gate is not accepted')
    if record.get('ScientificInterpretation') not in ('QUALIFIED','PASS'): raise ValueError('Step7 scientific interpretation missing')
    if not record.get('RunID') or not record.get('SourceMD5'): raise ValueError('Step7 run identity/source hashes missing')
    import hashlib
    for rel,expected in record['SourceMD5'].items():
        src=checked_path(root,rel)
        if not src.is_file() or hashlib.md5(src.read_bytes()).hexdigest()!=expected: raise ValueError(f'Stale Step7 source: {rel}')
    if record.get('ReleaseGate') not in record['SourceMD5']: raise ValueError('Step7 release gate is not covered by source hashes')
    gate=checked_path(root,record['ReleaseGate'])
    with gate.open(newline='') as h: rows=list(csv.DictReader(h))
    if not rows: raise ValueError('Empty Step7 gate')
    for r in rows:
        if str(r.get('Blocking','true')).lower() in ('true','1') and r.get('Status','').upper() not in ('PASS','QUALIFIED'):
            raise ValueError('Step7 has a failed blocking gate')
    return record

def validate(root,figure=None,require_final=False):
    errors=[];warnings=[];pending=blocked=0
    root=root.resolve();manifest=root/'8-Visualization/tables/figure_manifest.csv'
    for sub in ('scripts','tables','figures/main','figures/supplement','reports'):
        if not (root/'8-Visualization'/sub).is_dir():errors.append(f'Missing directory: {sub}')
    try:
        with manifest.open(newline='',encoding='utf-8-sig') as h:
            reader=csv.DictReader(h)
            required={'Figure','Panel','SourceTable','TargetOutput','Status','Notes'}
            if not required.issubset(reader.fieldnames or []):errors.append('Manifest missing required columns')
            all_rows=list(reader)
    except OSError as e:return [str(e)],[],0,0,0
    rows=[r for r in all_rows if figure is None or r['Figure']==figure or (figure.startswith('S') and r['Figure']=='Supplement' and r['Panel'].startswith(figure) and (r['Panel']==figure or r['Panel'][len(figure):len(figure)+1].isalpha()))]
    if not rows:errors.append(f'No manifest rows selected: {figure}')
    keys=set();targets=set()
    for r in all_rows:
        key=(r['Figure'],r['Panel']);target=r['TargetOutput']
        if key in keys:errors.append(f'Duplicate figure ID: {key}')
        if target in targets:errors.append(f'Duplicate target: {target}')
        keys.add(key);targets.add(target)
    for r in rows:
        label=f"{r['Figure']} {r['Panel']}";status=r['Status'];sources=split_paths(r['SourceTable'])
        if status not in ALLOWED_STATUS:errors.append(f'Unknown status: {label}: {status}')
        is_step7=r.get('AnalysisOwner')=='Step7'
        if status=='Blocked':
            blocked+=1
            if not r.get('Notes'):errors.append(f'Blocked reason missing: {label}')
            if require_final:errors.append(f'Blocked: {label}')
        elif is_step7:
            try:record=check_step7_release(root)
            except (ValueError,KeyError,OSError) as e:errors.append(str(e));record=None
            if record and any(s not in record['SourceMD5'] for s in sources):errors.append(f'Step7 source not covered by release hashes: {label}')
        if not sources:errors.append(f'No source: {label}')
        for rel in sources:
            try:
                p=checked_path(root,rel)
                if not p.is_file() or not p.stat().st_size:errors.append(f'Missing source: {rel}')
                elif p.suffix.lower()=='.csv':
                    with p.open(newline='',encoding='utf-8-sig') as h:
                        if next(csv.DictReader(h),None) is None:errors.append(f'Header-only CSV: {rel}')
            except (ValueError,OSError) as e:errors.append(str(e))
        try:
            target=checked_path(root,r['TargetOutput'])
            if not r['TargetOutput'].startswith(('8-Visualization/figures/main/','8-Visualization/figures/supplement/')) or target.suffix.lower()!='.pdf':
                errors.append(f'Invalid standalone PDF target: {label}');continue
        except ValueError as e:errors.append(str(e));continue
        if not target.is_file():
            pending+=1
            if status in ('Rendered','Finalized') or (require_final and status!='Blocked'):errors.append(f'Missing PDF: {label}')
            continue
        if status=='Blocked':errors.append(f'Blocked source has an active publication PDF: {label}');continue
        with target.open('rb') as h:
            if h.read(5)!=b'%PDF-':errors.append(f'Invalid PDF header: {label}')
        pdfinfo=shutil.which('pdfinfo')
        if pdfinfo:
            proc=subprocess.run([pdfinfo,str(target)],capture_output=True,text=True)
            pages=re.search(r'^Pages:\s+(\d+)',proc.stdout,re.M)
            size=re.search(r'^Page size:\s+([\d.]+) x ([\d.]+)',proc.stdout,re.M)
            if proc.returncode or not pages or int(pages[1])!=1:errors.append(f'PDF must contain one page: {label}')
            if size and r.get('WidthMM') and r.get('HeightMM'):
                actual=[float(size[i])*25.4/72 for i in (1,2)]
                if any(abs(a-float(b))>.36 for a,b in zip(actual,[r['WidthMM'],r['HeightMM']])):errors.append(f'PDF dimensions differ from manifest: {label}')
        elif require_final:errors.append('pdfinfo unavailable: final page/size check cannot run')
        else:warnings.append('pdfinfo unavailable; page/size checks skipped')
        if not r.get('SourceMD5'):errors.append(f'Rendered PDF lacks source fingerprint: {label}')
        elif all((root/s).is_file() for s in sources) and r['SourceMD5']!=source_fingerprint(root,sources):errors.append(f'Source changed after rendering: {label}')
        if require_final and status!='Finalized':errors.append(f'Visual QA not finalized: {label}')
        if status=='Finalized':
            try:
                qa=json.loads(checked_path(root,r.get('QAFile','')).read_text())
                import hashlib
                if qa.get('VisualQA')!='PASS' or qa.get('SourceMD5')!=r['SourceMD5'] or qa.get('PDFMD5')!=hashlib.md5(target.read_bytes()).hexdigest():errors.append(f'Missing/stale visual QA: {label}')
            except (ValueError,OSError,KeyError):errors.append(f'Invalid visual QA record: {label}')
    return errors,sorted(set(warnings)),len(rows),pending,blocked

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--figure',help='Figure 1, Figure 4, Supplement, or S10, etc.')
    p.add_argument('--require-final',action='store_true',help='Require selected standalone PDFs and matching visual-QA records')
    p.add_argument('--root',type=Path,default=ROOT)
    a=p.parse_args();e,w,n,pending,blocked=validate(a.root,a.figure,a.require_final)
    print('Step8 standalone PDF validation: '+('FAIL' if e else 'PASS (structural/source checks only)'))
    print(f'- Selected rows: {n}; pending PDFs: {pending}; blocked: {blocked}')
    for x in e:print('- '+x)
    for x in w:print('- WARNING: '+x)
    return bool(e)
if __name__=='__main__':sys.exit(main())
