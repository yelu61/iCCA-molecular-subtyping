"""Independent metrics/list overlap checks and a compact nine-setting result table."""
from pathlib import Path
from collections import defaultdict
import csv,json,math,statistics
ROOT=Path(__file__).resolve().parents[2]
OUT=ROOT/'5-ClassifierValidation/runs/bulk-rnaseq__20261031_step5_threshold_grid_v1'
PARENT=ROOT/'5-ClassifierValidation/runs/bulk-rnaseq__20261030_step5_consensus_v1'
CL=['C1','C2','C3'];checks={}
def read(p):
    with Path(p).open() as f:return list(csv.DictReader(f))
def check(n,ok):
    checks[n]=bool(ok)
    if not ok:raise AssertionError(n)
def auc(y,s):
    pairs=sorted(zip(s,y));total=0.;i=0
    while i<len(pairs):
        j=i+1
        while j<len(pairs) and pairs[j][0]==pairs[i][0]:j+=1
        total+=sum(p[1] for p in pairs[i:j])*(i+1+j)/2;i=j
    pos=sum(y);return (total-pos*(pos+1)/2)/(pos*(len(y)-pos))
def metrics(d):
    recall=[sum(r['Truth']==c and r['Prediction']==c for r in d)/sum(r['Truth']==c for r in d) for c in CL]
    f1=[2*sum(r['Truth']==c and r['Prediction']==c for r in d)/(sum(r['Truth']==c for r in d)+sum(r['Prediction']==c for r in d)) for c in CL]
    result={'Accuracy':sum(r['Truth']==r['Prediction'] for r in d)/len(d),'BalancedAccuracy':sum(recall)/3,'MacroF1':sum(f1)/3,'MacroAUROC':sum(auc([int(r['Truth']==c) for r in d],[float(r['Probability'+c]) for r in d]) for c in CL)/3}
    result.update({c+'_Recall':v for c,v in zip(CL,recall)});return result
specs=read(OUT/'tables/specifications.csv');names=[r['Spec'] for r in specs]
raw=read(OUT/'tables/predictions.csv');oof=read(OUT/'tables/oof_predictions.csv');perf=read(OUT/'tables/performance_ci.csv');delta=read(OUT/'tables/paired_differences.csv')
check('complete raw/aggregate predictions',len(raw)==9*497*5 and len(oof)==9*497)
groups=defaultdict(list)
for r in raw:groups[(r['Spec'],r['Sample'])].append(r)
for r in oof:
    z=groups[(r['Spec'],r['Sample'])];assert len(z)==5 and len({t['Repeat'] for t in z})==5
    p=[sum(float(t['Probability'+c]) for t in z)/5 for c in CL]
    assert max(abs(p[i]-float(r['Probability'+c])) for i,c in enumerate(CL))<1e-12
    assert r['Prediction']==CL[max(range(3),key=lambda i:p[i])]
check('five-repeat aggregation and labels',True)
check('finite normalized probabilities',all(all(math.isfinite(float(r['Probability'+c])) and 0<=float(r['Probability'+c])<=1 for c in CL) and abs(sum(float(r['Probability'+c]) for c in CL)-1)<1e-10 for r in raw+oof))
computed={sp:metrics([r for r in oof if r['Spec']==sp]) for sp in names}
for sp in names:
    p=[r for r in perf if r['Spec']==sp]
    check(sp+' seven metrics independently reproduced',len(p)==7 and all(abs(computed[sp][r['Metric']]-float(r['Estimate']))<1e-12 for r in p))
check('paired differences direction and estimates',len(delta)==8*7 and all(abs(computed[r['Spec']][r['Metric']]-computed['L30_RF200'][r['Metric']]-float(r['Difference']))<1e-12 for r in delta))
parent={r['Sample']:r for r in read(PARENT/'tables/nested_oof_predictions.csv') if r['Model']=='Consensus_SVM'}
base={r['Sample']:r for r in oof if r['Spec']=='L30_RF200'}
err=max(abs(float(base[s]['Probability'+c])-float(parent[s]['Probability'+c])) for s in base for c in CL)
check('reference reproduces parent probabilities',err<1e-10 and all(base[s]['Prediction']==parent[s]['Prediction'] for s in base))
full=read(OUT/'tables/full_development_features.csv');fs={sp:{r['Gene'] for r in full if r['Spec']==sp} for sp in names};ref=fs['L30_RF200']
check('reference full panel has114',len(ref)==114)
for r in read(OUT/'tables/full_development_overlap.csv'):
    g=fs[r['Spec']];shared=len(g&ref);union=len(g|ref)
    assert len(g)==int(r['N']) and shared==int(r['Shared']) and union==int(r['Union']) and abs(shared/union-float(r['Jaccard']))<1e-12
check('full-development overlaps independently reproduced',True)
outer=read(OUT/'tables/features.csv');sets=defaultdict(set)
for r in outer:sets[(r['Context'],r['Spec'])].add(r['Gene'])
overlap=read(OUT/'tables/outer_reference_overlap.csv')
for r in overlap:
    g=sets[(r['Context'],r['Spec'])];b=sets[(r['Context'],'L30_RF200')]
    assert len(g)==int(r['N']) and abs(len(g&b)/len(g|b)-float(r['Jaccard']))<1e-12
check('225 within-fold overlaps reproduced',len(overlap)==225)
freq=read(OUT/'tables/feature_selection_frequency.csv');freq_counts=defaultdict(int)
for (ctx,sp),g in sets.items():
    for gene in g:freq_counts[(sp,gene)]+=1
check('outer selection frequencies',len(freq)==len(freq_counts) and all(int(r['OuterSelections'])==freq_counts[(r['Spec'],r['Gene'])] and abs(float(r['Fraction'])-freq_counts[(r['Spec'],r['Gene'])]/25)<1e-12 for r in freq))
check('selector/membership audit passed',json.loads((OUT/'review/selection_verification.json').read_text())['passed'])
check('analysis completed',json.loads((OUT/'status.json').read_text())['state']=='COMPUTED_REVIEW')
summary=[]
for sp in specs:
    name=sp['Spec'];p={r['Metric']:r for r in perf if r['Spec']==name};dd={r['Metric']:r for r in delta if r['Spec']==name}
    ns=[len(g) for (ctx,s),g in sets.items() if s==name];js=[float(r['Jaccard']) for r in overlap if r['Spec']==name];g=fs[name]
    b=p['BalancedAccuracy'];db=dd.get('BalancedAccuracy',{'Difference':0,'Lower95':0,'Upper95':0})
    row=dict(sp,FullGenes=len(g),FullJaccard=len(g&ref)/len(g|ref),Added=len(g-ref),Removed=len(ref-g),OuterMedian=statistics.median(ns),OuterMin=min(ns),OuterMax=max(ns),OuterJaccardMedian=statistics.median(js),OuterJaccardMin=min(js),OuterJaccardMax=max(js),BalancedAccuracy=float(b['Estimate']),BALower95=float(b['Lower95']),BAUpper95=float(b['Upper95']),BADifference=float(db['Difference']),BADifferenceLower95=float(db['Lower95']),BADifferenceUpper95=float(db['Upper95']),MacroF1=computed[name]['MacroF1'],MacroAUROC=computed[name]['MacroAUROC'],OOFAssignmentDifferences=sum(r['Prediction']!=base[r['Sample']]['Prediction'] for r in oof if r['Spec']==name))
    row.update({c+'Recall':computed[name][c+'_Recall'] for c in CL});summary.append(row)
with (OUT/'tables/threshold_summary.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(summary[0]));w.writeheader();w.writerows(summary)
(OUT/'review/independent_metrics.json').write_text(json.dumps({'passed':all(checks.values()),'checks':checks,'reference_probability_error':err,'bootstrap_scope':'Point estimates and paired differences independently reproduced. R conditional bootstrap endpoints not independently resampled; pointwise exploratory intervals.'},indent=2))
print(f'PASS: {len(checks)} metric/list completeness checks; reference probability error {err:.3g}')
for r in summary:print(r['Spec'],r['FullGenes'],f"J={r['FullJaccard']:.3f}",f"BA={100*r['BalancedAccuracy']:.3f}%",f"delta={100*r['BADifference']:+.3f}pp")
