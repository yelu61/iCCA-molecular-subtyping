"""Conservative, response-blind name/alias crosswalk. Never changes raw data."""
import csv, hashlib, json, re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '5-ClassifierValidation/review/20261101_drug_evidence_refresh_v1'
(OUT / 'tables').mkdir(parents=True, exist_ok=True)
inputs = []
def read(p):
    p = ROOT / p; inputs.append(p)
    with p.open(newline='') as f: return list(csv.DictReader(f))
def write(name, rows, fields=None):
    with (OUT / 'tables' / (name + '.csv')).open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0])); w.writeheader(); w.writerows(rows)
def norm(s):
    # Keep stereochemical + / - labels; do not strip them into a false identity.
    return re.sub(r'[\s._-]', '', s.upper())

drugs = read('5-ClassifierValidation/publication/tables/step5_oncopredict_gdsc2_cv_quality.csv')
meta = read('oncoPredict_data/GLDS/GDSCv2/screened_compunds_rel_8.2.csv')
md = {r['DRUG_ID']: r for r in meta}
prpath = ROOT / '5-ClassifierValidation/Drug_sensitivity_AUC_(PRISM_Repurposing_Secondary_Screen).csv'
inputs.append(prpath)
with prpath.open() as f: assays = next(csv.reader(f))[1:]
assaymap = {}
for a in assays: assaymap.setdefault(a.split(' (BRD:')[0], []).append(a)
curated = {
    '5-Fluorouracil': [('FLUOROURACIL', 'https://www.cancer.gov/publications/dictionaries/cancer-drug/def/fluorouracil')],
    'Selumetinib': [('AZD6244', 'https://www.cancer.gov/publications/dictionaries/cancer-drug/def/selumetinib'), ('ARRY-142886', 'https://www.cancer.gov/publications/dictionaries/cancer-drug/def/selumetinib')],
}
blocked = {('Selumetinib', 'AZD1480'): 'Incorrect resource synonym: AZD6244 is selumetinib; AZD1480 is a separate drug.'}
salt_hold = {'Obatoclax Mesylate', 'Vinorelbine', 'Fludarabine'}
rows, inventory, targets = [], [], []
for d in drugs:
    drug = d['Drug']; name, drugid = drug.rsplit('_', 1); m = md[drugid]
    aliases = [(name, 'Exact or punctuation-normalized model name', 'GDSC2 model column'),
               (m['DRUG_NAME'], 'Resource drug name', 'GDSC release 8.2')]
    aliases += [(x.strip(), 'Resource synonym', 'GDSC release 8.2') for x in m['SYNONYMS'].split(',') if x.strip()]
    # Repeated GDSC model IDs of the same named compound share identity aliases,
    # while their model scores remain separate throughout the analysis.
    for other in drugs:
        oname, oid = other['Drug'].rsplit('_', 1)
        if norm(oname)==norm(name) and oid!=drugid:
            aliases += [(x.strip(), 'Resource synonym from same compound model', 'GDSC release 8.2; same named compound ID '+oid) for x in md[oid]['SYNONYMS'].split(',') if x.strip()]
    aliases += [(x, 'Curated primary-source synonym', u) for x, u in curated.get(name, [])]
    found = {}
    for alias, kind, src in aliases:
        for compound in assaymap:
            if norm(alias) != norm(compound): continue
            status, reason = 'Accepted name/alias link', 'Not independently verified by chemical structure'
            if (name, alias) in blocked: status, reason = 'Rejected erroneous alias', blocked[name, alias]
            elif name in salt_hold and norm(compound) != norm(name):
                status, reason = 'Held formulation/active-moiety link', 'Salt/formulation equivalence unresolved; not pooled'
            elif '(+)' in compound and '(+)' not in name:
                status, reason = 'Held stereochemistry link', 'Stereochemistry requires explicit structural confirmation'
            item = (kind, src, alias, status, reason)
            if compound not in found or kind == 'Curated primary-source synonym': found[compound] = item
    for compound, (kind, src, alias, status, reason) in found.items():
        for assay in assaymap[compound]: rows.append(dict(Drug=drug,ModelCompound=name,GDSCDrugID=drugid,PRISMCompound=compound,Assay=assay,MatchedAlias=alias,MatchType=kind,Source=src,Status=status,Qualification=reason))
    accepted = [c for c, v in found.items() if v[3] == 'Accepted name/alias link']
    inventory.append(dict(Drug=drug,ModelCompound=name,AcceptedCompoundNames=len(accepted),AcceptedAssays=sum(len(assaymap[c]) for c in accepted),HeldOrRejectedNames=sum(v[3]!='Accepted name/alias link' for v in found.values()),CoverageStatus='Mapped names' if accepted else 'No accepted name link; measurement availability unresolved'))
    targets.append(dict(Drug=drug,ModelCompound=name,GDSCDrugID=drugid,SourceTarget=m['TARGET'],SourcePathway=m['TARGET_PATHWAY'],SourceSynonyms=m['SYNONYMS']))

# Same PRISM assay linked to different normalized model compounds is ambiguous.
byassay = {}
for r in rows:
    if r['Status']=='Accepted name/alias link': byassay.setdefault(r['Assay'],set()).add(norm(r['ModelCompound']))
for r in rows:
    if len(byassay.get(r['Assay'],set())) > 1:
        r['Status']='Held ambiguous cross-model alias';r['Qualification']='Multiple distinct model compounds map to the same assay'
for inv in inventory:
    rr=[r for r in rows if r['Drug']==inv['Drug'] and r['Status']=='Accepted name/alias link']
    inv['AcceptedCompoundNames']=len({r['PRISMCompound'] for r in rr});inv['AcceptedAssays']=len(rr)
    inv['CoverageStatus']='Mapped names' if rr else 'No accepted name link; measurement availability unresolved'
write('chemical_crosswalk_all_links',rows)
write('chemical_model_coverage',inventory)
write('drug_target_resource_annotations',targets)
write('identity_input_hashes',[dict(Path=str(p.relative_to(ROOT)),SHA256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in inputs])
summary=dict(models=len(drugs),accepted_models=sum(r['AcceptedAssays']>0 for r in inventory),accepted_unique_assays=len({r['Assay'] for r in rows if r['Status']=='Accepted name/alias link'}),held_or_rejected_links=sum(r['Status']!='Accepted name/alias link' for r in rows),structural_completeness_claim=False)
(OUT/'identity_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary))
