"""Cache public NCBI gene metadata for the measured genes in displayed CNV bands."""
import csv
import hashlib
import json
import time
import urllib.request
import urllib.error
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / '4-Multiomic/refinements/fig4_cnv_gene_annotations_20261030/reference'
OUT.mkdir(parents=True, exist_ok=True)


def read(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f))


bands = {x['Cytoband'] for x in read(ROOT / '4-Multiomic/refinements/fig4_cnv_pairwise_20261030/pairwise_display.csv')}
genes = {x['Gene'] for x in read(ROOT / '4-Multiomic/refinements/fig4_cnv_20261030/gene_annotation_audit.csv')
         if x['Eligible'] == 'TRUE' and x['Cytoband'] in bands}
with open(ROOT / '0-Data/raw/2022Dong/WES_CNV_GISTIC2_raw.txt') as f:
    reader = csv.reader(f, delimiter='\t')
    next(reader)
    ids = sorted({x[1] for x in reader if x[0].split('|')[0].strip() in genes and int(x[1]) > 0}, key=int)
records = []
for i in range(0, len(ids), 150):
    chunk = ids[i:i + 150]
    url = 'https://api.ncbi.nlm.nih.gov/datasets/v2/gene/id/' + ','.join(chunk) + '?page_size=1000'
    path = OUT / f'ncbi_gene_batch_{i // 150 + 1:02d}.json'
    if not path.exists():
        request = urllib.request.Request(url, headers={'Accept': 'application/json', 'User-Agent': 'iCCA-region-gene-annotation/1.0'})
        try:
            with urllib.request.urlopen(request, timeout=45) as response:
                payload = response.read()
        except urllib.error.HTTPError as e:
            raise RuntimeError(e.read().decode()[:900]) from e
        json.loads(payload)
        path.write_bytes(payload)
        time.sleep(.4)
    data = json.loads(path.read_text())
    assert 'reports' in data
    records.append(dict(File=path.name, URL=url, RequestedIDs=';'.join(chunk),
                        Reports=len(data['reports']), MD5=hashlib.md5(path.read_bytes()).hexdigest()))
    print(f'Cached batch {i // 150 + 1}: {len(data["reports"])} gene records', flush=True)
with open(OUT / 'reference_manifest.csv', 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=list(records[0]))
    writer.writeheader()
    writer.writerows(records)
(OUT / 'retrieval.json').write_text(json.dumps(dict(provider='NCBI Datasets Gene API v2',
    retrieved_utc=datetime.now(timezone.utc).isoformat(), requested_ids=len(ids), source_gene_features=len(genes),
    source_bands=len(bands)), indent=2) + '\n')
