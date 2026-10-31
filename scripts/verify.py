#!/usr/bin/env python3
from pathlib import Path
from html.parser import HTMLParser
from urllib.parse import urlsplit, parse_qs, unquote
import hashlib,json
ROOT=Path(__file__).resolve().parents[1]
class Links(HTMLParser):
 def __init__(self):super().__init__();self.links=[]
 def handle_starttag(self,tag,attrs):
  d=dict(attrs)
  for k in ['href','src']:
   if k in d:self.links.append(d[k])
errors=[];checked=0;run_links=0
manifest=json.loads((ROOT/'studies/morris-2022/viewer/manifest.json').read_text());ids={r['id'] for r in manifest['runs']}
for p in ROOT.rglob('*.html'):
 doc=Links();doc.feed(p.read_text())
 for link in doc.links:
  u=urlsplit(link)
  if u.scheme or u.netloc or not u.path:continue
  if u.path.startswith('/'):errors.append('Root-relative path: '+link)
  target=(p.parent/unquote(u.path)).resolve()
  if target.is_dir():target=target/'index.html'
  if not target.is_file():errors.append('Missing '+str(target.relative_to(ROOT)))
  checked+=1
  if 'run' in parse_qs(u.query):
   run_links+=1
   if parse_qs(u.query)['run'][0] not in ids:errors.append('Unknown run '+link)
records=json.loads((ROOT/'data/studies.json').read_text())['studies']
assert len(records)==1 and len(records[0]['analyses'])==14
assert records[0]['n_cells']==1377 and records[0]['whole_ovary_n_cells']==54748
assert len(json.loads((ROOT/'data/integrations.json').read_text())['integrations'])==0
hashes={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in ROOT.rglob('*') if p.is_file() and '.git' not in p.parts and 'docs' not in p.parts and '.hallmark' not in p.parts}
result={'local_links_checked':checked,'valid_run_links':run_links,'studies':len(records),'saved_runs':14,'integrations':0,'errors':errors,'sha256':hashes}
(ROOT/'docs/LIBRARY_VALIDATION.json').write_text(json.dumps(result,indent=2))
print(json.dumps({k:v for k,v in result.items() if k!='sha256'}))
assert not errors
