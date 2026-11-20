from pathlib import Path
import hashlib,json,datetime
root=Path('/private/tmp/ghist-final-submissions-20261120');rows=[]
for n in [1,2,3,4,5,6,7,8,12]:
 p=root/f'task{n}'/'input_provenance.json';x=json.loads(p.read_text());f=Path(x['path']);h=hashlib.sha256();m=hashlib.md5()
 with f.open('rb') as handle:
  for b in iter(lambda:handle.read(1024*1024),b''):h.update(b);m.update(b)
 checks={'bytes':f.stat().st_size==x['bytes'],'sha256':h.hexdigest()==x['sha256'],'md5':m.hexdigest()==x['md5'],'live_metadata_md5':m.hexdigest()==x['metadata_md5'],'current_version':x['current_version_confirmed'] is True,'version_integer':isinstance(x['version'],int) and x['version']>0}
 assert all(checks.values()),(n,checks)
 rows.append({'task':n,'entity':x['entity'],'version':x['version'],'filename':x['filename'],'sha256':h.hexdigest(),'md5':m.hexdigest(),'bytes':f.stat().st_size,'provenance_sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'checks':checks,'status':'PASS'})
report={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Independent local byte identity versus root browser-collected current-version provenance for Tasks1–8 and12. Browser metadata attestation belongs to root; no new browser authentication, final inference or competition submission in this audit. Auxiliary references are separately bound by their owners.','inputs':rows,'all_inputs_pass':True}
out=root/'independent_input_identity_audit_v2.json';out.write_text(json.dumps(report,indent=2)+'\n');print(json.dumps({'path':str(out),'sha256':hashlib.sha256(out.read_bytes()).hexdigest(),'all_pass':True,'inputs':[(x['task'],x['entity'],x['version']) for x in rows]},indent=2))
