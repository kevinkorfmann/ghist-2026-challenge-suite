from pathlib import Path
import csv,gzip,hashlib,itertools,json,sys,datetime
base=Path('/Users/kevinkorfmann/.cache/ghist-2026/final-task12-20261120-v2')
snapshot=Path('/Users/kevinkorfmann/.cache/ghist-2026/secondary-contact-interval-v2-integration-20261120')
sys.path.insert(0,str(snapshot/'src'))
from ghist_submission.validate import validate_submission
from ghist_relatedness.labels import HUMAN_LABELS
central=Path('/private/tmp/ghist-final-submissions-20261120/task12')
source=base/'output/task12.final.csv';status=json.loads((base/'output/status.json').read_text())
assert status['status']=='COMPLETE'
h=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
assert h(source)==status['output_sha256']=='b94c2bef386ffcf062ad2fa1fe0c0593b2d2dadfc49efdaa916596897f718dac'
assert h(central/'relatedness4.final.vcf.gz')==status['input']['sha256']
with gzip.open(central/'relatedness4.final.vcf.gz','rt') as handle:
 for line in handle:
  if line.startswith('#CHROM\t'):
   samples=line.rstrip('\r\n').split('\t')[9:];break
expected=set(itertools.combinations(sorted(samples),2))
assert len(samples)==297 and len(expected)==43956
with source.open() as handle:
 reader=csv.DictReader(handle);assert reader.fieldnames==['indID1','indID2','relation'];rows=list(reader)
assert len(rows)==len(expected)
actual=[(r['indID1'],r['indID2']) for r in rows]
assert len(set(actual))==len(actual) and set(actual)==expected
assert set(r['relation'] for r in rows)<=set(HUMAN_LABELS.values())
validation=validate_submission(source,12,'final',expected)
destination=central/'task12.final.csv'
with destination.open('xb') as handle:handle.write(source.read_bytes())
assert h(destination)==h(source)
review={'created_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'Independent local CSV/VCF sample header and byte identity checks; no scientific accuracy claim or external submission.','candidate_sha256':h(destination),'input_sha256':status['input']['sha256'],'sample_count':len(samples),'expected_pair_count':len(expected),'actual_pair_count':len(rows),'unique_pair_count':len(set(actual)),'label_domain_verified':True,'schema':validation,'checks_passed':True}
(base/'output/independent_local_validation.json').write_text(json.dumps(review,indent=2)+'\n')
manifest={'task':12,'dataset':'final','ready_for_central_uploader':True,'candidate_path':str(destination),'candidate_sha256':h(destination),'candidate_bytes':destination.stat().st_size,'input_entity':'syn75203676','input_version':1,'input_sha256':status['input']['sha256'],'input_md5':status['input']['md5'],'input_provenance_sha256':h(central/'input_provenance.json'),'method_lock_sha256':status['lock_sha256'],'source_manifest_sha256':status['source_manifest_sha256'],'validation':review,'slurm_job_id':'10859506','slurm_status':'COMPLETED exit0 in2m28s','all_pairs':43956,'sample_count':297,'zero_segment_pairs':0,'segment_present_pairs':43956,'model':'fixed historical empirical_background_v1','correction':'legacy task12-top20-v1; quantile0.25;minimum5cM;corrected weight1','preparation':'Major-allele completion of148382missing alleles;3multiallelic sites skipped;249062sites retained. Original observed alleles and phase guards preserved.','first_attempt_preserved':'10859426 strict missing-allele guard failed before detector/model','scientifically_validated':False,'scientific_limitations':status['scientific_limitations'],'external_submission_performed':False,'uploader_owner':'Other user chat central uploader only; root coordinates.'}
(central/'candidate_manifest.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
(base/'handoff_complete.json').write_text(json.dumps(manifest,indent=2,allow_nan=False)+'\n')
print(json.dumps({'candidate':str(destination),'sha256':h(destination),'bytes':destination.stat().st_size,'pairs':len(rows),'manifest':str(central/'candidate_manifest.json'),'manifest_sha256':h(central/'candidate_manifest.json')},indent=2))
