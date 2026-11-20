"""Read-only integrity audit of existing development tree/output; no predictor runs or truth scoring."""
import hashlib,json,re,sys
from pathlib import Path
import numpy as np
import tskit
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
h=Path('/tmp/ghist-task7-development-retry-v2-terminal-handoff.json')
r=json.loads(h.read_text());s=Path(r['staging']);prefix=Path('experiments/readiness_2026_11_20/task7_development_replay');p=s/prefix/'retry_v2'
o=json.loads(Path(r['previous_terminal_handoff']).read_text());old=Path(o['staging'])/prefix
for hand in [r,o]:
 assert sha(Path(hand['archive']))==hand['archive_sha256']
 for name,rec in hand['owned_manifest'].items():
  f=Path(hand['staging'])/name
  assert sha(f)==rec['sha256'] and f.stat().st_size==rec['bytes']
protocol=json.loads((p/'protocol.json').read_text());previous=json.loads((old/'protocol.json').read_text())
for key in ['source_sha256','input_sha256','input_tree','software','method','settings','seed','original_generator']:
 assert protocol[key]==previous[key],key
completion=json.loads((p/'output/completion.json').read_text())
assert completion['state']=='COMPLETED'
assert completion['protocol_sha256']==sha(p/'protocol.json')
assert completion['software']==protocol['software']
for name,digest in completion['output_sha256'].items():assert sha(p/'output'/name)==digest
failure=json.loads((old/'output/completion.json').read_text())
assert failure['state']=='FAILED' and 'reference sequence' in failure['failure']
assert failure['slurm_job_id']=='10858405'
assert not (old/'output/canonical_input.vcf').exists()
for name,digest in failure['output_sha256'].items():assert sha(old/'output'/name)==digest
external=Path(r['external_tree']['local_preservation_path'])
assert sha(external)==r['external_tree']['sha256'] and external.stat().st_size==r['external_tree']['bytes']
ts=tskit.load(external); tables=ts.dump_tables();cleared=ts.dump_tables();cleared.reference_sequence.clear()
for name in ['individuals','nodes','edges','migrations','sites','mutations','populations','provenances']:
 assert getattr(tables,name).equals(getattr(cleared,name)),name
assert hashlib.sha256(ts.reference_sequence.data.encode()).hexdigest()==protocol['reference_handling']['expected_data_sha256']
assert len(ts.reference_sequence.data)==protocol['reference_handling']['expected_data_length']==ts.sequence_length
left=protocol['coordinate_transform']['original_left'];right=protocol['coordinate_transform']['original_right']
local=cleared.tree_sequence().keep_intervals([[left,right]],simplify=False).trim()
assert local.sequence_length==right-left==1000000
assert np.array_equal(local.samples(),ts.samples())
assert np.array_equal(local.genotype_matrix(),ts.genotype_matrix())
assert np.array_equal(local.tables.sites.position,ts.tables.sites.position-left)
for name in ['individuals','nodes','migrations','mutations','populations']:
 assert getattr(local.tables,name).equals(getattr(tables,name)),name
for name in ['parent','child','metadata','metadata_offset']:
 assert np.array_equal(getattr(local.tables.edges,name),getattr(tables.edges,name)),name
for name in ['left','right']:
 assert np.array_equal(getattr(local.tables.edges,name),getattr(tables.edges,name)-left),name
for name in ['ancestral_state','ancestral_state_offset','metadata','metadata_offset']:
 assert np.array_equal(getattr(local.tables.sites,name),getattr(tables.sites,name)),name
# Parse retained VCF directly, independently from production reader/exporter.
vcf=(p/'output/canonical_input.vcf').read_text().splitlines()
assert '##contig=<ID=7,length=1000000>' in vcf
header=next(x for x in vcf if x.startswith('#CHROM')).split('\t')
assert len(header[9:])==25
records=[line.split('\t') for line in vcf if not line.startswith('#')]
nodes=list(local.samples());groups={}
for n in nodes:groups.setdefault(local.node(n).individual,[]).append(n)
assert all(i>=0 and len(group)==2 for i,group in groups.items())
nodes=[n for group in groups.values() for n in group]
accepted=[];skipped=0
for variant in local.variants(samples=nodes):
 calls=variant.genotypes
 if len(variant.alleles)==2 and np.all((calls==0)|(calls==1)) and 0<int(calls.sum())<len(calls):
  accepted.append((int(variant.site.position),calls))
 else:skipped+=1
assert len(records)==len(accepted)==777 and skipped==1
for fields,(position,calls) in zip(records,accepted,strict=True):
 assert fields[0]=='7' and int(fields[1])==position+1
 assert fields[3:5]==['A','C'] and fields[8]=='GT'
 assert fields[9:]==[str(calls[i])+'|'+str(calls[i+1]) for i in range(0,len(calls),2)]
masked=[line for line in (p/'output/prediction/masked.vcf').read_text().splitlines() if not line.startswith('#')]
assert masked==['\t'.join(fields) for fields in records]
assert (p/'inputs/accessible.bed').read_text()=='7\t0\t1000000\n'
report=np.loadtxt(p/'output/prediction/RAiSD_Report.task7',comments='/',ndmin=2)
assert report.shape[0]==728 and np.isfinite(report).all()
assert np.all(report[:,0]>=1) and np.all(report[:,0]<=1000000)
assert re.findall(r'\|\s+snps\s+(\d+)',(p/'output/prediction/RAiSD_Info.task7').read_text())==['777']
pred=json.loads((p/'output/prediction/prediction.json').read_text())
assert pred['settings']==protocol['settings'] and pred['method']==protocol['method']
assert pred['consensus_eligible'] is False and pred['signal_ratio'] is None and pred['distant_peak_score'] is None
assert 1000000<protocol['settings']['secondary_peak_minimum_distance']
assert pred['input_provenance']['vcf_sha256']==sha(p/'output/canonical_input.vcf')
assert pred['input_provenance']['accessibility_sha256']==sha(p/'inputs/accessible.bed')
for name,digest in pred['source_sha256'].items():assert sha(p/'source_snapshot'/name)==digest
for name,digest in pred['artifact_sha256'].items():assert sha(p/'output/prediction'/name)==digest
lo,hi=pred['interval'];assert 0<=lo<hi<=1000000
assert (p/'output/prediction/candidate.bed').read_text()==f'7\t{lo}\t{hi}\n'
assert sha(external)==r['external_tree']['sha256']
print(json.dumps({'status':'passed','local_audit_versions':{'python':sys.version,'numpy':np.__version__,'tskit':tskit.__version__},'frozen_betty_runtime':protocol['software'], 'retry_owned_files_verified':len(r['owned_manifest']),'original_owned_files_verified':len(o['owned_manifest']),'source_files_identical_across_attempts':len(protocol['source_sha256']),'output_hashes_verified':len(completion['output_sha256']),'tree_sha256':sha(external),'original_failure_retained':True,'reference_clear_tables_unchanged':True,'exact_edge_and_site_translation':True,'sample_mutation_and_genotype_invariance':True,'direct_vcf_roundtrip':{'accepted_sites':777,'skipped_sites':1,'diploid_samples':25,'length':1000000},'masked_vcf_record_parity':True,'raisd_windows':len(report),'raisd_usable_snps':777,'expected_short_segment_abstention':True,'retained_bed_coordinate_parity':True,'predictor_reexecuted':False,'raisd_reexecuted':False,'accuracy_evaluations':0,'scope':'Existing development evidence audit; no new inference, models, biology, truth comparisons or qualification.'},indent=2,allow_nan=False))
