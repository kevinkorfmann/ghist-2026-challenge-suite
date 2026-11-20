"""Read-only resource/hash audit; no simulation, prediction, outcome scoring or reexecution."""
from pathlib import Path
import csv,datetime,hashlib,importlib.util,io,json,math,re,subprocess,tarfile
base=Path('/Users/kevinkorfmann/.cache/ghist-2026/complete-integration-20261120/experiments/readiness_2026_11_20/sweeps')
p=base/'task5_interval100_feasibility';old=base/'task5_interval20_feasibility';t=p/'terminal'
sha=lambda path:hashlib.sha256(path.read_bytes()).hexdigest()
protocol=json.loads((p/'protocol.json').read_text());c=json.loads((t/'terminal_collection.json').read_text())
assert sha(p/'protocol.json')=='9004cf2349334767c1bb9777d3ad9ef7063b961876e8dbd92ac7222dc1901307'
assert sha(t/'terminal_collection.json')=='8e9142561bd8e88622ae2290e2a9e96811c2c4a75d8939f0431ca8b6ec0b5099'
assert sha(p/'collect_terminal.py')==c['collector_sha256']
assert not c['post_run_mismatches']
for name,digest in protocol['source_sha256'].items():assert sha(p/name)==digest,name
for name,rec in c['output_files'].items():
 f=t/name;assert sha(f)==rec['sha256'] and f.stat().st_size==rec['bytes'],name
with tarfile.open(t/'generated_scripts.tar.gz') as ar:
 for member in ar.getmembers():
  if member.isfile():
   leaf=Path(member.name).name
   assert leaf in ('generated_original.slim','generated_interval100.slim')
   assert ar.extractfile(member).read()==(t/'diagnostics'/leaf).read_bytes()
resources=(t/'diagnostics/resources.jsonl').read_bytes();phase=(t/'diagnostics/phase.tsv').read_bytes();mid=json.loads((p/'midpoint_observation.json').read_text())
for name,raw in [('resources',resources),('phase',phase)]:assert hashlib.sha256(raw[:mid[name+'_prefix_bytes']]).hexdigest()==mid[name+'_prefix_sha256']
rows=[json.loads(line) for line in resources.splitlines()];times=[r['unix_time'] for r in rows];gaps=[b-a for a,b in zip(times,times[1:])]
assert all(math.isfinite(v) for v in times) and min(gaps)>0 and len(rows)==271
jobgroups=[next(v for k,v in r['cgroups'].items() if k.endswith('/job_10857763')) for r in rows]
limits={int(v['memory.limit_in_bytes']) for v in jobgroups};assert limits=={22*1024**3}
for r in rows:
 for v in r['cgroups'].values():
  assert int(v['memory.failcnt'])==0
  oom=dict(line.split() for line in v['memory.oom_control'].splitlines())
  assert int(oom['under_oom'])==int(oom['oom_kill'])==0
assert not (t/'diagnostics/monitor_failure.json').exists()
hi=max(int(v['memory.max_usage_in_bytes']) for v in jobgroups)
usage=max(int(v['memory.usage_in_bytes']) for v in jobgroups)
assert hi==c['resource_trace']['maximum_recorded_high_water_bytes']==8365776896
assert usage==c['resource_trace']['maximum_recorded_usage_bytes']==8287920128
assert max(gaps)==c['resource_trace']['maximum_sample_gap_seconds']
assert {r['stage'] for r in rows}=={'initializing','slim_running'}
phases=list(csv.reader(io.StringIO(phase.decode()),delimiter='\t'))
assert phases.pop(0)==['tick','diploid_population_size','saves','restores']
phases=[list(map(int,row)) for row in phases];ticks=[r[0] for r in phases]
assert ticks==[1]+list(range(100,33601,100))
assert all(r[2:]==[0,0] for r in phases)
assert phases[-1]==[33600,52855,0,0]
assert max(r[1] for r in phases)==71194
scheduler=list(csv.DictReader(io.StringIO((t/'sacct_utc.txt').read_text()),delimiter='|'));job=scheduler[0];batch=scheduler[1]
assert job['JobID']=='10857763' and job['State']=='TIMEOUT' and job['Elapsed']=='00:45:10'
assert job['Timelimit']=='00:45:00' and job['AllocCPUS']=='4' and job['ReqMem']=='22G'
assert batch['State']=='CANCELLED' and batch['ExitCode']=='0:15'
parse=lambda v:datetime.datetime.fromisoformat(v).replace(tzinfo=datetime.timezone.utc).timestamp()
start_gap=times[0]-parse(job['Start']);tail_job=parse(job['End'])-times[-1];tail_batch=parse(batch['End'])-times[-1]
assert 0<=start_gap<10 and 0<tail_job<3
assert json.loads((t/'diagnostics/adapter_status.json').read_text())['state']=='RUNNING'
assert json.loads((t/'case_005/status.json').read_text())['state']=='RUNNING'
# Compare complete source strings, normalizing only explicitly registered filesystem destinations/interval.
pattern=re.compile(r'^    defineConstant\("trees_file", "([^"\n]+)"\);$',re.MULTILINE)
a=(old/'terminal/diagnostics/generated_original.slim').read_text();b=(t/'diagnostics/generated_original.slim').read_text()
assert len(pattern.findall(a))==len(pattern.findall(b))==1
normal=lambda s:pattern.sub('    defineConstant("trees_file", "NORMALIZED_TEMP_TREE");',s)
assert normal(a)==normal(b)
fa=(old/'terminal/diagnostics/generated_interval20.slim').read_text();fb=(t/'diagnostics/generated_interval100.slim').read_text()
oldout=json.loads((old/'protocol.json').read_text())['output_directory'];newout=protocol['output_directory']
assert fa.replace('simplificationInterval=20,','simplificationInterval=100,').replace(oldout+'/diagnostics/raw_slim.trees',newout+'/diagnostics/raw_slim.trees').replace(oldout+'/diagnostics/phase.tsv',newout+'/diagnostics/phase.tsv')==fb
spec=importlib.util.spec_from_file_location('audit_adapter',p/'adapter.py');adapter=importlib.util.module_from_spec(spec);spec.loader.exec_module(adapter)
assert adapter.transform_script(b,Path(newout)/'diagnostics/raw_slim.trees',Path(newout)/'diagnostics/phase.tsv')==fb
assert protocol['seed']==202612005 and protocol['case']==5
assert not (t/'diagnostics/raw_slim.trees').exists()
# Re-attest remote source/runtime/locks independently of the collector; only byte hashes are returned.
remote_code='''from pathlib import Path
import datetime,hashlib,json,importlib.metadata,platform
p=Path('/vast/projects/smathi/cohort/korfmann/ghist2026/reviews/readiness-sweeps-interval100-20261120')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
r=json.loads((p/'protocol.json').read_text())
assert sha(p/'protocol.json')=='9004cf2349334767c1bb9777d3ad9ef7063b961876e8dbd92ac7222dc1901307'
for name,digest in r['source_sha256'].items():assert sha(p/name)==digest,name
for name,digest in r['external_sha256'].items():assert sha(Path(name))==digest,name
versions={name:importlib.metadata.version(name) for name in r['versions']}
assert versions==r['versions'] and platform.python_version()==r['python']['version']
print(json.dumps({'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'source_hashes':len(r['source_sha256']),'external_hashes':len(r['external_sha256']),'versions':versions,'python':platform.python_version(),'all_match':True}))
'''
remote=subprocess.run(['ssh','-o','BatchMode=yes','betty','source /etc/profile && module load python/3.11.11 && /vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python -'],input=remote_code,text=True,capture_output=True,check=True)
payloads=[line for line in remote.stdout.splitlines() if line.startswith('{')]
assert len(payloads)==1
attest=json.loads(payloads[0])
print(json.dumps({'status':'passed','scope':'Read-only resource/source audit, no accuracy inspection or run.','hashes':{name:sha(p/name) for name in ['protocol.json','collect_terminal.py','terminal/terminal_collection.json','terminal/diagnostics/resources.jsonl','terminal/diagnostics/phase.tsv','terminal/sacct_utc.txt','terminal/generated_scripts.tar.gz','generated_script_equivalence.json']},'local_source_hashes_verified':len(protocol['source_sha256']),'terminal_files_verified':len(c['output_files']),'remote_independent_attestation':attest,'monitoring':{'records':len(rows),'span_seconds':times[-1]-times[0],'max_gap_seconds':max(gaps),'startup_to_first_sample_seconds':start_gap,'last_sample_to_job_end_seconds':tail_job,'last_sample_to_batch_end_seconds':tail_batch,'all_sampled_cgroup_failcnt_and_oom_zero':True,'observed_high_water_bytes':hi,'observed_high_water_GiB':hi/1024**3,'maximum_sampled_usage_bytes':usage,'job_limit_GiB':22,'observed_headroom_GiB':22-hi/1024**3,'midpoint_prefixes_preserved':True,'monitor_failure_file_absent':True,'unobserved_tail_and_future_memory_not_certified':True},'phases':{'rows':len(phases),'all_ticks_1_then_every100_to33600':True,'last':phases[-1],'maximum_logged_population':max(r[1] for r in phases),'source_end_tick':c['source_schedule']['end_tick'],'saves_and_restores_observed':0,'only_initializing_and_slim_running_stages_observed':True,'final_raw_tree_absent':True,'stale_running_status_is_interrupted_cleanup':True},'same_seed_source_comparison':{'original_scripts_identical_after_one_temporary_loader_path_normalization':True,'instrumented_scripts_identical_after_only_registered_interval_and_diagnostic_paths':True,'exact_adapter_transform_verified':True,'biology_and_event_statements_unchanged':True},'scheduler':{'job_state':job['State'],'job_elapsed':job['Elapsed'],'batch_state':batch['State'],'batch_exit_code':batch['ExitCode'],'batch_max_rss':batch['MaxRSS']},'scientifically_validated':False},indent=2,allow_nan=False))
