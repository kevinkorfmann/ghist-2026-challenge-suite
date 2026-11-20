import hashlib,json,os,pathlib,platform,subprocess,sys,time,xml.etree.ElementTree as ET
root=pathlib.Path(sys.argv[1]).resolve(); out=root/sys.argv[2]; out.mkdir(exist_ok=False)
env=dict(os.environ);env['PYTHONPATH']=str(root/'src')+os.pathsep+str(root)+(os.pathsep+env['PYTHONPATH'] if env.get('PYTHONPATH') else '')
env['OMP_NUM_THREADS']='1';env['OPENBLAS_NUM_THREADS']='1'
steps=[('source_before',[sys.executable,'collect_evidence.py',str(root),str(root/'source_snapshot.json')]),('fixtures_before',[sys.executable,'verify_fixtures.py']),('ruff',[sys.executable,'-m','ruff','check','src','scripts','tests']),('pytest',[sys.executable,'-m','pytest','-q','--junitxml='+str(out/'pytest.xml')]),('source_after',[sys.executable,'collect_evidence.py',str(root),str(root/'source_snapshot.json')]),('fixtures_after',[sys.executable,'verify_fixtures.py'])]
records={}
for name,cmd in steps:
 start=time.monotonic()
 with (out/(name+'.log')).open('wb') as f:
  try:
   r=subprocess.run(cmd,cwd=root,env=env,stdout=f,stderr=subprocess.STDOUT,timeout=600);code=r.returncode
  except subprocess.TimeoutExpired: code=124
 log=out/(name+'.log');records[name]={'exit_code':code,'elapsed_seconds':round(time.monotonic()-start,3),'log_sha256':hashlib.sha256(log.read_bytes()).hexdigest()}
 if code and name!='pytest': break
summary={'platform':platform.platform(),'python':sys.version,'source_manifest_sha256':hashlib.sha256((root/'source_snapshot.json').read_bytes()).hexdigest(),'fixture_manifest_sha256':hashlib.sha256((root/'fixture_snapshot.json').read_bytes()).hexdigest(),'steps':records,'all_checks_passed':len(records)==len(steps) and all(x['exit_code']==0 for x in records.values())}
if (out/'pytest.xml').exists():
 suites=ET.parse(out/'pytest.xml').getroot();summary['pytest_suites']=[x.attrib for x in suites.iter('testsuite')]
(out/'summary.json').write_text(json.dumps(summary,indent=2)+'\n');print(json.dumps(summary,indent=2),flush=True)
if not summary['all_checks_passed']:
 for name,r in records.items():
  if r['exit_code']: print(name,(out/(name+'.log')).read_text(errors='replace')[-6000:])
 sys.exit(1)
