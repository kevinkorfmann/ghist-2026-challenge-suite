from pathlib import Path
import contextlib,importlib.util,io,json,sys,tempfile
import pytest
root=Path(sys.argv[1]);sys.path[:0]=[str(root/'src'),str(root)]
spec=importlib.util.spec_from_file_location('fixture',root/'tests/test_demography_freeze_cross_review.py');fixture=importlib.util.module_from_spec(spec);spec.loader.exec_module(fixture)
keys=('mainland_population_size','island_population_size','generations_since_split','generations_since_migration','migration_rate')
base={key:[.5,3.] for key in keys}
cases={'complete':base,'empty':{},'reversed':dict(base,migration_rate=[2.,1.]),'negative':dict(base,mainland_population_size=[-3.,-1.])}
out={}
for name,intervals in cases.items():
 with tempfile.TemporaryDirectory() as d,pytest.MonkeyPatch.context() as m,contextlib.redirect_stdout(io.StringIO()):
  src=fixture.secondary_contact_sources();src['bootstrap']['observed_sampling_interval_95']=intervals
  try:
   result=fixture.run_secondary_contact_freeze(src,Path(d),m);out[name]={'accepted':True,'exported':result['conservative_sampling_interval_95']}
  except Exception as error:out[name]={'accepted':False,'exception':type(error).__name__,'message':str(error)}
print(json.dumps(out,indent=2))
