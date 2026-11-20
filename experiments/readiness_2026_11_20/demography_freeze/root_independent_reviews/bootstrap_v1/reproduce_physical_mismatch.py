"""Read-only arithmetic integrity reproduction; synthetic fits, no model evaluation."""
from pathlib import Path
import sys,json,tempfile,io,contextlib
root=Path(sys.argv[1]).resolve()
sys.path[:0]=[str(root/'src'),str(root),str(root/'tests')]
import pytest
from test_demography_bootstrap_contract import bootstrap_fixture
from test_demography_freeze_integrity import sources_for,freeze
from ghist_inference.bootstrap import summarize_payloads,require_summary
results=[]
for task in [3,4]:
 folder=Path(tempfile.mkdtemp(prefix='bootstrap-physical-main-repro-'))
 payload=bootstrap_fixture(folder,task)
 for record in payload['results']:
  fit=record['fit'];fit['sequence_length']=100000.;fit['theta']=100.;fit['restart_diagnostics'][0]['theta']=100.
  n=fit['theta']/(4*fit['mutation_rate']*fit['sequence_length'])
  if task==3:
   fit.update(ancestral_effective_size=n,modern_effective_size=n*fit['post_growth_scaling'],scaled_duration_2N=fit['generations']/(2*n))
  else:
   fit.update(ancestral_population_size=n,scaled_parameters={'AFR_size_ratio':fit['AFR_population_size']/n,'EUR_size_ratio':fit['EUR_population_size']/n,'split_duration_2N':fit['generations_since_split']/(2*n),'migration_2N':fit['migration_rate']*2*n})
 parameter='generations' if task==3 else 'generations_since_split'
 original=require_summary(summarize_payloads([payload],task),task)[parameter]
 payload['results'][0]['fit'][parameter]*=10
 summary=summarize_payloads([payload],task)
 module,sources,flags=sources_for(task,folder);sources['bootstrap']=summary
 with pytest.MonkeyPatch.context() as mp,contextlib.redirect_stdout(io.StringIO()):
  output=freeze(module,sources,flags,folder,mp)
 results.append({'task':task,'original_interval':original,'changed_field':parameter,'change':'replica0 times10; selected theta and scaled values held fixed','summary_state':summary['state'],'actual_freeze_interval':output['linked_block_interval_95'][parameter]})
print(json.dumps(results,indent=2))
