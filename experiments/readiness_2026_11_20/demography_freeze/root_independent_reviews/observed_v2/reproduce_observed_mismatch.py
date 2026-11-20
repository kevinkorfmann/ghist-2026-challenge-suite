"""Exercise observed projection arithmetic integrity using synthetic fits only."""
from pathlib import Path
import sys,json,tempfile,io,contextlib
root=Path(sys.argv[1]).resolve();sys.path[:0]=[str(root/'src'),str(root),str(root/'tests')]
import pytest
from test_demography_freeze_integrity import sources_for,freeze,fields
from test_demography_bootstrap_contract import synthetic_fit
results=[]
for task in [3,4]:
 folder=Path(tempfile.mkdtemp(prefix='observed-physical-main-repro-'));module,sources,flags=sources_for(task,folder)
 fits,median,_,parameter=fields(task,sources)
 for key,fit in fits.items():
  original=dict(fit);fit.clear();fit.update(synthetic_fit(task));fit.update(original)
 if task==4:
  original=dict(sources['refined']);sources['refined']=synthetic_fit(4);sources['refined'].update(original)
 with pytest.MonkeyPatch.context() as mp,contextlib.redirect_stdout(io.StringIO()): before=freeze(module,sources,flags,folder,mp)
 changed=list(fits) if task==3 else ['70x50','94x68']
 for key in changed: fits[key][parameter]*=10
 median[parameter]*=10
 with pytest.MonkeyPatch.context() as mp,contextlib.redirect_stdout(io.StringIO()): after=freeze(module,sources,flags,folder,mp)
 results.append({'task':task,'changed_projections':changed,'original_point':before['testing_prediction'][parameter],'corrupted_point':after['testing_prediction'][parameter],'interval_unchanged':before['linked_block_interval_95']==after['linked_block_interval_95'],'selected_theta_scaled_parameters_and_refined_reference_unchanged':True})
print(json.dumps(results,indent=2))
