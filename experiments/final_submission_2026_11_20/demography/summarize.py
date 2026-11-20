"""Package locked best-effort point predictions; does not certify science or submit."""
import argparse
import json
import math
from pathlib import Path
import statistics
import sys

ROOT = Path(__file__).resolve().parent
sys.path[:0] = [str(ROOT / 'source_snapshot/src'), str(ROOT / 'source_snapshot'), str(ROOT)]
import yaml
from run import sha, verify_source, input_record
from ghist_inference.checkpoint import atomic_json
from ghist_inference.freeze import validate_physical_identity
from ghist_bottleneck.optimization import validate_optimizer_contract
from ghist_submission.validate import validate_submission, DEMOGRAPHY_FIELDS


def positive(value):
    return not isinstance(value, bool) and isinstance(value, (int, float)) and math.isfinite(value) and value > 0


def check_fit(record):
    case, fit = record['case'], record['fit']
    task = case['task']
    expected_module = ('src/ghist_bottleneck/sfs_baseline.py' if case.get('engine') == 'moments' else 'src/ghist_bottleneck/dadi_reference.py') if task == 1 else {2: 'src/ghist_secondary_contact/spectrum.py', 3: 'src/ghist_growth_bgs/spectrum.py', 4: 'src/ghist_split_migration_bgs/spectrum.py'}[task]
    imported = record.get('imported_fitter', {})
    manifest = json.loads((ROOT / 'source_manifest.json').read_text())
    if imported.get('relative_path') != expected_module or imported.get('sha256') != manifest[expected_module]:
        raise ValueError('Actual imported fitter identity differs from locked source')
    if fit.get('mutation_rate') != case['settings']['mutation_rate']:
        raise ValueError('Fit mutation rate differs from locked setting')
    for key in ('projection', 'folded', 'convention'):
        if key in case['settings'] and fit.get(key) != case['settings'][key]:
            raise ValueError('Fit setting differs from lock: ' + key)
    if not math.isfinite(fit['multinomial_log_likelihood']):
        raise ValueError('Nonfinite likelihood')
    if type(fit.get('best_restart_converged')) is not bool:
        raise ValueError('Selected convergence absent')
    if fit['sequence_length'] != record['fit_sequence_length_argument']:
        raise ValueError('Fit exposure differs from exact fitter argument')
    if task == 1:
        validate_optimizer_contract(fit, require_converged=False)
    elif task in (3, 4):
        validate_physical_identity(fit, task)
    else:
        theta, rate, length = fit['theta'], fit['mutation_rate'], fit['sequence_length']
        if not all(positive(x) for x in (theta, rate, length)):
            raise ValueError('Invalid scaling')
        n = theta / (4 * rate * length)
        scaled = fit['scaled_parameters']
        expected = {'ancestral_effective_size': n,
                    'population1_effective_size': n * scaled['population1_size_ratio'],
                    'population2_effective_size': n * scaled['population2_size_ratio'],
                    'migration_rate': scaled['migration_2N'] / (2*n)}
        if fit['convention'] == 'generator_matched':
            expected['generations_since_split'] = scaled['split_duration_2N'] * 2*n
            expected['generations_since_contact'] = expected['generations_since_split'] * scaled['contact_fraction']
        else:
            expected['generations_since_split'] = (scaled['isolation_duration_2N'] + scaled['contact_duration_2N']) * 2*n
            expected['generations_since_contact'] = scaled['contact_duration_2N'] * 2*n
        for key, value in expected.items():
            if not positive(fit[key]) or not math.isclose(fit[key], value, rel_tol=1e-12):
                raise ValueError('Scaled/physical arithmetic mismatch: '+key)
        rows = [r for r in fit['restart_diagnostics'] if r['restart'] == fit['best_restart']]
        if len(rows) != 1 or rows[0]['theta'] != theta or rows[0]['optimizer']['converged'] != fit['best_restart_converged']:
            raise ValueError('Selected theta/convergence differs from retained trace')
    return fit


def prediction(fit, task):
    if task == 2:
        rename = {'mainland_population_size': 'population1_effective_size',
                  'island_population_size': 'population2_effective_size',
                  'generations_since_migration': 'generations_since_contact'}
    else:
        rename = {}
    return {key: float(fit[rename.get(key,key)]) for key in DEMOGRAPHY_FIELDS[task]}


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--name', required=True);args=parser.parse_args()
    if not args.name.replace('_','').replace('-','').isalnum(): raise ValueError('Invalid unique package name')
    lock_hash=sha(ROOT/'method_lock.json');protocol=verify_source(lock_hash)
    output=ROOT/'packages'/args.name
    if output.exists(): raise ValueError('Refusing package overwrite')
    output.mkdir(parents=True)
    summaries={}
    for task in (1,2,3,4):
        eligible=[];states=[];all_cases=[c for c in protocol['fits'] if c['task']==task]
        try: inputs=input_record(task)
        except Exception as error:
            summaries[str(task)]={'state':'NO_VERIFIED_INPUT','reason':str(error)};continue
        spectrum_path=ROOT/'results'/f'task{task}'/'spectrum.json'
        if not spectrum_path.exists():
            summaries[str(task)]={'state':'EXTRACTION_PENDING','planned_fit_cases':len(all_cases)};continue
        spec=json.loads(spectrum_path.read_text())
        if spec['input_record'] != inputs or spec['method_lock_sha256'] != lock_hash:
            raise ValueError('Spectrum provenance mismatch')
        for case in all_cases:
            path=ROOT/'results'/f'task{task}'/(case['id']+'.json')
            state={'id':case['id'],'case':case,'state':'MISSING'}
            if path.exists():
                record=json.loads(path.read_text());state.update(state=record['state'],artifact=str(path.relative_to(ROOT)),sha256=sha(path))
                if record['case'] != case or record['method_lock_sha256'] != lock_hash or record['input_vcf_sha256'] != inputs['files']['vcf']['sha256'] or record['spectrum_sha256'] != spec['spectrum_sha256']:
                    raise ValueError('Fit input/method mismatch')
                if record['state']=='COMPLETED' and record.get('postfit_identity_verified') is True:
                    try:
                        fit=check_fit(record)
                        if fit['sequence_length'] != spec['sequence_length']:
                            raise ValueError('Fitted exposure differs from extraction')
                        values=prediction(fit,task)
                        # Validate each endpoint in the exact output contract before selection.
                        temporary=output/(case['id']+'.check.yaml');temporary.write_text(yaml.safe_dump(values,sort_keys=False))
                        try: validate_submission(temporary,task,'final')
                        finally: temporary.unlink()
                        eligible.append((case,fit,state))
                        state['eligible_for_best_effort']=True
                        state['selected_restart_converged']=fit['best_restart_converged']
                        state['multinomial_log_likelihood']=fit['multinomial_log_likelihood']
                    except Exception as error:
                        state['eligible_for_best_effort']=False;state['reason']=f'{type(error).__name__}: {error}'
                elif 'failure' in record:state['failure']=record['failure']
            states.append(state)
        selected=[];selection=None;complete_projection_set=None
        if task==1:
            for engine in ('moments','dadi'):
                selected=[r for r in eligible if r[0]['engine']==engine]
                if selected:selection=engine+' primary' if engine=='moments' else 'dadi deadline fallback';break
        elif task==2:
            eligible.sort(key=lambda r:(r[0]['settings']['convention']!='simple',r[0]['origin']))
            if eligible:
                selected=[max(eligible,key=lambda r:r[1]['multinomial_log_likelihood'])];selection='maximum finite likelihood on common unfolded40x30 SFS, equal5parameter models; theta optimized'
        else:
            projections=(84,60,40) if task==3 else ((50,40),(70,50),(94,68))
            by_projection={}
            for projection in projections:
                rows=[r for r in eligible if (r[0]['settings']['projection']==projection if task==3 else tuple(r[0]['settings']['projection'])==projection) and (task==4 or r[0]['settings']['convention']=='exponential')]
                if rows:
                    rows.sort(key=lambda r:r[0].get('origin',0));by_projection[projection]=max(rows,key=lambda r:r[1]['multinomial_log_likelihood'])
            complete_projection_set=len(by_projection)==3
            if complete_projection_set:
                selected=[by_projection[p] for p in projections];selection='componentwise physical projection median from all3fixed projections'
            elif by_projection:
                selected=[by_projection[next(p for p in projections if p in by_projection)]];selection='prespecified incomplete-panel projection fallback; not a full projection median'
        summary={'task':task,'planned_fit_cases':len(all_cases),'eligible_fit_cases':len(eligible),'case_states':states,'method_lock_sha256':lock_hash,'complete_projection_set':complete_projection_set,'scientifically_validated':False,'release_qualified':False,'submission_count':0,'intervals_generated':False}
        if selected:
            values=[prediction(r[1],task) for r in selected]
            final={k:float(statistics.median(v[k] for v in values)) for k in DEMOGRAPHY_FIELDS[task]}
            candidate=output/f'task{task}.final.yaml';candidate.write_text(yaml.safe_dump(final,sort_keys=False))
            format_check=validate_submission(candidate,task,'final')
            summary.update(state='BEST_EFFORT_CANDIDATE',selection=selection,prediction=final,
                           candidate=str(candidate),candidate_sha256=sha(candidate),format_check=format_check,
                           selected_cases=[r[0]['id'] for r in selected],selected_convergence=[r[1]['best_restart_converged'] for r in selected],
                           selected_fits=[{'case':r[0],'fit':r[1],'artifact':r[2]['artifact'],'sha256':r[2]['sha256']} for r in selected],
                           input_vcf_sha256=inputs['files']['vcf']['sha256'],input_record=inputs,spectrum_sha256=spec['spectrum_sha256'],sequence_length=spec['sequence_length'],
                           uncertainty='Point estimate only; scientific recovery, interval coverage, global optimum and misspecification robustness not established.')
        else:summary['state']='NO_ELIGIBLE_FINITE_FIT'
        summaries[str(task)]=summary
        atomic_json(output/f'task{task}.audit.json',summary)
    verify_source(lock_hash)
    atomic_json(output/'summary.json',{'method_lock_sha256':lock_hash,'packager_sha256':sha(Path(__file__)),'tasks':summaries,'no_external_submission':True})
    print(json.dumps({t:{k:v for k,v in row.items() if k in ('state','prediction','selected_cases','selected_convergence','candidate','candidate_sha256','planned_fit_cases','eligible_fit_cases')} for t,row in summaries.items()},indent=2))

if __name__=='__main__':main()
