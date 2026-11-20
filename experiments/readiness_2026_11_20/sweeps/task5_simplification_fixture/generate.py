"""Generate the paired software fixture without invoking SLiM."""
from contextlib import redirect_stdout
import hashlib
import io
import json
from pathlib import Path

import msprime
import stdpopsim

HERE = Path(__file__).resolve().parent
SEED = 2026112091
ROOT = Path('/vast/projects/smathi/cohort/korfmann/ghist2026')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


demography = msprime.Demography()
demography.add_population(name='generic', initial_size=100)
for time, size in [(80, 160), (160, 120), (240, 120)]:
    demography.add_population_parameters_change(time=time, initial_size=size, population='generic')
model = stdpopsim.DemographicModel(
    id='BookkeepingFixture', description='Paired bookkeeping software fixture',
    long_description='Fixed small synthetic fixture; not biological recovery evidence.',
    generation_time=1, model=demography, mutation_rate=1e-8,
)
contig = stdpopsim.Contig.basic_contig(length=100_000, mutation_rate=1e-8, recombination_rate=1e-7)
contig.add_single_site(id='target', coordinate=50_000)
buffer = io.StringIO()
with redirect_stdout(buffer):
    stdpopsim.get_engine('slim').simulate(
        model, contig, samples={'generic': 12}, seed=SEED,
        extended_events=stdpopsim.selective_sweep(
            single_site_id='target', population='generic', mutation_generation_ago=60,
            selection_coeff=.2,
        ),
        slim_scaling_factor=1, slim_burn_in=1, slim_script=True,
    )
script = buffer.getvalue()
lines = script.splitlines()
tree_lines = [i for i, line in enumerate(lines) if 'defineConstant("trees_file",' in line]
assert len(tree_lines) == 1
lines[tree_lines[0]] = '    defineConstant("trees_file", "result.trees");'
script = '\n'.join(lines) + '\n'
script += '''
// Software diagnostics do not draw random numbers.
1 late() {
    sim.setValue("forced_restore_done", F);
    forced_tick = time_to_tick(drawn_mutations[0,0]) + 2;
    community.registerLateEvent(NULL,
        "{if (!sim.getValue('forced_restore_done')) {sim.setValue('forced_restore_done', T); catn('FIXTURE_FORCED_RESTORE'); restore();}}",
        forced_tick, forced_tick);
    community.registerLateEvent(NULL,
        "{mut=sim.mutations[sim.mutations.position==50000]; f=0.0; if (length(mut)>0) f=sum(sim.mutationFrequencies(p0, mut)); catn(paste(c('FIXTURE_TRAJECTORY', community.tick, p0.individualCount, f, sim.getValue('n_restores')), sep='\\t'));}",
        2, G0);
}
'''
needle = 'initializeTreeSeq(timeUnit="generations");'
assert script.count(needle) == 1
scripts = {
    'default.slim': script,
    'interval20.slim': script.replace(needle, 'initializeTreeSeq(timeUnit="generations", simplificationInterval=20);'),
}
for name, text in scripts.items():
    path = HERE / name
    if path.exists():
        raise FileExistsError(path)
    path.write_text(text)
engine = Path(stdpopsim.slim_engine.__file__)
protocol = {
    'schema': 'paired_slim_simplification_software_fixture_v1',
    'seed': SEED,
    'biology': {'length': 100000, 'recombination_rate': 1e-7, 'mutation_rate': 1e-8,
                'present_diploid_N': 100, 'backward_size_events': [[80, 160], [160, 120], [240, 120]],
                'burn_in': 1, 'scaling_factor': 1, 'diploid_samples': 12,
                'selected_position': 50000, 'selected_origin_generations': 60, 'selection_coefficient': .2,
                'conditioning': 'stdpopsim allele-survival conditioning; one additional forced checkpoint restore'},
    'arms': {'default': 'SLiM4.3 default ratio10 adaptive interval', 'interval20': 'explicit fixed simplificationInterval20'},
    'comparisons': ['logged tick/population/target-frequency/restore trajectory', 'forced restore exactly once per arm',
                    'remembered sample pedigree and genome identities', 'canonical marginal ancestry after sample simplification',
                    'binary selected-site genotypes'],
    'resource_limits': {'slurm_minutes': 5, 'cpus': 1, 'memory_MiB': 5632, 'per_arm_timeout_seconds': 100},
    'source_sha256': {name: digest(HERE / name) for name in ['generate.py', 'run.py', 'run.sbatch', *scripts]},
    'engine_path': str(engine), 'engine_sha256': digest(engine),
    'slim_path': str(ROOT / 'external/SLiM-v4.3/build/slim'),
    'slim_sha256': digest(ROOT / 'external/SLiM-v4.3/build/slim'),
    'versions': {'stdpopsim': stdpopsim.__version__, 'msprime': msprime.__version__},
    'scientifically_validated': False,
    'scope': 'Bookkeeping equivalence software fixture only; no full-chromosome resource or recovery claim.',
}
path = HERE / 'protocol.json'
if path.exists():
    raise FileExistsError(path)
path.write_text(json.dumps(protocol, indent=2) + '\n')
print(json.dumps({'protocol_sha256': digest(path), 'source_sha256': protocol['source_sha256']}))
