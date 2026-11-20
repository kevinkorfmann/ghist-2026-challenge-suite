"""Check future pedigree definitions against Ped-sim parentage without genetics."""
import argparse
import json
import os
import subprocess
from pathlib import Path

from nested_qualification_design import compile_definition

from scripts.develop_task10_detected_ibd import family_truth, parents_from_fam, save_json, sha256


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-root', type=Path, required=True)
    parser.add_argument('--protocol', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if not os.environ.get('SLURM_JOB_ID'):
        parser.error('Pedigree structural checks require Slurm')
    protocol = json.loads(args.protocol.read_text())
    args.output.mkdir(parents=True, exist_ok=False)
    entries = protocol['cohorts']
    definition = args.output / 'nested.def'
    definition.write_text(''.join(row['definition'] for row in entries))
    assets = args.source_root / 'assets'
    with (args.output / 'dry_run.log').open('w') as log:
        subprocess.run(list(map(str, [assets / 'external/ped-sim/ped-sim', '-d', definition,
            '-m', assets / 'external/ped-sim-map/refined_mf.simmap', '-o', args.output / 'nested',
            '--intf', assets / 'external/ped-sim/interfere/nu_p_campbell.tsv',
            '--dry_run', '--seed', 2026200010])), stdout=log, stderr=subprocess.STDOUT, check=True)
    parents = parents_from_fam(args.output / 'nested-everyone.fam')
    verified = []
    for row in entries:
        _, roster, expected = compile_definition(row['terminal_paths'], row['name'], row['people'])
        actual = family_truth(parents, {sample: sample for sample in roster})
        if {(item['indID1'], item['indID2']): item['relation'] for item in actual} != expected:
            raise ValueError(f"Parentage differs from declared nested design {row['name']}")
        truth_path = args.output / f"{row['name']}_truth.json"
        save_json(truth_path, {'cohort_id': row['name'], 'roster': roster, 'truth': actual})
        verified.append({'name': row['name'], 'people': len(roster), 'pairs': len(actual),
            'labels': len({item['relation'] for item in actual}), 'all_declared_pairs_verified': True,
            'truth_path': truth_path.name, 'truth_sha256': sha256(truth_path)})
    result = {'state': 'STRUCTURES_VERIFIED_NO_GENETIC_DATA_GENERATED', 'cohorts': verified,
        'protocol_sha256': sha256(args.protocol), 'pedigree_sha256': sha256(args.output / 'nested-everyone.fam'),
        'definition_sha256': sha256(definition), 'validation_source_sha256': sha256(Path(__file__))}
    save_json(args.output / 'result.json', result)
    print(json.dumps(result))


if __name__ == '__main__':
    main()
