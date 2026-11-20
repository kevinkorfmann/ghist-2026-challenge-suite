#!/bin/bash
#SBATCH --job-name=ghist-lowcov-window-full
#SBATCH --partition=genoa-std-mem
#SBATCH --cpus-per-task=4
#SBATCH --mem-per-cpu=5632M
#SBATCH --time=01:00:00
#SBATCH --array=0-1%2
#SBATCH --output=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120/winsfs-full-%A_%a.log
source /etc/profile
set -euo pipefail
TASK_ROOT=/vast/projects/smathi/cohort/korfmann/ghist2026/experiments/lowcoverage-20261120
TASK_NAMES=(bottleneck split_isolation)
TASK_NAME=${TASK_NAMES[$SLURM_ARRAY_TASK_ID]}
cd "$TASK_ROOT/source_v7"
export PYTHONPATH="$TASK_ROOT/source_v7/src:$TASK_ROOT/source_v7"
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
verify_execution() {
 /vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python - "$TASK_ROOT" "$TASK_NAME" "$1" <<'PYVERIFY'
import json,sys,hashlib,os
from pathlib import Path
from ghist_lowcoverage.pipeline import audit_inputs
root,task,phase=Path(sys.argv[1]),sys.argv[2],sys.argv[3]
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest_path=root/'source_v7/source_manifest.json'
assert sha(manifest_path)=='43a608f4bde3334fc8239f48b4dd4834ee3a87ee547acb15b379bfc4838e21f4'
manifest=json.loads(manifest_path.read_text())
for name,digest in manifest.items():
    assert sha(root/'source_v7'/name)==digest,name
protocol=root/'winsfs_full_development_v1.json'
assert sha(protocol)=='6bb7811d87d0028dfa5d817886e5d26dbe171585e7939b4fc027e8e28d4415f6'
directory=root/f'full_{task}/remapped'
audit=audit_inputs(directory,directory/'reference.fasta',task,10_000_000)
binding={'source_manifest_sha256':sha(manifest_path),'source_hashes_verified':manifest,
         'protocol_sha256':sha(protocol),'launch_script_sha256':sha(root/'winsfs_full_v7.sbatch'),
         'input_audit':audit,
         'input_index_sha256':{str(index):sha(index) for index in [Path(str(directory/'reference.fasta')+'.fai'),*[Path(sample['path']+'.bai') for sample in audit['samples'].values()]]},
         'external_tool_sha256':{name:sha(root/name) for name in ['external/angsd/angsd','external/winsfs-release-v0.7.0/winsfs']},
         'environment_record_sha256':sha(root/'environment_v7.json'),
         'slurm_job_id':os.environ['SLURM_JOB_ID'],
         'slurm_array_task_id':os.environ['SLURM_ARRAY_TASK_ID']}
initial=root/f'full_execution_initial_{task}_v7.json'
if phase=='initial':
    assert not initial.exists()
    initial.write_text(json.dumps(binding,indent=2)+'\n')
else:
    assert json.loads(initial.read_text())==binding,'inputs, code, protocol, or launch wrapper changed'
    output=root/f'window_full_{task}_v7'
    assert json.loads((output/'safs/inputs.json').read_text())==audit,'actual SAF inputs differed'
    binding['source_and_input_pre_post_verification_passed']=True
    binding['initial_record_sha256']=sha(initial)
    binding['scored_result_sha256']=sha(output/'scored.json')
    binding['likelihood_protocol_sha256']=sha(output/'safs/protocol.json')
    binding['window_protocol_sha256']=sha(output/'window_comparison/protocol.json')
    (output/'execution_binding.json').write_text(json.dumps(binding,indent=2)+'\n')
PYVERIFY
}
verify_execution initial
/vast/projects/smathi/cohort/korfmann/ghist2026/.venv/bin/python scripts/lowcoverage_winsfs.py \
 --case "$TASK_ROOT/full_$TASK_NAME" --output "$TASK_ROOT/window_full_${TASK_NAME}_v7" \
 --angsd "$TASK_ROOT/external/angsd/angsd" --winsfs "$TASK_ROOT/external/winsfs-release-v0.7.0/winsfs" \
 --threads 4 --max-epochs 50

verify_execution final
