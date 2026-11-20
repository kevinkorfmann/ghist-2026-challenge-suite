"""Isolated, one-case ancestry-bookkeeping feasibility adapter."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
import os
from pathlib import Path
import re
import runpy
import shutil
import sys
import threading
import time

HERE = Path(__file__).resolve().parent
DECLARATION = 'initializeTreeSeq(timeUnit="generations");'
FIXED_DECLARATION = 'initializeTreeSeq(simplificationInterval=100, timeUnit="generations");'
TREES = re.compile(r'^    defineConstant\("trees_file", "([^"\n]+)"\);$', re.MULTILINE)


def sha(path):
    with Path(path).open('rb') as handle:
        return hashlib.file_digest(handle, 'sha256').hexdigest()


def transform_script(script, checkpoint, phase_log):
    """Change one bookkeeping control and output locations; add no RNG calls."""
    if script.count(DECLARATION) != 1 or len(TREES.findall(script)) != 1:
        raise ValueError('Expected exactly one tree-recording declaration and output path')
    for path in (checkpoint, phase_log):
        if any(char in str(path) for char in ('"', '\\', '\n', '\r')):
            raise ValueError('Diagnostic output path cannot be represented safely in Eidos')
    changed = script.replace(DECLARATION, FIXED_DECLARATION)
    changed = TREES.sub(lambda _: f'    defineConstant("trees_file", "{checkpoint}");', changed)
    diagnostics = f'''
// Passive phase logging: no random draws or population/allele modifications.
1: late() {{
    if ((community.tick % 100 == 0) | (community.tick == 1) | (community.tick == G0)) {{
        writeFile("{phase_log}", paste(c(community.tick,
            sum(sim.subpopulations.individualCount), sim.getValue("n_saves"),
            sim.getValue("n_restores")), sep="\\t"), append=T);
    }}
}}
'''
    return changed + diagnostics


def restore_loader_input(retained, original):
    """The unchanged stdpopsim loader still reads its original temporary path."""
    if Path(retained).resolve() == Path(original).resolve():
        raise ValueError('Retained and original tree paths must be distinct')
    shutil.copyfile(retained, original)


def cgroup_paths():
    paths = []
    for line in Path('/proc/self/cgroup').read_text().splitlines():
        _, controllers, relative = line.split(':', 2)
        if 'memory' in controllers.split(','):
            current = Path('/sys/fs/cgroup/memory') / relative.lstrip('/')
            while current != Path('/sys/fs/cgroup/memory'):
                paths.append(current)
                if current.name.startswith('job_'):
                    break
                current = current.parent
        elif controllers == '':
            paths.append(Path('/sys/fs/cgroup') / relative.lstrip('/'))
    return paths


def process_rows():
    pending, rows, seen = [os.getpid()], [], set()
    while pending:
        pid = pending.pop()
        if pid in seen:
            continue
        seen.add(pid)
        try:
            status = Path(f'/proc/{pid}/status').read_text().splitlines()
            row = {'pid': pid}
            row.update({line.split(':', 1)[0]: line.split(':', 1)[1].strip()
                        for line in status if line.startswith(('Name:', 'VmRSS:', 'VmHWM:', 'PPid:'))})
            rows.append(row)
            pending.extend(int(value) for value in Path(f'/proc/{pid}/task/{pid}/children').read_text().split())
        except (FileNotFoundError, ProcessLookupError):
            pass
    return rows


def resource_row(paths, stage):
    groups = {}
    for path in paths:
        values = {}
        for name in ('memory.usage_in_bytes', 'memory.max_usage_in_bytes', 'memory.limit_in_bytes',
                     'memory.failcnt', 'memory.oom_control', 'memory.current', 'memory.peak',
                     'memory.max', 'memory.events'):
            try:
                values[name] = (path / name).read_text().strip()
            except FileNotFoundError:
                pass
        groups[str(path)] = values
    if not any(('memory.usage_in_bytes' in values and 'memory.limit_in_bytes' in values)
               or ('memory.current' in values and 'memory.max' in values) for values in groups.values()):
        raise ValueError('No readable memory usage/limit counters in the allocated cgroup')
    return {'unix_time': time.time(), 'stage': stage['value'], 'cgroups': groups, 'processes': process_rows()}


def monitor_resources(path, paths, stage, stop, errors):
    try:
        with Path(path).open('w', buffering=1) as handle:
            while True:
                handle.write(json.dumps(resource_row(paths, stage)) + '\n')
                if stop.wait(10):
                    handle.write(json.dumps(resource_row(paths, stage)) + '\n')
                    break
    except Exception as error:
        errors.append(f'{type(error).__name__}: {error}')
        try:
            Path(path).with_name('monitor_failure.json').write_text(json.dumps({'errors': errors}) + '\n')
        except OSError:
            pass  # Errors are also retained in adapter_status by the main thread.


def finalize_status(status, status_path, provenance_error, monitor_errors):
    """Persist evidence failures and make the process fail after cleanup."""
    status['frozen_inputs_unchanged'] = provenance_error is None
    status['resource_monitor_complete'] = not monitor_errors
    if provenance_error is not None:
        status.update(state='FAILED', provenance_error=str(provenance_error))
    if monitor_errors:
        status.update(state='FAILED', resource_monitor_errors=list(monitor_errors))
    Path(status_path).write_text(json.dumps(status, indent=2) + '\n')
    print(json.dumps(status))
    if provenance_error is not None:
        raise RuntimeError('Final frozen-input verification failed; diagnostics retained') from provenance_error
    if monitor_errors:
        raise RuntimeError('Resource monitor failed; complete feasibility evidence unavailable')


def main():
    protocol_path = HERE / 'protocol.json'
    protocol_hash = sha(protocol_path)
    protocol = json.loads(protocol_path.read_text())
    locks = {HERE / name: digest for name, digest in protocol['source_sha256'].items()}
    locks.update({Path(name): digest for name, digest in protocol['external_sha256'].items()})
    locks[protocol_path] = protocol_hash

    def verify():
        for path, digest in locks.items():
            if sha(path) != digest:
                raise ValueError(f'Frozen input/source changed: {path}')

    verify()
    snapshot = HERE / 'snapshot'
    sys.path[:0] = [str(snapshot / 'src'), str(snapshot)]
    os.chdir(snapshot)
    import stdpopsim
    versions = {name: importlib.metadata.version(name) for name in protocol['versions']}
    if versions != protocol['versions']:
        raise ValueError('Runtime versions differ from the frozen protocol')
    for name in ('ghist_single_sweep.scan', 'ghist_single_sweep.export', 'ghist_single_sweep.predict',
                 'ghist_submission.validate', 'ghist_relatedness.labels'):
        source = Path(importlib.import_module(name).__file__).resolve()
        if not source.is_relative_to(snapshot):
            raise ValueError(f'Import escaped frozen source snapshot: {name}: {source}')
    output = Path(protocol['output_directory'])
    output.mkdir(exist_ok=False)
    diagnostics = output / 'diagnostics'
    diagnostics.mkdir()
    phase = diagnostics / 'phase.tsv'
    phase.write_text('tick\tdiploid_population_size\tsaves\trestores\n')
    (diagnostics / 'cgroup_membership.txt').write_text(Path('/proc/self/cgroup').read_text())
    stage = {'value': 'initializing'}
    stop = threading.Event()
    paths = cgroup_paths()
    if not paths:
        raise ValueError('No memory cgroup identified; refusing unmonitored feasibility attempt')
    resource_row(paths, stage)  # Fail before launching SLiM if counters are unavailable.
    monitor_errors = []
    thread = threading.Thread(target=monitor_resources,
                              args=(diagnostics / 'resources.jsonl', paths, stage, stop, monitor_errors), daemon=True)
    thread.start()
    logger = logging.getLogger('stdpopsim.slim_engine')
    logger.setLevel(logging.DEBUG)
    handler = logging.FileHandler(diagnostics / 'slim_engine.log')
    logger.addHandler(handler)
    engine = stdpopsim.get_engine('slim')
    original_run = engine._run_slim
    invocation = []

    def run_slim(script_file, **kwargs):
        if invocation:
            raise ValueError('This protocol permits only one SLiM invocation')
        invocation.append(dict(kwargs))
        if kwargs.get('seed') != protocol['seed'] or kwargs.get('slim_path') != protocol['slim_path']:
            raise ValueError('Engine invocation changed its frozen seed or binary')
        original = Path(script_file).read_text()
        original_tree_paths = TREES.findall(original)
        if len(original_tree_paths) != 1:
            raise ValueError('Generated tree destination is missing or ambiguous')
        original_tree_path = Path(original_tree_paths[0])
        (diagnostics / 'generated_original.slim').write_text(original)
        modified = transform_script(original, diagnostics / 'raw_slim.trees', phase)
        (diagnostics / 'generated_interval100.slim').write_text(modified)
        Path(script_file).write_text(modified)
        (diagnostics / 'generated_script_hashes.json').write_text(json.dumps({
            'original_sha256': sha(diagnostics / 'generated_original.slim'),
            'interval100_sha256': sha(diagnostics / 'generated_interval100.slim'),
            'adapter_sha256': sha(Path(__file__)), 'protocol_sha256': protocol_hash,
            'seed': kwargs['seed'], 'biological_parameters_modified': False,
        }, indent=2) + '\n')
        verify()
        stage['value'] = 'slim_running'
        try:
            value = original_run(script_file, **kwargs)
            stage['value'] = 'copying_final_raw_tree_for_original_loader'
            restore_loader_input(diagnostics / 'raw_slim.trees', original_tree_path)
            stage['value'] = 'engine_postprocessing'
            return value
        finally:
            if stage['value'] != 'engine_postprocessing':
                stage['value'] = 'slim_returned_or_failed'
            handler.flush()

    engine._run_slim = run_slim
    started = time.monotonic()
    status = {'state': 'RUNNING', 'protocol_sha256': protocol_hash, 'case': 5, 'seed': protocol['seed'],
              'scientifically_validated': False, 'scope': 'One same-seed resource feasibility attempt only'}
    status_path = diagnostics / 'adapter_status.json'
    status_path.write_text(json.dumps(status, indent=2) + '\n')
    sys.argv = ['scripts.run_sweep_readiness_case', '--design', protocol['design_relative_path'], '--case', '5',
                '--accessible', protocol['accessibility_path'], '--slim', protocol['slim_path'],
                '--raisd', protocol['raisd_path'], '--out-dir', str(output)]
    try:
        runpy.run_module('scripts.run_sweep_readiness_case', run_name='__main__')
        verify()
        status['state'] = 'COMPLETED'
        stage['value'] = 'workflow_complete'
    except BaseException as error:
        status.update(state='FAILED', error=f'{type(error).__name__}: {error}')
        stage['value'] = 'workflow_failed'
        raise
    finally:
        status['elapsed_seconds'] = time.monotonic() - started
        status['slim_invocations'] = len(invocation)
        provenance_error = None
        try:
            verify()
        except Exception as error:
            provenance_error = error
        stop.set()
        thread.join(timeout=2)
        if thread.is_alive():
            monitor_errors.append('Resource monitor did not stop during cleanup')
        engine._run_slim = original_run
        logger.removeHandler(handler)
        handler.close()
        finalize_status(status, status_path, provenance_error, monitor_errors)


if __name__ == '__main__':
    main()
