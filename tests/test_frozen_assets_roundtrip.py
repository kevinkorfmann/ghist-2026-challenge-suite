"""The release bundle must contain file bytes even when a source is a link."""
import hashlib
import importlib.util
import json
from pathlib import Path
import sys
import tarfile


def load_script(name):
    path = Path(__file__).resolve().parents[1] / 'scripts' / f'{name}.py'
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_symlink_first_object_roundtrip(tmp_path, monkeypatch):
    source = tmp_path / 'source'
    source.mkdir()
    library = source / 'library.so.1'
    library.write_bytes(b'example library bytes')
    library.chmod(0o755)
    link = source / 'library.so'
    link.symlink_to(library.name)
    model = source / 'model.bin'
    model.write_bytes(b'example model bytes')
    selection = tmp_path / 'selection.json'
    selection.write_text(json.dumps({
        'original_root': str(source),
        'files': [str(link), str(library), str(model)],
        'expected_sha256': {str(link): hashlib.sha256(library.read_bytes()).hexdigest()},
        'purpose': 'Archive round-trip regression test',
    }))
    archive = tmp_path / 'assets.tar.gz'
    monkeypatch.setattr(sys, 'argv', ['create_frozen_assets.py', str(selection), str(archive)])
    load_script('create_frozen_assets').main()
    with tarfile.open(archive) as bundle:
        members = bundle.getmembers()
        assert len(members) == 3
        assert all(member.isfile() for member in members)
    restored = tmp_path / 'restored'
    restore = load_script('restore_frozen_assets').restore
    restore(archive, restored, verify_only=True)
    assert not restored.exists()
    restore(archive, restored)
    for path in [link, library, model]:
        output = restored / path.name
        assert output.read_bytes() == path.read_bytes()
        assert not output.is_symlink()
    assert (restored / library.name).stat().st_mode & 0o777 == 0o755
