import gzip
import importlib.util
from pathlib import Path

import pytest

spec = importlib.util.spec_from_file_location("map_stage_test", Path(__file__).with_name("stage_task11_map_bounds_v1.py"))
stage = importlib.util.module_from_spec(spec)
spec.loader.exec_module(stage)
HEADER = b"##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ta\tb\tc\td\te\tf\tg\th\n"


def row(chromosome, position):
    return f"{chromosome}\t{position}\t.\tA\tG\t.\tPASS\t.\tGT\t0|1\t1|0\t.|.\t0|0\t1|1\t0|1\t1|0\t0|0\n".encode()


def inputs(tmp_path, data):
    map_path, vcf = tmp_path / "map", tmp_path / "input.vcf"
    map_path.write_text("1 . 0 10\n1 . 1 100\n2 . 0 1\n2 . 1 20\n")
    vcf.write_bytes(HEADER + data)
    return map_path, vcf


def test_inclusive_bounds_preserve_retained_genotypes_and_header_bytes(tmp_path):
    kept = row(1, 10) + row(1, 100) + row(2, 1) + row(2, 20)
    map_path, vcf = inputs(tmp_path, row(1, 9) + row(1, 10) + row(1, 100) + row(1, 101) + row(2, 1) + row(2, 20) + row(2, 21))
    before = stage.sha(vcf)
    target = tmp_path / "cropped.gz"
    audit = stage.stage(vcf, map_path, target, stage.sha(map_path), 1)
    assert gzip.decompress(target.read_bytes()) == HEADER + kept
    assert audit["excluded_records"] == 3
    assert audit["chromosomes"]["1"]["excluded_below"] == 1
    assert audit["chromosomes"]["1"]["excluded_above"] == 1
    assert audit["chromosomes"]["2"]["retained_raw_markers"] == 2
    assert stage.sha(vcf) == before


def test_identity_and_deterministic_gzip(tmp_path):
    map_path, vcf = inputs(tmp_path, row(1, 10) + row(1, 100) + row(2, 1) + row(2, 20))
    one, two = tmp_path / "one.gz", tmp_path / "two.gz"
    a = stage.stage(vcf, map_path, one, stage.sha(map_path), 1)
    b = stage.stage(one, map_path, two, stage.sha(map_path), 1)
    assert one.read_bytes() == two.read_bytes()
    assert gzip.decompress(one.read_bytes()) == vcf.read_bytes()
    assert a["excluded_records"] == b["excluded_records"] == 0
    assert a["input_variant_lines_sha256"] == a["retained_variant_lines_sha256"]


@pytest.mark.parametrize("data,minimum", [
    (row(1, 10) + row(3, 1), 1),
    (row(1, 10), 1),
    (row(1, 10) + row(2, 1), 2),
    (row(1, 100) + row(1, 10) + row(2, 1), 1)])
def test_reject_unmapped_missing_underpopulated_or_unsorted_chromosomes(tmp_path, data, minimum):
    map_path, vcf = inputs(tmp_path, data)
    with pytest.raises(ValueError):
        stage.stage(vcf, map_path, tmp_path / "out.gz", stage.sha(map_path), minimum)


def test_reject_changed_map_and_existing_output(tmp_path):
    map_path, vcf = inputs(tmp_path, row(1, 10) + row(2, 1))
    with pytest.raises(ValueError):
        stage.stage(vcf, map_path, tmp_path / "out.gz", "wrong", 1)
    target = tmp_path / "existing.gz"
    target.write_bytes(b"preserve")
    with pytest.raises(FileExistsError):
        stage.stage(vcf, map_path, target, stage.sha(map_path), 1)
    assert target.read_bytes() == b"preserve"
