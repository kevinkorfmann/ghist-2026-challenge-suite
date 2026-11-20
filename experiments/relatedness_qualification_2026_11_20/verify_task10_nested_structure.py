"""Check proposed truth independently from Ped-sim's actual parentage graph."""
import hashlib
import json
from collections import Counter
from datetime import UTC, datetime
from itertools import combinations
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT / "experiments/task10_development_2026_11_20"
OUT = ROOT / "experiments/relatedness_qualification_2026_11_20/task10_nested_structure_audit.json"
PROTOCOL_SHA = "8455e205bbd2c257225a1c0bd07c365b78f0e6260211d0e7ba7626986c94e18a"
LABELS = {
    (1, 1): "siblings", (1, 2): "avuncular", (1, 3): "avuncular_once_removed",
    (1, 4): "avuncular_twice_removed", (2, 2): "full_1st_cousin",
    (2, 3): "full_1st_cousin_once_removed", (2, 4): "full_1st_cousin_twice_removed",
    (2, 5): "full_1st_cousin_thrice_removed", (3, 3): "full_2nd_cousin",
    (3, 4): "full_2nd_cousin_once_removed", (3, 5): "full_2nd_cousin_twice_removed",
    (4, 4): "full_3rd_cousin", (4, 5): "full_3rd_cousin_once_removed", (5, 5): "full_4th_cousin",
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load(path):
    return json.loads(path.read_text())


def main():
    manifest_path = BASE / "nested_structure_manifest.json"
    manifest = load(manifest_path)
    expected_hashes = manifest["files_sha256"]
    before = {name: digest(BASE / name) for name in expected_hashes}
    assert before == expected_hashes
    assert digest(BASE / "nested_qualification_protocol.json") == PROTOCOL_SHA
    protocol = load(BASE / "nested_qualification_protocol.json")
    parents, family = {}, {}
    pedigree = BASE / "nested_structure_evidence/nested-everyone.fam"
    for line in pedigree.read_text().splitlines():
        fields = line.split()
        assert len(fields) == 6
        fid, sample, father, mother, *_ = fields
        assert sample not in parents
        assert (father == "0") == (mother == "0")
        parents[sample] = () if father == "0" else (father, mother)
        family[sample] = fid
    for sample, pair in parents.items():
        assert all(parent in parents and family[parent] == family[sample] for parent in pair)
    couples = {frozenset(pair) for pair in parents.values() if pair}

    def ancestors(sample):
        distances, pending = {sample: 0}, [sample]
        while pending:
            descendant = pending.pop()
            for parent in parents[descendant]:
                # Repeated ancestry would violate this design's outbred domain.
                assert parent not in distances, (sample, parent)
                distances[parent] = distances[descendant] + 1
                pending.append(parent)
        return distances

    results = []
    for cohort in protocol["cohorts"]:
        truth = load(BASE / f"nested_structure_evidence/{cohort['name']}_truth.json")
        roster = truth["roster"]
        assert roster == cohort["declared_roster"]
        assert len(roster) == len(set(roster)) == cohort["people"]
        assert {family[sample] for sample in roster} == {cohort["name"] + "1"}
        ancestral = {sample: ancestors(sample) for sample in roster}
        rows = truth["truth"]
        claimed = {tuple(sorted((row["indID1"], row["indID2"]))): row["relation"] for row in rows}
        assert len(claimed) == len(rows) == len(roster) * (len(roster) - 1) // 2 == cohort["pairs"]
        expected_pairs = set(combinations(sorted(roster), 2))
        assert set(claimed) == expected_pairs
        counts = Counter()
        for a, b in expected_pairs:
            aa, bb = ancestral[a], ancestral[b]
            common = aa.keys() & bb.keys()
            assert common
            distance = min(aa[node] + bb[node] for node in common)
            closest = {node for node in common if aa[node] + bb[node] == distance}
            assert len(closest) == 2 and frozenset(closest) in couples
            paths = {(aa[node], bb[node]) for node in closest}
            assert len(paths) == 1
            label = LABELS[tuple(sorted(next(iter(paths))))]
            assert label == claimed[(a, b)], (a, b, label, claimed[(a, b)])
            counts[label] += 1
        assert set(counts) == set(LABELS.values())
        assert dict(counts) == cohort["label_pair_counts"]
        results.append({"name": cohort["name"], "people": len(roster), "pairs": len(rows),
                        "independent_label_counts": dict(sorted(counts.items()))})
    assert sum(row["pairs"] for row in results) == manifest["total_pairs"] == 53764
    assert before == {name: digest(BASE / name) for name in expected_hashes}
    report = {"verified_at_utc": datetime.now(UTC).isoformat(), "protocol_sha256": PROTOCOL_SHA,
              "manifest_sha256": digest(manifest_path), "inputs_sha256_before_and_after": before,
              "audit_source_sha256": digest(Path(__file__)), "cohorts": results,
              "pair_labels_verified": 53764, "scope": "Structure only; no genotypes or accuracy qualification",
              "method": "Walk actual parent edges, reject repeated ancestry, require one closest ancestral couple, derive labels from the two path lengths.",
              "release_authorized": False}
    assert not OUT.exists()
    OUT.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps({"cohorts": len(results), "pairs_verified": 53764, "all_fourteen_labels_each": True}))


if __name__ == "__main__":
    main()
