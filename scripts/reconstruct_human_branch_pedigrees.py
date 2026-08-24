from __future__ import annotations

import argparse
import json
from collections import defaultdict
from pathlib import Path

from reconstruct_bonsai_pedigree import (
    convert_segments,
    segment_summaries,
    serializable_topology,
)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd", type=Path, required=True)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--threshold-cm", type=float, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    from bonsaitree import bonsai

    clusters = json.loads(args.clusters.read_text())
    bonsai.validate_bio = lambda bio, idx: None
    relationship_counts: dict[str, int] = defaultdict(int)
    relationships = []
    branch_results = []
    for branch_index, members in enumerate(clusters["families"]):
        if len(members) < 2:
            branch_results.append(
                {"branch": branch_index, "members": members, "relationships": 0}
            )
            continue
        segments, bio_info, id_to_sample = convert_segments(
            args.ibd,
            args.map,
            args.threshold_cm,
            allowed_samples=set(members),
        )
        result = bonsai.build_pedigree(
            segments,
            bio_info,
            ibd_stats=segment_summaries(segments, sorted(id_to_sample)),
            max_add_degree=4,
            max_rel_append_types=3,
            ibd_threshold=0,
            remove_distant_threshold=float("inf"),
        )
        pedigree = result["ped_obj"]
        identifiers = sorted(id_to_sample)
        branch_rows = []
        for offset, identifier1 in enumerate(identifiers):
            for identifier2 in identifiers[offset + 1 :]:
                relationship = pedigree.rels.get(identifier1, {}).get(identifier2)
                if relationship is None:
                    continue
                relationship = tuple(relationship)
                relationship_counts[str(relationship)] += 1
                row = {
                    "branch": branch_index,
                    "indID1": id_to_sample[identifier1],
                    "indID2": id_to_sample[identifier2],
                    "relationship_tuple": relationship,
                }
                branch_rows.append(row)
                relationships.append(row)
        branch_results.append(
            {
                "branch": branch_index,
                "members": members,
                "input_segments": len(segments),
                "relationships": len(branch_rows),
                "topology": serializable_topology(pedigree.up_pedigree_dict),
            }
        )
    summary = {
        "state": "COMPLETED",
        "individuals": clusters["samples"],
        "branches": branch_results,
        "relationship_counts": dict(sorted(relationship_counts.items())),
        "relationships": relationships,
        "metadata_policy": "age and sex unknown",
        "submission_count": 0,
    }
    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(summary, indent=2) + "\n")
    print(
        json.dumps(
            {key: value for key, value in summary.items() if key != "relationships"},
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
