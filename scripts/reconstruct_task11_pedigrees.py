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


BONSAI_AUTOSOMAL_CM = 3544.57


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd", type=Path, required=True)
    parser.add_argument("--map", type=Path, required=True)
    parser.add_argument("--clusters", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--genome-cm", type=float, required=True)
    parser.add_argument("--threshold-cm", type=float, default=10.0)
    args = parser.parse_args()

    from bonsaitree import bonsai

    cluster_summary = json.loads(args.clusters.read_text())
    bonsai.validate_bio = lambda bio, idx: None
    family_results = []
    relationship_counts: dict[str, int] = defaultdict(int)
    all_relationships = []
    for family_index, members in enumerate(cluster_summary["families"]):
        segments, bio_info, id_to_sample = convert_segments(
            args.ibd,
            args.map,
            args.threshold_cm,
            allowed_samples=set(members),
            cm_scale=BONSAI_AUTOSOMAL_CM / args.genome_cm,
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
        family_relationships = []
        for offset, identifier1 in enumerate(identifiers):
            for identifier2 in identifiers[offset + 1 :]:
                relationship = pedigree.rels.get(identifier1, {}).get(identifier2)
                if relationship is None:
                    continue
                relationship = tuple(relationship)
                relationship_counts[str(relationship)] += 1
                row = {
                    "family": family_index,
                    "indID1": id_to_sample[identifier1],
                    "indID2": id_to_sample[identifier2],
                    "relationship_tuple": relationship,
                }
                family_relationships.append(row)
                all_relationships.append(row)
        family_results.append(
            {
                "family": family_index,
                "members": members,
                "input_segments": len(segments),
                "placed_pairs": len(family_relationships),
                "topology": serializable_topology(pedigree.up_pedigree_dict),
            }
        )

    summary = {
        "state": "COMPLETED",
        "families": family_results,
        "relationship_counts": dict(sorted(relationship_counts.items())),
        "relationships": all_relationships,
        "length_scale": BONSAI_AUTOSOMAL_CM / args.genome_cm,
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
