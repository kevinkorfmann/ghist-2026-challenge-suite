from __future__ import annotations

import argparse
import bisect
import csv
import gzip
import heapq
from collections import defaultdict
from pathlib import Path


BIN_EDGES = (0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 60.0, 80.0, 120.0, float("inf"))
ORDERED = 20


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--ibd", type=Path, required=True)
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--out", type=Path, required=True)
    args = parser.parse_args()

    histograms: dict[tuple[str, str], list[int]] = defaultdict(
        lambda: [0] * (len(BIN_EDGES) - 1)
    )
    longest: dict[tuple[str, str], list[float]] = defaultdict(list)
    with gzip.open(args.ibd, "rt") as handle:
        for line in handle:
            left, _, right, _, _, _, _, length_text = line.split()
            pair = tuple(sorted((left, right)))
            length = float(length_text)
            bin_index = bisect.bisect_right(BIN_EDGES, length) - 1
            histograms[pair][min(bin_index, len(BIN_EDGES) - 2)] += 1
            heap = longest[pair]
            if len(heap) < ORDERED:
                heapq.heappush(heap, length)
            elif length > heap[0]:
                heapq.heapreplace(heap, length)

    with args.base.open() as handle:
        base_rows = list(csv.DictReader(handle))
    rows = []
    for row in base_rows:
        pair = tuple(sorted((row["indID1"], row["indID2"])))
        output = dict(row)
        output.update(
            {
                f"length_bin_{index}": count
                for index, count in enumerate(histograms[pair])
            }
        )
        ordered = sorted(longest[pair], reverse=True)
        output.update(
            {
                f"ordered_length_{rank + 1}": ordered[rank]
                if rank < len(ordered)
                else 0.0
                for rank in range(ORDERED)
            }
        )
        rows.append(output)
    args.out.parent.mkdir(parents=True, exist_ok=True)
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    print(f"pairs={len(rows)} output={args.out}")


if __name__ == "__main__":
    main()
