from __future__ import annotations

import argparse
import json

from .features import FeatureBundle, extract_vcf
from .infer import predict
from .simulate import generate_dataset
from .train import train_ensemble


def main() -> None:
    p = argparse.ArgumentParser(prog="ghist-bottleneck")
    sub = p.add_subparsers(dest="command", required=True)
    ex = sub.add_parser("extract")
    ex.add_argument("--vcf", required=True)
    ex.add_argument("--out", required=True)
    sim = sub.add_parser("simulate")
    sim.add_argument("--observed", required=True)
    sim.add_argument("--out", required=True)
    sim.add_argument("--n-sims", type=int, default=5000)
    sim.add_argument("--seed", type=int, default=20260821)
    sim.add_argument("--quick", action="store_true")
    tr = sub.add_parser("train")
    tr.add_argument("--dataset", required=True)
    tr.add_argument("--out-dir", required=True)
    tr.add_argument("--quick", action="store_true")
    inf = sub.add_parser("infer")
    inf.add_argument("--observed", required=True)
    inf.add_argument("--checkpoint-dir", required=True)
    inf.add_argument("--out", required=True)
    args = p.parse_args()
    if args.command == "extract":
        b = extract_vcf(args.vcf)
        b.save(args.out)
        print(json.dumps({"n_haplotypes": b.n_haplotypes, "n_variants": b.n_variants, "sequence_length": b.sequence_length}, indent=2))
    elif args.command == "simulate":
        generate_dataset(
            FeatureBundle.load(args.observed), args.out,
            n_sims=args.n_sims, seed=args.seed, quick=args.quick,
        )
    elif args.command == "train":
        train_ensemble(args.dataset, args.out_dir, quick=args.quick)
    else:
        print(json.dumps(predict(args.observed, args.checkpoint_dir, args.out), indent=2))


if __name__ == "__main__":
    main()
