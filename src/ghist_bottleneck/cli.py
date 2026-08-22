from __future__ import annotations

import argparse
import json
from pathlib import Path

import yaml

from .features import FeatureBundle, extract_vcf
from .infer import predict
from .simulate import generate_dataset
from .train import train_ensemble
from .sfs_baseline import fit_two_epoch, write_fit


def main() -> None:
    p = argparse.ArgumentParser(prog="ghist-bottleneck")
    sub = p.add_subparsers(dest="command", required=True)
    ex = sub.add_parser("extract")
    ex.add_argument("--vcf", required=True)
    ex.add_argument("--out", required=True)
    ex.add_argument("--trees", help="Optional inferred tree sequence")
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
    tr.add_argument("--config", help="YAML model and training configuration")
    inf = sub.add_parser("infer")
    inf.add_argument("--observed", required=True)
    inf.add_argument("--checkpoint-dir", required=True)
    inf.add_argument("--out", required=True)
    base = sub.add_parser("sfs-reference")
    base.add_argument("--vcf", required=True)
    base.add_argument("--out", required=True)
    base.add_argument("--restarts", type=int, default=24)
    args = p.parse_args()
    if args.command == "extract":
        b = extract_vcf(args.vcf, tree_sequence=args.trees)
        b.save(args.out)
        print(json.dumps({"n_haplotypes": b.n_haplotypes, "n_variants": b.n_variants, "sequence_length": b.sequence_length}, indent=2))
    elif args.command == "simulate":
        generate_dataset(
            FeatureBundle.load(args.observed), args.out,
            n_sims=args.n_sims, seed=args.seed, quick=args.quick,
        )
    elif args.command == "train":
        configuration = yaml.safe_load(Path(args.config).read_text()) if args.config else {}
        training = configuration.get("training", {})
        train_ensemble(
            args.dataset,
            args.out_dir,
            quick=args.quick,
            model_config=configuration.get("model", {}),
            seeds=tuple(training.get("seeds", (17, 29, 43, 71, 101))),
            epochs=int(training.get("epochs", 80)),
            learning_rate=float(training.get("learning_rate", 2e-4)),
            weight_decay=float(training.get("weight_decay", 2e-3)),
            batch_size=int(training.get("batch_size", 64)),
        )
    elif args.command == "infer":
        print(json.dumps(predict(args.observed, args.checkpoint_dir, args.out), indent=2))
    else:
        result = fit_two_epoch(args.vcf, restarts=args.restarts)
        write_fit(result, args.out)
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
