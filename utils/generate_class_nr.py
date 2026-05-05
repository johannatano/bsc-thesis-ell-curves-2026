#!/usr/bin/env python3
"""Generate data/class_numbers/ bucket JSONs: h(O_K) for all |D| <= D_MAX.

Saves one JSON per 10k discriminants:
    class_numbers_00001_10000.json, class_numbers_10001_20000.json, ...

Requires SageMath. Run as:
    sage utils/generate_class_nr.py [--max-abs-D 10000]
"""
import json
import argparse
from pathlib import Path

from sage.all import *

_BUCKET = 10_000
_DEFAULT_OUT_DIR = Path(__file__).parent / "data" / "class_numbers"


def _bucket_path(out_dir: Path, idx: int) -> Path:
    lo = idx * _BUCKET + 1
    hi = (idx + 1) * _BUCKET
    return out_dir / f"class_numbers_{lo:05d}_{hi:05d}.json"


def generate(max_abs_D: int, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    current: dict[str, int] = {}
    cur_bucket = 0

    for abs_D in range(1, max_abs_D + 1):
        bucket = (abs_D - 1) // _BUCKET
        if bucket != cur_bucket:
            with open(_bucket_path(out_dir, cur_bucket), "w") as fh:
                json.dump(current, fh)
            print(f"Saved {_bucket_path(out_dir, cur_bucket).name} ({len(current)} entries)")
            current = {}
            cur_bucket = bucket

        D = -abs_D
        K = QuadraticField(D, "a")
        current[str(D)] = int(K.class_number())

        if abs_D % 1000 == 0:
            print(f"  D={D} ...")

    if current:
        with open(_bucket_path(out_dir, cur_bucket), "w") as fh:
            json.dump(current, fh)
        print(f"Saved {_bucket_path(out_dir, cur_bucket).name} ({len(current)} entries)")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--max-abs-D", type=int, default=10_000)
    ap.add_argument("--out-dir", default=str(_DEFAULT_OUT_DIR))
    args = ap.parse_args()

    print(f"Computing h(O_K) for |D| <= {args.max_abs_D} ...")
    generate(args.max_abs_D, Path(args.out_dir))
