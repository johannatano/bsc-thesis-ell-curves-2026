import argparse
import sys
import os
import time

from lib.curves_classifier import CurvesClassifier_Fq
from lib.nr_fields_classifier import NumberFieldsClassifier_Fq

from utils.common import Logger, Colors, Data, Config
from sympy import primerange
from sage.all import *
import numpy as np
import requests

def parse_args():
    p = argparse.ArgumentParser(description="Classify curves over F_q.")
    p.add_argument("-p", "--p", type=int, required=False, default=-1, help="Field char p")
    p.add_argument("-n", "--n", type=int, required=False, default=1, help="Field extension degree n")
    p.add_argument("-N", "--N", type=int, required=False, default=-1, help="Level ℓ")
    p.add_argument("-k", "--k", type=int, required=False, default=2, help="Weight k")
    p.add_argument("--use-hcp", action="store_true", default=False, help="Use HCP (Hilbert class polynomial) enumeration instead of direct method")
    p.add_argument("--use-cn", action="store_true", default=False, help="Use Class Numbers ie no j invariants instead of direct method")
    p.add_argument(
        "--rank-method",
        choices=["auto", "div_poly", "mod_poly", "invariants"],
        default="mod_poly",
        help="Method for above-floor rank detection (default: auto — div_poly for ℓ<13, mod_poly otherwise)",
    )
    p.add_argument(
        "--compare",
        action="store_true",
        default=False,
        help="Use Sage's built-in Hecke operator for comparison (only for small p and N)",
    )
    p.add_argument("--true-height", action="store_true", default=False, help="Use exact BFS height in isogeny volcano instead of floor test")
    return p.parse_args()

def run(p: int, l:int, k:int, n:int, use_HCP=False, use_CN=False):
    primes = list(primerange(5, 100)) if p == -1 else [p]
    # primes = list(primerange(10**6, 10**6+100)) if p == -1 else [p]
    # 1000033
    # 10093
    # p=1091
    # 100043
    # 1000003
    # 10050013, 10050017, 10050023, 10050049, 10050059, 10050071, 10050083, 10050101, 10050133, 10050137, 10050167, 10050181, 10050191, 10050197, 10050203, 10050217, 10050223, 10050233, 10050253, 10050283, 10050317, 10050319, 10050331, 10050353, 10050367, 10050377, 10050389, 10050407, 10050413, 10050419, 10050427, 10050437, 10050463, 10050493
    p_powers = [n]
    dsize = len(primes)
    q_max = 10**20
    levels = [l] if l != -1 else list(primerange(2, 15))
    diffs = {}  # (ell, k) -> list of diffs per prime p
    for i in range(dsize):
        p = primes[i]
        nf = None
        if use_HCP or use_CN:
            NFC = NumberFieldsClassifier_Fq(p)
            nf = NFC.generate(p_powers, q_max=q_max)
        q = p**n
        if q > q_max:
            print(f"Skipping F_{q} due to size > {q_max}")
            continue
        CC = CurvesClassifier_Fq(p, n, NF=nf)
        CC.enumerate_curves(use_HCP=use_HCP, use_CN=use_CN, add_curves=True, add_SS=True)

        for ell in levels:
            T, NC, NSS, traces, hk_evals, vals, full_r = CC.compute_hecke(k=args.k, level=ell, use_CN=use_CN)
            trace_val = CuspForms(Gamma1(ell), k + 2).hecke_operator(q).trace()
            diff = T - trace_val
            print(sorted(traces))
            diffs.setdefault((ell, k, (ell - 1) // 2), []).append(
                (q, diff, T, trace_val, NC, NSS, vals, traces, hk_evals)
            )
    print("\n" + "="*80)
    print("Diff summary (ell, k, dim) -> [(p, diff)]")
    print("="*80)
    for key, entries in sorted(diffs.items()):
        ell, k_, dim = key
        print(f"  (ell={ell}, k={k_}):")
        for (
            p_val,
            d,
            T_val,
            trace_val,
            NC_val,
            NSS_val,
            vals,
            traces_val,
            hk_evals_val,
        ) in entries:
            if p_val % ell == 1:
                color = Colors.GREEN
                label = " [q≡1 mod ell]"
            elif p_val % ell == ell - 1:
                color = Colors.WARNING if k_ % 2 == 1 else Colors.GREEN
                label = " [q≡-1 mod ell]"
            else:
                color = Colors.ENDC
                label = ""
            print(
                f"{color}q={p_val}: diff={d}, T={T_val}, sage_trace={trace_val}, NC={NC_val}, NSS={NSS_val}, {label} [q≡{p_val % ell} mod {ell}]{Colors.ENDC}"
             )
if __name__ == "__main__":
    args = parse_args()
    Config.rank_method = args.rank_method
    print("\n")
    print("="*80 + "")
    print(f"Using rank detection method: {Config.rank_method}")
    print("="*80 + "\n")
    start_hcp = time.time()
    run(args.p, args.N, args.k, args.n, use_HCP=args.use_hcp, use_CN=args.use_cn)
    end_hcp = time.time()
    print(f"{Colors.HEADER}Hecke Trace computed in {end_hcp - start_hcp:.2f} seconds{Colors.ENDC}")
