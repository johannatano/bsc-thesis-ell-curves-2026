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
    p.add_argument("-n", "--n", type=int, required=False, default=-1, help="Field extension degree n")
    p.add_argument("-l", "--l", type=int, required=False, default=-1, help="Level ℓ")
    p.add_argument("-k", "--k", type=int, required=False, default=2, help="Weight k")
    p.add_argument("--use-hcp", action="store_true", default=False, help="Use HCP (Hilbert class polynomial) enumeration instead of direct method")
    p.add_argument("--use-cn", action="store_true", default=False, help="Use Class Numbers ie no j invariants instead of direct method")
    p.add_argument("--rank-method", choices=["auto", "div_poly", "mod_poly", "invariants"], default="auto", help="Method for above-floor rank detection (default: auto — div_poly for ℓ<13, mod_poly otherwise)")
    p.add_argument("--true-height", action="store_true", default=False, help="Use exact BFS height in isogeny volcano instead of floor test")
    return p.parse_args()
            
def run(p: int, l:int, k:int, use_HCP=False, use_CN=False):
    primes = list(primerange(5, 50)) if p == -1 else [p]
    p_powers = [1]
    dsize = len(primes)
    q_max = 10**20
    levels = [l] if l != -1 else list(primerange(2, 20))
    for ell in levels:
        print(f"\n{Colors.BLUE}=== Computing Hecke operator T_{ell} for weight {k} ==={Colors.ENDC}\n")
        for i in range(dsize):
            p = primes[i]
            nf = None
            if use_HCP or use_CN:
                NFC = NumberFieldsClassifier_Fq(p)
                nf = NFC.generate(p_powers, q_max=q_max)
            q = p
            if q > q_max:
                print(f"Skipping F_{q} due to size > {q_max}")
                continue
            CC = CurvesClassifier_Fq(p, 1, NF=nf)
            CC.enumerate_curves(use_HCP=use_HCP, use_CN=use_CN)
            T = CC.compute_hecke(k=args.k, level=ell, use_CN=use_CN)
            print(f"Total Hecke trace for level {ell} and weight {args.k}: {T}")
            
if __name__ == "__main__":
    args = parse_args()
    Config.rank_method = args.rank_method
    print("\n")
    print("="*80 + "")
    print(f"Using rank detection method: {Config.rank_method}")
    print("="*80 + "\n")
    start_hcp = time.time()
    run(args.p, args.l, args.k, use_HCP=args.use_hcp, use_CN=args.use_cn)
    end_hcp = time.time()
    print(f"{Colors.HEADER}Hecke Trace computed in {end_hcp - start_hcp:.2f} seconds{Colors.ENDC}")
    
    
        
        
    
    
    