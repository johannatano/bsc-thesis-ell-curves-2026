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
    p.add_argument("--use-hcp", action="store_true", default=False, help="Use HCP (Hilbert class polynomial) enumeration instead of direct method")
    p.add_argument("--use-cn", action="store_true", default=False, help="Use Class Numbers ie no j invariants instead of direct method")
    p.add_argument("--rank-method", choices=["auto", "div_poly", "mod_poly", "invariants"], default="auto", help="Method for above-floor rank detection (default: auto — div_poly for ℓ<13, mod_poly otherwise)")
    p.add_argument("--true-height", action="store_true", default=False, help="Use exact BFS height in isogeny volcano instead of floor test")
    return p.parse_args()

def run(p: int, n: int, l: int, use_HCP: bool = False, use_CN: bool = False):
    primes = list(primerange(5, 50)) if p == -1 else [p]
    p_powers = [i for i in range(1, 3)] if n == -1 else [n]
    levels = [l] if l != -1 else list(primerange(2, 50))
    dsize = len(primes)
    q_max = 10**20
    for i in range(dsize):
        p = primes[i]
        nf = None
        if use_HCP or use_CN:
            NFC = NumberFieldsClassifier_Fq(p)
            nf = NFC.generate(p_powers, q_max=q_max)
        for n in p_powers:
            q = p**n
            if q > q_max:
                print(f"Skipping F_{q} due to size > {q_max}")
                continue
            CC = CurvesClassifier_Fq(p, n, NF=nf)
            CC.enumerate_curves(use_HCP=use_HCP, use_CN=use_CN)
            
            for ell in levels:
                #CC.compute_volcano(ell=ell, edges=False)
                N_EP = CC.count_EP(ell=ell, use_CN=use_CN)
                print(f"{Colors.GREEN}F_{q}: Total count of (E,P) for ell={ell}: {N_EP}{Colors.ENDC}")
    
            
if __name__ == "__main__":
    args = parse_args()
    Config.rank_method = args.rank_method
    Config.use_true_height = args.true_height
    print("\n")
    print("="*80 + "")
    print(f"Using rank detection method: {Config.rank_method}, use_true_conductor: {Config.use_true_height}, use_HCP: {args.use_hcp}, use_CN: {args.use_cn}")
    print("="*80 + "\n")
    start_hcp = time.time()
    run(args.p, args.n, args.l, use_HCP=args.use_hcp, use_CN=args.use_cn)
    end_hcp = time.time()
    print(f"{Colors.HEADER}Classification completed in {end_hcp - start_hcp:.2f} seconds{Colors.ENDC}")
    
    
        
        
    
    
    