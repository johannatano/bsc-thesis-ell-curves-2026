import argparse
import sys
import os
import time

from utils.ell_curves import CurvesClassifier_Fq
from utils.ell_nr_field import NumberFieldsClassifier_Fq
from utils.common import Logger, Colors, Data, Config
from math import gcd
from sympy import primerange
from sage.all import *
import numpy as np
import requests
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
from sage.all import DirichletGroup

def parse_args():
    p = argparse.ArgumentParser(description="Classify curves over F_q.")
    p.add_argument("-p", "--p", type=int, required=False, default=-1, help="Field char p")
    p.add_argument("-n", "--n", type=int, required=False, default=-1, help="Field extension degree n")
    p.add_argument("-l", "--l", type=int, required=False, default=-1, help="Level ℓ")
    p.add_argument("-k", "--k", type=int, required=False, default=2, help="Weight k")
    p.add_argument("--use-hpc", action="store_true", default=False, help="Use HPC (Hilbert class polynomial) enumeration instead of direct method")
    p.add_argument("--rank-method", choices=["auto", "div_poly", "mod_poly", "invariants"], default="auto", help="Method for above-floor rank detection (default: auto — div_poly for ℓ<13, mod_poly otherwise)")
    p.add_argument("--true-height", action="store_true", default=False, help="Use exact BFS height in isogeny volcano instead of floor test")
    p.add_argument("--max-ell", type=int, default=50, help="Maximum torsion prime ℓ for compute_torsion (default: 50)")
    return p.parse_args()

def enum(p: int, n: int, use_HPC=False, max_ell=50):
    primes = list(primerange(5, 50)) if p == -1 else [p]
    p_powers = [i for i in range(1, 3)] if n == -1 else [n]
    dsize = len(primes)
    q_max = 10**20  # Set a maximum q to avoid long computations
    for i in range(dsize):
        p = primes[i]
        nf = None
        if use_HPC:
            data = Data.loadJSON(f"./data/{p}/nr_fields.json")
            if data is not None:
                NFC = NumberFieldsClassifier_Fq.fromJson(data)  # Load the number fields data into the classifier
                nf = NFC.nr_fields  # Access the number fields data
            else:
                NFC = NumberFieldsClassifier_Fq(p)
                nf = NFC.generate(p_powers, q_max=q_max)
        for n in p_powers:
            q = p**n
            if q > q_max:
                print(f"Skipping F_{q} due to size > 100000")
                continue
            CC = CurvesClassifier_Fq(p, n, NF=nf)
            CC.enumerate_curves(use_HPC=use_HPC) #here we create curves by j invariant, and generate the twists
            CC.compute_torsion(max_ell=max_ell) #
            
if __name__ == "__main__":
    args = parse_args()
    Config.rank_method = args.rank_method
    Config.use_true_height = args.true_height
    print("\n")
    print("="*80 + "")
    print(f"Using rank detection method: {Config.rank_method}, use_true_height: {Config.use_true_height}, max_ell: {args.max_ell}")
    print("="*80 + "\n")
    start_hcp = time.time()
    enum(args.p, args.n, use_HPC=args.use_hpc, max_ell=args.max_ell)
    end_hcp = time.time()
    print(f"{Colors.HEADER}Classification completed in {end_hcp - start_hcp:.2f} seconds{Colors.ENDC}")
    
    
        
        
    
    
    