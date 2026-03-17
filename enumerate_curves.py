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


catalogue_by_Q = {}
catalogue_by_HCP = {}


def enum(p: int, n: int, use_HPC=False, max_ell=50):
    primes = list(primerange(5, 50)) if p == -1 else [p]
    p_powers = [i for i in range(1, 3)] if n == -1 else [n]
    dsize = len(primes)
    q_max = 10**20  # Set a maximum q to avoid long computations
    for i in range(dsize):
        p = primes[i]
        nf = None
        
        # we prepopulate a catalogue of possible traces, hence discriminants ie nr fields for this char, up to a max n
        if use_HPC:
            print(f"\n=== Processing prime p={p}, loading json.... ===")
            data = Data.loadJSON(f"./data/{p}/nr_fields.json")
            NFC = NumberFieldsClassifier_Fq.fromJson(data)  # Load the number fields data into the classifier
            nf = NFC.nr_fields  # Access the number fields data
            
            #NFC = NumberFieldsClassifier_Fq(p)
            #nf = NFC.generate(p_powers, q_max=q_max)

        #data = Data.loadJSON(f"./data/{p}/nr_fields.json")
        #nrFieldsData = data['nr_fields']['nf']  # Extract the list of number fields
        for n in p_powers:
            q = p**n
            if q > q_max:
                print(f"Skipping F_{q} due to size > 100000")
                continue
            # if we do not load precomputed json, first generate the number fields catalogue for this p and n
            # create the classifier for curves over F_q, add the precomputed nr fields data, and enumerate curves using HCP
            CC = CurvesClassifier_Fq(p, n, NF=nf)
            CC.enumerate_curves(use_HPC=use_HPC) #here we create curves by j invariant, and generate the twists
            CC.compute_torsion(max_ell=max_ell) #
            #CC.enumerate_curves()
            #print(f"\nTotal unique j-invariants found: {len(added_js)}")
            catalogue_by_HCP[q] = CC
            #CC.compute_volcanoes()
            #print(f"Saving results for F_{q}...")
            #Data.saveJSON(f"./data/{CFq.p}", f"curves_HCP_{q}.json", CFq.toJSON(), readable=False)
            #print(f"Done with F_{q}. {CC.catalogue.size} curves found.\n") 
            

def debug():
    for q in catalogue_by_Q:
        dataQ = catalogue_by_Q[q]
        dataHCP = catalogue_by_HCP[q]
        print(f"\n=== Comparing results for F_{q} ===")
        # Get j-invariants from both methods
        # catalogue has number_fields -> tree -> isogeny_classes -> orders -> curves
        j_invs_Q = set()
        j_invs_HCP = set()
        
        for nf_info in dataQ.catalogue.NFC.data.values():
            for tree in nf_info.tree:
                for ic in tree.isogeny_classes.values():
                    for f, curves in ic.curves_by_order.items():
                        for curve in curves:
                            j_invs_Q.add(curve.j)
                            print(f"Q method: Found curve with j={curve.j} and trace={curve.t}")
        
        for nf_info in dataHCP.catalogue.NFC.data.values():
            for tree in nf_info.tree:
                for ic in tree.isogeny_classes.values():
                    for f, curves in ic.curves_by_order.items():
                        for curve in curves:
                            j_invs_HCP.add(curve.j)
                            print(f"HCP method: Found curve with j={curve.j} and trace={curve.t}")
        
        print(f"J-invariants found by Q method: {len(j_invs_Q)}")
        print(f"J-invariants found by HCP method: {len(j_invs_HCP)}")
        
        # Check for differences
        only_in_Q = j_invs_Q - j_invs_HCP
        only_in_HCP = j_invs_HCP - j_invs_Q
        common = j_invs_Q & j_invs_HCP
        
        if only_in_Q:
            print(f"  ⚠ {len(only_in_Q)} j-invariants ONLY in Q method: {list(only_in_Q)[:5]}{'...' if len(only_in_Q) > 5 else ''}")
        if only_in_HCP:
            print(f"  ⚠ {len(only_in_HCP)} j-invariants ONLY in HCP method: {list(only_in_HCP)[:5]}{'...' if len(only_in_HCP) > 5 else ''}")
        if not only_in_Q and not only_in_HCP:
            print(f"  ✓ Perfect match! All {len(common)} j-invariants agree")
        else:
            print(f"  {len(common)} j-invariants in common")
        
        # Compare total curve counts
        total_curves_Q = sum(len(curves) for nf_info in dataQ.catalogue.NFC.data.values() 
                            for tree in nf_info.tree for ic in tree.isogeny_classes.values() 
                            for f, curves in ic.curves_by_order.items())
        total_curves_HCP = sum(len(curves) for nf_info in dataHCP.catalogue.NFC.data.values() 
                              for tree in nf_info.tree for ic in tree.isogeny_classes.values() 
                              for f, curves in ic.curves_by_order.items())
        
        print(f"\n  Total curves by Q method: {total_curves_Q}")
        print(f"  Total curves by HCP method: {total_curves_HCP}")
        if total_curves_Q == total_curves_HCP:
            print(f"  ✓ Curve counts match!")
        else:
            print(f"  ⚠ Curve count mismatch!")
            
if __name__ == "__main__":
    args = parse_args()
    Config.rank_method = args.rank_method
    Config.use_true_height = args.true_height
    print("\n")
    print("="*80 + "")
    print(f"Using rank detection method: {Config.rank_method}, use_true_height: {Config.use_true_height}, max_ell: {args.max_ell}")
    print("="*80 + "\n")
    #print("\n" + "="*60)
    #print("Starting enumeration...")
    start_hcp = time.time()
    enum(args.p, args.n, use_HPC=args.use_hpc, max_ell=args.max_ell)
    end_hcp = time.time()
    print(f"{Colors.HEADER}Classification completed in {end_hcp - start_hcp:.2f} seconds{Colors.ENDC}")
    #print("="*60 + "\n")
    
    '''print("\n" + "="*60)
    print("Starting Q enumeration...")
    start_q = time.time()
    #enum(args.p, args.n, use_HPC=False)
    end_q = time.time()
    print(f"Q enumeration completed in {end_q - start_q:.2f} seconds")
    print("="*60 + "\n")
    
    print(f"\nTotal execution time: {end_q - start_hcp:.2f} seconds")
    
    print("\n" + "="*60)
    print("COMPARING RESULTS")
    print("="*60)
    #debug()'''
    
    
        
        
    
    
    