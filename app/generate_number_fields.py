import argparse
import sys
import os
from lib.nr_fields_classifier import NumberFieldsClassifier_Fq
from utils.common import Logger, Colors, Data
from math import gcd
from sympy import primerange
from sage.all import *
import numpy as np
import requests
from sage.schemes.elliptic_curves.ell_finite_field import supersingular_j_polynomial
from sage.all import DirichletGroup

def parse_args():
    p = argparse.ArgumentParser(description="Classify curves over F_q.")
    p.add_argument("-n", "--n", type=int, required=False, default=-1, help="Field extension degree n")
    p.add_argument("-p", "--p", type=int, required=False, default=-1, help="Field char p")
    p.add_argument("-q_max", "--q_max", type=int, required=False, default=1000000, help="Field extension degree n")
    return p.parse_args()

def run(p: int, n:int, q_max:int):
    primes = list(primerange(5, 100)) if p == -1 else [p]
    p_powers = [i for i in range(1, 10)] if n == -1 else [n]
    if n != -1:
        q_max = p ** n + 1
    for i in range(len(primes)): #skip char=2,3
        p = primes[i]
        CNF = NumberFieldsClassifier_Fq(p)
        CNF.generate(p_powers, q_max=q_max)
        Data.saveJSON(f"data/{p}", f"nr_fields.json", CNF.toJSON(), readable=False) 
     
if __name__ == "__main__":
    args = parse_args()
    print(f"Generating number fields for p={args.p}, n={args.n}, q_max={args.q_max}...")
    run(args.p, args.n, args.q_max)