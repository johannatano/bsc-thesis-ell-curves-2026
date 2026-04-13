'''D_K_l = kronecker(K.discriminant(), l)
            D_K = K.discriminant()
            #if D_K > -5:
            #    continue  # skip small discriminants
            #print(K)
            #print(K.maximal_order())
            H = pari(4*self.q - t*t).qfbhclassno()
            #print(H)
            #f_pi_factors = f_pi.conductor.factor()'''


'''Dpi = ZZ(t*t - 4*q)
            K = QuadraticField(Dpi)          # same field as Q(sqrt(Dpi))
            DK = ZZ(K.discriminant())        # fundamental discriminant
            fpi2 = ZZ(Dpi // DK)
            fpi = fpi2.isqrt()
            assert fpi*fpi == fpi2'''
'''# add for c = 1
            #n_tot = O_K.class_number()
            #O_pi = K.order_of_conductor(f_pi)
            #if(f_pi != 1):
            #    n_tot += O_pi.class_number()
            print(f"D_pi={t**2-4*self.q}, f_pi: {f_pi}, D_K={D_K}, D_K_factors={D_K.factor()}")
            
            divs = sorted(f_pi.divisors())
            n_tot = 0
            for f in divs:
                Of = K.order_of_conductor(f)               # Of = Z + f*O_K
                h_O = Of.class_number()
                n_tot += h_O
                level = f_pi // f
                l_volcano_floor = level % l != 0
                r = 1 if l_volcano_floor else 2
                N_EP_ += h_O*(l**r-1)//2
                print(f, Of.discriminant(), Of.class_number(), l_volcano_floor)'''


import argparse
import sys
import os
import time
from sympy import primerange
from sage.all import *
import requests

def parse_args():
    p = argparse.ArgumentParser(description="Classify curves over F_q.")
    p.add_argument("-p", "--p", type=int, required=False, default=-1, help="Field char p")
    p.add_argument("-n", "--n", type=int, required=False, default=1, help="Field extension degree n")
    p.add_argument("-l", "--l", type=int, required=False, default=-1, help="Level ℓ")
    p.add_argument("-k", "--k", type=int, required=False, default=2, help="Weight k")
    p.add_argument("--use-hcp", action="store_true", default=False, help="Use HCP (Hilbert class polynomial) enumeration instead of direct method")
    p.add_argument("--use-cn", action="store_true", default=False, help="Use Class Numbers ie no j invariants instead of direct method")
    p.add_argument(
        "--rank-method",
        choices=["auto", "div_poly", "mod_poly", "invariants"],
        default="mod_poly",
        help="Method for above-floor rank detection (default: auto — div_poly for ℓ<13, mod_poly otherwise)",
    )
    p.add_argument("--true-height", action="store_true", default=False, help="Use exact BFS height in isogeny volcano instead of floor test")
    return p.parse_args()


def compute_nf(D_pi):
    # Step 1: Extract square-free part (d_sf)
    d_sf = D_pi
    i = 2
    while i * i <= abs(d_sf):
        while d_sf % (i * i) == 0:
            d_sf //= i * i
        i += 1

    # Step 2: Identify fundamental discriminant D_K
    if d_sf % 4 == 1:
        D_K = d_sf
    else:
        D_K = 4 * d_sf

    f_pi = int((D_pi / D_K) ** 0.5)

    K = QuadraticField(D_pi)
    DK = ZZ(K.discriminant())
    fpi2 = ZZ(D_pi // DK)
    fpi = fpi2.isqrt()
    assert fpi == f_pi
    assert DK == D_K

    return D_K, f_pi, DK, fpi


def enum_ell(p: int, n:int, ell: int):
    q = p**n
    HB = math.isqrt(4 * q)
    N_EP = 0
    i_min = (q + 1 - HB + ell - 1) // ell
    i_max = (q + 1 + HB) // ell
    for i in range(i_min, i_max + 1):
        t = q + 1 - i*ell
        if t % p != 0:
            N_EP += count_t(q, t, ell)
    return N_EP


def count_t(q:int, t:int, ell:int):
    D_pi = t**2-4*q
    K = QuadraticField(D_pi)
    D_K = K.discriminant()
    f_pi = int((D_pi / D_K) ** 0.5)
    S = 0
    if D_K < -4:
        divs = ZZ(f_pi).divisors()
        for f in divs:
            Of = K.order_of_conductor(f)               # Of = Z + f*O_K
            h_O = 1#Of.class_number()
            D = f**2*D_K
            level = f_pi // f
            l_volcano_floor = level % ell != 0
            r = 1 if l_volcano_floor else 2
            S += h_O * (ell**r - 1) // 2
    return S

def count_t_multi(q:int, t:int, ell:list):
    D_pi = t**2-4*q
    K = QuadraticField(D_pi)
    D_K = K.discriminant()
    f_pi = int((D_pi / D_K) ** 0.5)
    S = [0]*len(ell)
    if D_K < -4:
        divs = ZZ(f_pi).divisors()
        for f in divs:
            Of = K.order_of_conductor(f)               # Of = Z + f*O_K
            h_O = 1#Of.class_number()
            D = f**2*D_K
            level = f_pi // f
            for i in range(0, len(ell)):
                l = ell[i]
                l_volcano_floor = level % l != 0
                r = 1 if l_volcano_floor else 2
                S[i] += h_O * (l**r - 1) // 2
    return S

def run(p: int, l:int, k:int, n:int, use_HCP=False, use_CN=False):
    primes = list(primerange(5, 20)) if p == -1 else [p]
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
        q = p**n
        if q > q_max:
            print(f"Skipping F_{q} due to size > {q_max}")
            continue
        if len(levels) == 1:
            ell = levels[0]
            N_EP = enum_ell(p, n, ell)
            print(f"q={q}, ell={ell}, #EP: {N_EP}")
            continue
        HB = math.isqrt(4*q)
        for t in range(1, HB):
            if t % p != 0:
                N_EP_LIST = count_t_multi(q, t, levels)
                print(f"q={q}, t={t}, #EP-list: {N_EP_LIST}")
    print("\n" + "="*80)
    print("="*80)

if __name__ == "__main__":
    args = parse_args()
    print("\n")
    print("="*80 + "")
    print("="*80 + "\n")
    start_hcp = time.time()
    run(args.p, args.l, args.k, args.n, use_HCP=args.use_hcp, use_CN=args.use_cn)
    end_hcp = time.time()
    print(f"Hecke Trace computed in {end_hcp - start_hcp:.2f} seconds")
